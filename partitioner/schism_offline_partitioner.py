from __future__ import annotations

import argparse
import json
import pickle
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
from tqdm import tqdm

from parse_schism_inputs import SchismGrid, load_schism_grid

try:
    import pymetis
except ImportError:  # pragma: no cover - optional backend
    pymetis = None


@dataclass(slots=True)
class PartitionGraph:
    vwgt: np.ndarray
    xadj: np.ndarray
    adjncy: np.ndarray
    adjwgt: np.ndarray
    centers: np.ndarray


def parse_layout_spec(spec: str) -> list[int]:
    layout: list[int] = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "x" in chunk.lower():
            value_text, repeat_text = chunk.lower().split("x", maxsplit=1)
            layout.extend([int(value_text)] * int(repeat_text))
        else:
            layout.append(int(chunk))
    if not layout:
        raise ValueError("Layout spec is empty")
    return layout


def load_layout(args: argparse.Namespace) -> list[int]:
    if args.layout_json is not None:
        with Path(args.layout_json).open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if isinstance(payload, dict) and "nodes" in payload:
            return [int(item["compute_ranks"]) for item in payload["nodes"]]
        if isinstance(payload, list):
            return [int(value) for value in payload]
        raise ValueError("layout JSON must be a list of compute ranks or an object with a nodes array")

    if args.layout_spec is not None:
        return parse_layout_spec(args.layout_spec)

    return parse_layout_spec("28x28,16")


def _graph_cache_path(grid: SchismGrid, cache_dir: Path) -> Path:
    stem = Path(grid.hgrid_path).stem
    vstem = Path(grid.vgrid_path).stem
    return cache_dir / f"graph_{stem}_{vstem}.pkl"


def build_node_incidence(grid: SchismGrid) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    counts = np.zeros(grid.nnode, dtype=np.int32)
    for ie in tqdm(range(grid.nelem), desc="Counting node incidence", unit="elem"):
        etype = int(grid.element_types[ie])
        for local_node in range(etype):
            counts[grid.elements[ie, local_node]] += 1

    offsets = np.empty(grid.nnode + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])

    cursor = offsets[:-1].copy()
    node_elements = np.empty(int(offsets[-1]), dtype=np.int32)
    for ie in tqdm(range(grid.nelem), desc="Building node incidence", unit="elem"):
        etype = int(grid.element_types[ie])
        for local_node in range(etype):
            node = grid.elements[ie, local_node]
            slot = cursor[node]
            node_elements[slot] = ie
            cursor[node] += 1

    return counts, offsets, node_elements


def build_side_data(grid: SchismGrid) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    total_sides = int(np.sum(grid.element_types, dtype=np.int64))
    edge_keys = np.empty(total_sides, dtype=np.int64)
    edge_owner = np.empty(total_sides, dtype=np.int32)

    cursor = 0
    stride = np.int64(grid.nnode)
    for ie in tqdm(range(grid.nelem), desc="Encoding sides", unit="elem"):
        etype = int(grid.element_types[ie])
        conn = grid.elements[ie]
        for local_side in range(etype):
            n0 = int(conn[local_side])
            n1 = int(conn[(local_side + 1) % etype])
            if n0 > n1:
                n0, n1 = n1, n0
            edge_keys[cursor] = np.int64(n0) * stride + np.int64(n1)
            edge_owner[cursor] = ie
            cursor += 1

    order = np.argsort(edge_keys, kind="mergesort")
    edge_keys = edge_keys[order]
    edge_owner = edge_owner[order]

    boundary_side_count = np.zeros(grid.nelem, dtype=np.int8)
    pair_count = 0
    total_sides = edge_keys.size
    start = 0
    while start < total_sides:
        stop = start + 1
        while stop < total_sides and edge_keys[stop] == edge_keys[start]:
            stop += 1
        group_size = stop - start
        if group_size == 1:
            boundary_side_count[edge_owner[start]] += 1
        else:
            pair_count += group_size * (group_size - 1) // 2
        start = stop

    pair_i = np.empty(pair_count, dtype=np.int32)
    pair_j = np.empty(pair_count, dtype=np.int32)
    cursor = 0
    start = 0
    while start < total_sides:
        stop = start + 1
        while stop < total_sides and edge_keys[stop] == edge_keys[start]:
            stop += 1
        owners = edge_owner[start:stop]
        if owners.size > 1:
            for lhs, rhs in combinations(owners.tolist(), 2):
                pair_i[cursor] = lhs
                pair_j[cursor] = rhs
                cursor += 1
        start = stop

    counts = np.zeros(grid.nelem, dtype=np.int32)
    np.add.at(counts, pair_i, 1)
    np.add.at(counts, pair_j, 1)

    offsets = np.empty(grid.nelem + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])

    cursor = offsets[:-1].copy()
    neighbors = np.empty(int(offsets[-1]), dtype=np.int32)
    for lhs, rhs in zip(pair_i, pair_j, strict=True):
        slot = cursor[lhs]
        neighbors[slot] = rhs
        cursor[lhs] += 1
        slot = cursor[rhs]
        neighbors[slot] = lhs
        cursor[rhs] += 1

    return boundary_side_count, offsets, neighbors


def compute_vertex_weights(
    grid: SchismGrid, nne: np.ndarray, boundary_side_count: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nelem = grid.nelem
    vwgt = np.empty((nelem, 4), dtype=np.int32)
    side_cost = np.empty(nelem, dtype=np.int32)
    node_cost = np.empty(nelem, dtype=np.int32)

    for ie in tqdm(range(nelem), desc="Computing SCHISM-like weights", unit="elem"):
        etype = int(grid.element_types[ie])
        conn = grid.elements[ie, :etype]
        kbetmp = int(np.min(grid.kbp[conn]))
        nlev = max(1, grid.nvrt - kbetmp)

        ptmp = float(np.sum(1.0 / nne[conn], dtype=np.float64))
        n_boundary = int(boundary_side_count[ie])
        stmp = n_boundary + 0.5 * (etype - n_boundary)

        vwgt[ie, 0] = int(np.rint(float(nlev)))
        vwgt[ie, 1] = int(np.rint(ptmp * nlev))
        vwgt[ie, 2] = int(np.rint(stmp * nlev))
        vwgt[ie, 3] = 1

        side_cost[ie] = int(nlev * (2 * etype - 3))
        node_cost[ie] = int(nlev * (2 * etype - 1))

    return vwgt, side_cost, node_cost


def build_dual_graph(
    grid: SchismGrid,
    node_offsets: np.ndarray,
    node_elements: np.ndarray,
    side_offsets: np.ndarray,
    side_neighbors: np.ndarray,
    side_cost: np.ndarray,
    node_cost: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nelem = grid.nelem
    seen = np.zeros(nelem, dtype=np.int32)
    counts = np.zeros(nelem, dtype=np.int32)

    for ie in tqdm(range(nelem), desc="Counting dual-graph edges", unit="elem"):
        marker = ie + 1
        etype = int(grid.element_types[ie])

        for idx in range(int(side_offsets[ie]), int(side_offsets[ie + 1])):
            je = int(side_neighbors[idx])
            if seen[je] != marker:
                seen[je] = marker
                counts[ie] += 1

        for local_node in range(etype):
            node = int(grid.elements[ie, local_node])
            for idx in range(int(node_offsets[node]), int(node_offsets[node + 1])):
                je = int(node_elements[idx])
                if je == ie:
                    continue
                if seen[je] != marker:
                    seen[je] = marker
                    counts[ie] += 1

    xadj = np.empty(nelem + 1, dtype=np.int64)
    xadj[0] = 0
    np.cumsum(counts, out=xadj[1:])

    adjncy = np.empty(int(xadj[-1]), dtype=np.int32)
    adjwgt = np.empty(int(xadj[-1]), dtype=np.int32)

    seen.fill(0)
    for ie in tqdm(range(nelem), desc="Building dual graph", unit="elem"):
        marker = ie + 1
        etype = int(grid.element_types[ie])
        cursor = int(xadj[ie])

        for idx in range(int(side_offsets[ie]), int(side_offsets[ie + 1])):
            je = int(side_neighbors[idx])
            if seen[je] == marker:
                continue
            seen[je] = marker
            adjncy[cursor] = je
            adjwgt[cursor] = int(side_cost[ie] + side_cost[je])
            cursor += 1

        for local_node in range(etype):
            node = int(grid.elements[ie, local_node])
            for idx in range(int(node_offsets[node]), int(node_offsets[node + 1])):
                je = int(node_elements[idx])
                if je == ie or seen[je] == marker:
                    continue
                seen[je] = marker
                adjncy[cursor] = je
                adjwgt[cursor] = int(node_cost[ie] + node_cost[je])
                cursor += 1

        if cursor != int(xadj[ie + 1]):
            raise RuntimeError(f"Dual graph fill mismatch for element {ie}: {cursor} != {xadj[ie + 1]}")

    return xadj, adjncy, adjwgt


def build_partition_graph(
    grid: SchismGrid,
    cache_dir: Path | None = None,
    force_rebuild: bool = False,
) -> PartitionGraph:
    cache_path = None
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = _graph_cache_path(grid, cache_dir)

    if cache_path is not None and cache_path.exists() and not force_rebuild:
        with cache_path.open("rb") as handle:
            return pickle.load(handle)

    nne, node_offsets, node_elements = build_node_incidence(grid)
    boundary_side_count, side_offsets, side_neighbors = build_side_data(grid)
    vwgt, side_cost, node_cost = compute_vertex_weights(grid, nne, boundary_side_count)
    xadj, adjncy, adjwgt = build_dual_graph(
        grid,
        node_offsets,
        node_elements,
        side_offsets,
        side_neighbors,
        side_cost,
        node_cost,
    )

    centers = np.empty((grid.nelem, 2), dtype=np.float64)
    for ie in range(grid.nelem):
        etype = int(grid.element_types[ie])
        conn = grid.elements[ie, :etype]
        centers[ie] = grid.nodes_xy[conn].mean(axis=0)

    graph = PartitionGraph(
        vwgt=vwgt,
        xadj=xadj,
        adjncy=adjncy,
        adjwgt=adjwgt,
        centers=centers,
    )

    if cache_path is not None:
        with cache_path.open("wb") as handle:
            pickle.dump(graph, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return graph


def call_pymetis(
    xadj: np.ndarray,
    adjncy: np.ndarray,
    adjwgt: np.ndarray,
    vwgt: np.ndarray,
    nparts: int,
    tpwgts: np.ndarray | None = None,
) -> np.ndarray:
    if pymetis is None:
        raise RuntimeError("pymetis is not installed in this environment")
    if nparts == 1:
        return np.zeros(vwgt.shape[0], dtype=np.int32)

    kwargs: dict[str, object] = {
        "xadj": xadj.tolist(),
        "adjncy": adjncy.tolist(),
        "eweights": adjwgt.tolist(),
        "vweights": vwgt.reshape(-1).tolist(),
        "nparts": nparts,
    }
    if tpwgts is not None:
        kwargs["tpwgts"] = tpwgts.astype(np.float32).ravel().tolist()

    result = pymetis.part_graph(**kwargs)
    if isinstance(result, tuple):
        _, part = result
    else:
        part = result
    return np.asarray(part, dtype=np.int32)


def recursive_coordinate_bisection(
    centers: np.ndarray,
    scalar_weight: np.ndarray,
    part_weights: np.ndarray,
) -> np.ndarray:
    if centers.shape[0] != scalar_weight.size:
        raise ValueError("centers and scalar_weight must have the same number of vertices")

    part_assignment = np.empty(centers.shape[0], dtype=np.int32)
    part_ids = np.arange(part_weights.size, dtype=np.int32)

    def recurse(vertex_ids: np.ndarray, current_parts: np.ndarray) -> None:
        if current_parts.size == 1:
            part_assignment[vertex_ids] = int(current_parts[0])
            return

        target = np.cumsum(part_weights[current_parts])
        split_pos = int(np.searchsorted(target, target[-1] / 2.0, side="left")) + 1
        split_pos = min(max(split_pos, 1), current_parts.size - 1)
        left_parts = current_parts[:split_pos]
        right_parts = current_parts[split_pos:]

        local_centers = centers[vertex_ids]
        spreads = np.ptp(local_centers, axis=0)
        axis = int(np.argmax(spreads))
        order = np.argsort(local_centers[:, axis], kind="mergesort")
        sorted_vertices = vertex_ids[order]
        sorted_weights = scalar_weight[sorted_vertices]
        cumulative = np.cumsum(sorted_weights, dtype=np.float64)
        target_left = cumulative[-1] * (part_weights[left_parts].sum() / part_weights[current_parts].sum())
        cut = int(np.searchsorted(cumulative, target_left, side="left")) + 1
        cut = min(max(cut, 1), sorted_vertices.size - 1)

        recurse(sorted_vertices[:cut], left_parts)
        recurse(sorted_vertices[cut:], right_parts)

    recurse(np.arange(centers.shape[0], dtype=np.int32), part_ids)
    return part_assignment


def partition_with_backend(
    backend: str,
    centers: np.ndarray,
    xadj: np.ndarray,
    adjncy: np.ndarray,
    adjwgt: np.ndarray,
    vwgt: np.ndarray,
    nparts: int,
    tpwgts: np.ndarray | None = None,
) -> np.ndarray:
    if nparts == 1:
        return np.zeros(vwgt.shape[0], dtype=np.int32)

    chosen = backend
    if chosen == "auto":
        chosen = "pymetis" if pymetis is not None else "rcb"

    if chosen == "pymetis":
        return call_pymetis(xadj, adjncy, adjwgt, vwgt, nparts, tpwgts)
    if chosen == "rcb":
        scalar_weight = np.maximum(1, vwgt.sum(axis=1).astype(np.int64))
        if tpwgts is None:
            part_weights = np.ones(nparts, dtype=np.float64)
        else:
            reshaped = tpwgts.reshape(nparts, -1)
            part_weights = reshaped[:, 0]
        return recursive_coordinate_bisection(centers, scalar_weight, part_weights)

    raise ValueError(f"Unknown backend: {backend}")


def induced_subgraph(
    xadj: np.ndarray,
    adjncy: np.ndarray,
    adjwgt: np.ndarray,
    vertices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mapping = np.full(xadj.size - 1, -1, dtype=np.int32)
    mapping[vertices] = np.arange(vertices.size, dtype=np.int32)

    counts = np.zeros(vertices.size, dtype=np.int32)
    for local_idx, global_idx in enumerate(vertices):
        for idx in range(int(xadj[global_idx]), int(xadj[global_idx + 1])):
            mapped = mapping[adjncy[idx]]
            if mapped >= 0:
                counts[local_idx] += 1

    local_xadj = np.empty(vertices.size + 1, dtype=np.int64)
    local_xadj[0] = 0
    np.cumsum(counts, out=local_xadj[1:])

    local_adjncy = np.empty(int(local_xadj[-1]), dtype=np.int32)
    local_adjwgt = np.empty(int(local_xadj[-1]), dtype=np.int32)
    for local_idx, global_idx in enumerate(vertices):
        cursor = int(local_xadj[local_idx])
        for idx in range(int(xadj[global_idx]), int(xadj[global_idx + 1])):
            mapped = mapping[adjncy[idx]]
            if mapped >= 0:
                local_adjncy[cursor] = mapped
                local_adjwgt[cursor] = adjwgt[idx]
                cursor += 1

    mapping[vertices] = -1
    return local_xadj, local_adjncy, local_adjwgt


def hierarchical_partition(
    graph: PartitionGraph, layout: list[int], backend: str
) -> tuple[np.ndarray, np.ndarray]:
    coarse_weights = np.asarray(layout, dtype=np.float64)
    coarse_weights /= coarse_weights.sum()

    coarse_part = partition_with_backend(
        backend,
        graph.centers,
        graph.xadj,
        graph.adjncy,
        graph.adjwgt,
        graph.vwgt,
        len(layout),
        coarse_weights,
    )

    final_part = np.empty(graph.vwgt.shape[0], dtype=np.int32)
    rank_offset = 0
    for coarse_id, local_parts in enumerate(layout):
        vertices = np.flatnonzero(coarse_part == coarse_id)
        if vertices.size == 0:
            raise RuntimeError(f"Coarse partition {coarse_id} is empty")
        if local_parts == 1:
            final_part[vertices] = rank_offset
            rank_offset += 1
            continue

        local_xadj, local_adjncy, local_adjwgt = induced_subgraph(
            graph.xadj, graph.adjncy, graph.adjwgt, vertices
        )
        local_vwgt = graph.vwgt[vertices]
        local_centers = graph.centers[vertices]
        local_part = partition_with_backend(
            backend,
            local_centers,
            local_xadj,
            local_adjncy,
            local_adjwgt,
            local_vwgt,
            local_parts,
        )
        final_part[vertices] = local_part + rank_offset
        rank_offset += local_parts

    return coarse_part, final_part


def summarize_partition(
    coarse_part: np.ndarray,
    final_part: np.ndarray,
    xadj: np.ndarray,
    adjncy: np.ndarray,
    adjwgt: np.ndarray,
) -> dict[str, int | float]:
    cut_edges = 0
    weighted_cut = 0
    offnode_cut_edges = 0
    offnode_weighted_cut = 0

    for ie in range(final_part.size):
        for idx in range(int(xadj[ie]), int(xadj[ie + 1])):
            je = int(adjncy[idx])
            if je <= ie:
                continue
            if final_part[ie] != final_part[je]:
                cut_edges += 1
                weighted_cut += int(adjwgt[idx])
                if coarse_part[ie] != coarse_part[je]:
                    offnode_cut_edges += 1
                    offnode_weighted_cut += int(adjwgt[idx])

    unique_parts, counts = np.unique(final_part, return_counts=True)
    return {
        "nparts": int(unique_parts.size),
        "cut_edges": int(cut_edges),
        "weighted_cut": int(weighted_cut),
        "offnode_cut_edges": int(offnode_cut_edges),
        "offnode_weighted_cut": int(offnode_weighted_cut),
        "min_elements_per_rank": int(counts.min()),
        "max_elements_per_rank": int(counts.max()),
        "mean_elements_per_rank": float(counts.mean()),
    }


def write_partition_prop(path: Path, part: np.ndarray) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for ie, rank in enumerate(part, start=1):
            handle.write(f"{ie} {int(rank)}\n")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Offline hierarchical SCHISM partitioner with node-aware layout."
    )
    parser.add_argument("hgrid", type=Path, help="Path to hgrid.gr3")
    parser.add_argument("vgrid", type=Path, help="Path to vgrid.in")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("partition.prop"),
        help="Output partition.prop path",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("partitioner/.cache"),
        help="Directory for parsed/grid caches",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=Path("partition_summary.json"),
        help="Where to write partition summary JSON",
    )
    parser.add_argument(
        "--layout-spec",
        default="28x28,16",
        help="Node compute layout, e.g. '28x28,16' for 28 nodes with 28 ranks and one node with 16",
    )
    parser.add_argument(
        "--layout-json",
        type=Path,
        default=None,
        help="Optional JSON file containing per-node compute ranks",
    )
    parser.add_argument(
        "--force-reparse",
        action="store_true",
        help="Ignore cached parsed hgrid/vgrid data",
    )
    parser.add_argument(
        "--force-rebuild-graph",
        action="store_true",
        help="Ignore cached dual-graph data",
    )
    parser.add_argument(
        "--backend",
        choices=("auto", "pymetis", "rcb"),
        default="auto",
        help="Partition backend. 'auto' prefers pymetis when available, otherwise uses recursive coordinate bisection",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    layout = load_layout(args)
    if sum(layout) <= 0:
        raise ValueError("Layout must contain at least one compute rank")

    grid = load_schism_grid(
        args.hgrid,
        args.vgrid,
        cache_dir=args.cache_dir,
        force_reparse=args.force_reparse,
    )
    graph = build_partition_graph(
        grid,
        cache_dir=args.cache_dir,
        force_rebuild=args.force_rebuild_graph,
    )

    coarse_part, final_part = hierarchical_partition(graph, layout, args.backend)
    summary = summarize_partition(
        coarse_part,
        final_part,
        graph.xadj,
        graph.adjncy,
        graph.adjwgt,
    )
    summary["layout"] = layout

    write_partition_prop(args.output, final_part)
    with args.summary_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
