from __future__ import annotations

import hashlib
import pickle
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import Iterator

import numpy as np
from tqdm import tqdm


@dataclass(slots=True)
class SchismGrid:
    hgrid_path: str
    vgrid_path: str
    nelem: int
    nnode: int
    nvrt: int
    nodes_xy: np.ndarray
    nodes_depth: np.ndarray
    elements: np.ndarray
    element_types: np.ndarray
    kbp: np.ndarray


def _file_signature(path: Path) -> tuple[str, int, int]:
    stat = path.stat()
    return (str(path.resolve()), stat.st_size, int(stat.st_mtime_ns))


def _iter_tokens(path: Path) -> Iterator[str]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            for token in line.split():
                yield token


def _default_cache_path(hgrid_path: Path, vgrid_path: Path, cache_dir: Path) -> Path:
    digest = hashlib.sha256()
    digest.update(str(hgrid_path.resolve()).encode("utf-8"))
    digest.update(str(vgrid_path.resolve()).encode("utf-8"))
    digest.update(str(hgrid_path.stat().st_size).encode("utf-8"))
    digest.update(str(vgrid_path.stat().st_size).encode("utf-8"))
    return cache_dir / f"parsed_{digest.hexdigest()[:16]}.pkl"


def read_hgrid(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    path = Path(path)

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        handle.readline()
        nelem, nnode = map(int, handle.readline().split())

        nodes_xy = np.empty((nnode, 2), dtype=np.float64)
        nodes_depth = np.empty(nnode, dtype=np.float64)

        for _ in tqdm(range(nnode), desc="Reading hgrid nodes", unit="node"):
            parts = handle.readline().split()
            idx = int(parts[0]) - 1
            nodes_xy[idx, 0] = float(parts[1])
            nodes_xy[idx, 1] = float(parts[2])
            nodes_depth[idx] = float(parts[3])

        elements = np.full((nelem, 4), -1, dtype=np.int32)
        element_types = np.empty(nelem, dtype=np.int8)

        for ie in tqdm(range(nelem), desc="Reading hgrid elements", unit="elem"):
            parts = handle.readline().split()
            etype = int(parts[1])
            element_types[ie] = etype
            for local_node in range(etype):
                elements[ie, local_node] = int(parts[2 + local_node]) - 1

    return nodes_xy, nodes_depth, elements, element_types


def read_vgrid_kbp(path: str | Path, nnode: int) -> tuple[int, np.ndarray]:
    path = Path(path)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        ivcor = int(handle.readline().split()[0])
        if ivcor != 1:
            raise ValueError(
                f"Only ivcor=1 is supported for offline partitioning; got ivcor={ivcor}"
            )
        nvrt = int(handle.readline().split()[0])

    tokens = _iter_tokens(path)
    next(tokens)
    next(tokens)
    kbp = np.fromiter(islice(tokens, nnode), dtype=np.int32, count=nnode) - 1
    if kbp.size != nnode:
        raise ValueError(f"Expected {nnode} kbp entries in {path}, found {kbp.size}")

    return nvrt, kbp


def load_schism_grid(
    hgrid_path: str | Path,
    vgrid_path: str | Path,
    cache_dir: str | Path | None = None,
    force_reparse: bool = False,
) -> SchismGrid:
    hgrid_path = Path(hgrid_path)
    vgrid_path = Path(vgrid_path)

    cache_path = None
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = _default_cache_path(hgrid_path, vgrid_path, cache_dir)

    if cache_path is not None and cache_path.exists() and not force_reparse:
        with cache_path.open("rb") as handle:
            payload = pickle.load(handle)
        if payload["hgrid_sig"] == _file_signature(hgrid_path) and payload["vgrid_sig"] == _file_signature(vgrid_path):
            return payload["grid"]

    nodes_xy, nodes_depth, elements, element_types = read_hgrid(hgrid_path)
    nvrt, kbp = read_vgrid_kbp(vgrid_path, nodes_xy.shape[0])

    grid = SchismGrid(
        hgrid_path=str(hgrid_path.resolve()),
        vgrid_path=str(vgrid_path.resolve()),
        nelem=elements.shape[0],
        nnode=nodes_xy.shape[0],
        nvrt=nvrt,
        nodes_xy=nodes_xy,
        nodes_depth=nodes_depth,
        elements=elements,
        element_types=element_types,
        kbp=kbp,
    )

    if cache_path is not None:
        payload = {
            "hgrid_sig": _file_signature(hgrid_path),
            "vgrid_sig": _file_signature(vgrid_path),
            "grid": grid,
        }
        with cache_path.open("wb") as handle:
            pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return grid
