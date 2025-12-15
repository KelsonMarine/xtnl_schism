module scribe_io
   use schism_glbl, only : rkind,errmsg,natrm,max_ncoutvar
   use schism_msgp, only : comm_schism,comm_scribe,nproc_schism,nproc_scribe,nscribes, &
   &myrank_scribe,myrank_schism,rtype,itype,parallel_abort
   use netcdf
   implicit none
   include 'mpif.h'
   private

   !===============================================================================
   ! New Type: OutputVariable
   ! Encapsulates all information needed to handle one output variable
   !===============================================================================
   type :: OutputVariable
      ! Identification
      character(len=20) :: name                ! Variable name (e.g., 'temperature')
      integer :: id                            ! Unique ID for this variable
      integer :: global_index                  ! Position in overall output sequence

      ! Dimensionality and location
      integer :: dimension                     ! 2 or 3 for 2D/3D
      integer :: location                      ! 1=node, 2=elem, 3=side
      integer :: num_components                ! 1=scalar, 2=vector
      integer :: i23d_flag                     ! Flag for output type

      ! Work assignment
      integer :: assigned_scribe               ! Which scribe handles this variable
      integer :: mpi_tag_base                  ! Base tag for MPI communication

      ! Size information
      integer :: horizontal_size               ! np_global, ne_global, or ns_global
      integer :: vertical_size                 ! 1 for 2D, nvrt for 3D

      ! NetCDF information
      integer :: ncid                          ! NetCDF file ID (for 3D vars)
      integer :: varid                         ! NetCDF variable ID
      logical :: file_open                     ! Whether file is currently open

      ! Data storage (pointer to appropriate global array)
      ! These point to shared arrays, not allocated per-variable
      !   real(4), pointer :: data_2d(:,:) => null()      ! (horizontal, 1) for 2D
      !   real(4), pointer :: data_3d(:,:) => null()      ! (vertical, horizontal) for 3D

      ! Module/feature info
      character(len=20) :: module_name         ! e.g., 'hydro', 'wwm', 'sed'
      integer :: module_index                  ! Index within module's output array

   contains
      procedure :: initialize => var_initialize
      procedure :: assign_scribe => var_assign_scribe
      procedure :: receive_data => var_receive_data
      procedure :: write_output => var_write_output
      procedure :: cleanup => var_cleanup
   end type OutputVariable

   !===============================================================================
   ! Module-level variables
   !===============================================================================

   ! Core simulation parameters
!    integer,save :: nvrt,nproc_compute,np_global,ne_global,ns_global
!    integer,save :: np_max,ne_max,ns_max
!    real(rkind), save :: dt,h0
!    integer,save :: ihfskip,nspool,nc_out,ifile,ics,iof_ugrid

   ! Output variable management
   type(OutputVariable), allocatable, save :: output_vars(:)
   integer, save :: num_output_vars
   integer, save :: num_2d_vars, num_3d_vars

   ! Shared data arrays (allocated once, reused by variables)
!    real(4), save, allocatable :: shared_2d_node(:,:)      ! (np_global, 1)
!    real(4), save, allocatable :: shared_2d_elem(:,:)      ! (ne_global, 1)
!    real(4), save, allocatable :: shared_2d_side(:,:)      ! (ns_global, 1)
!    real(4), save, allocatable :: shared_3d_node(:,:)      ! (nvrt, np_global)
!    real(4), save, allocatable :: shared_3d_elem(:,:)      ! (nvrt, ne_global)
!    real(4), save, allocatable :: shared_3d_side(:,:)      ! (nvrt, ns_global)

   ! Temporary receive buffers (per-rank data before assembly)
!    real(4), save, allocatable :: recv_buffer_2d(:,:,:)    ! (max_vals, max_size, nproc)
!    real(4), save, allocatable :: recv_buffer_3d(:,:,:)    ! (nvrt, max_size, nproc)

   ! Mesh data (unchanged from original)
!    integer,save,allocatable :: np(:),ne(:),ns(:),iplg(:,:),ielg(:,:),islg(:,:),kbp00(:)
!    real(rkind),save,allocatable :: xnd(:),ynd(:),dp(:),xel(:),yel(:),xsd(:),ysd(:)
!    integer,save,allocatable :: i34(:),elnode(:,:),isidenode(:,:)

   integer,save :: node_dim,nele_dim,nedge_dim,four_dim,nv_dim, &
   &one_dim,two_dim,time_dim,itime_id,elnode_id, iside_id, i34_id,ix_id,iy_id,ih_id
!    integer, save:: ixel_id2, iyel_id2, ixsd_id2, iysd_id2
   integer,save :: node_dim2,nele_dim2,nedge_dim2,four_dim2,nv_dim2, &
   &one_dim2,two_dim2,time_dim2,itime_id2,i34_id2
   integer,save :: var2d_dims(2),var3d_dims(3),var4d_dims(4),dummy_dim(1), &
   &data_start_1d(1),data_start_2d(2),data_start_3d(3),data_start_4d(4), &
   &data_count_1d(1),data_count_2d(2),data_count_3d(3),data_count_4d(4)

   integer,save :: ifile,ihfskip,nspool,nc_out,nvrt,nproc_compute,np_global,ne_global,ns_global, &
   &np_max,ne_max,ns_max,ncount_2dnode,ncount_2delem,ncount_2dside,ncount_3dnode,ncount_3delem,ncount_3dside, &
   &iths0,ncid_schism_2d,istart_sed_3dnode,start_year,start_month,start_day, ics,iof_ugrid
   !Output flag dim must be same as schism_init!
   integer,save :: ntrs(natrm),iof_hydro(40),iof_wwm(40),iof_cos(20),iof_fib(5), &
   &iof_sed2d(14),iof_ice(10),iof_ana(20),iof_marsh(2),counter_out_name,nout_icm_3d(2)
   real(rkind), save :: dt,h0,start_hour,utc_start
   character(len=20), save :: out_name(max_ncoutvar)
   integer, save :: iout_23d(max_ncoutvar)
   character(len=1000),save :: out_dir
   character(len=48),save :: start_time
   character(len=48), save :: isotimestring

   integer,save,allocatable :: np(:),ne(:),ns(:),iplg(:,:),ielg(:,:),islg(:,:),kbp00(:), &
   &i34(:),elnode(:,:),rrqst2(:),ivar_id2(:),iof_gen(:),iof_age(:),iof_sed(:),iof_eco(:), &
   &iof_dvd(:),isidenode(:,:)
   real(rkind),save,allocatable :: xnd(:),ynd(:),dp(:),xel(:),yel(:),xsd(:),ysd(:)
   real(4),save,allocatable :: var2dnode(:,:,:),var2dnode_gb(:,:),var2delem(:,:,:),var2delem_gb(:,:), &
   &var2dside(:,:,:),var2dside_gb(:,:),var3dnode(:,:,:),var3dnode_gb(:,:),var3dside(:,:,:),var3dside_gb(:,:), &
   &var3delem(:,:,:),var3delem_gb(:,:)

   ! MPI communication
!    integer,save,allocatable :: rrqst2(:)

!    ! Time and file management
!    character(len=1000),save :: out_dir
!    character(len=48),save :: start_time, isotimestring
!    integer,save :: start_year,start_month,start_day,iths0
!    real(rkind), save :: start_hour,utc_start

!    ! 2D output file management (all 2D vars in one file)
!    integer,save :: ncid_schism_2d
!    integer,save :: node_dim2,nele_dim2,nedge_dim2,time_dim2,itime_id2

   public :: scribe_init
   public :: scribe_step
   public :: scribe_finalize

contains

   !===============================================================================
   ! OutputVariable Type Methods
   !===============================================================================

   subroutine var_initialize(self, name, id, dim, loc, ncomp, i23d, module_name, mod_idx)
      class(OutputVariable), intent(inout) :: self
      character(len=*), intent(in) :: name, module_name
      integer, intent(in) :: id, dim, loc, ncomp, i23d, mod_idx

      self%name = trim(name)
      self%id = id
      self%dimension = dim
      self%location = loc
      self%num_components = ncomp
      self%i23d_flag = i23d
      self%module_name = trim(module_name)
      self%module_index = mod_idx
      self%file_open = .false.

      ! Set horizontal size based on location
      select case(loc)
       case(1)
         self%horizontal_size = np_global
       case(2)
         self%horizontal_size = ne_global
       case(3)
         self%horizontal_size = ns_global
      end select

      ! Set vertical size
      if (dim == 2) then
         self%vertical_size = 1
      else
         self%vertical_size = nvrt
      endif

      ! Assign data pointer to appropriate shared array
      !   if (dim == 2) then
      !      select case(loc)
      !       case(1)
      !         self%data_2d => shared_2d_node
      !       case(2)
      !         self%data_2d => shared_2d_elem
      !       case(3)
      !         self%data_2d => shared_2d_side
      !      end select
      !   else  ! 3D
      !      select case(loc)
      !       case(1)
      !         self%data_3d => shared_3d_node
      !       case(2)
      !         self%data_3d => shared_3d_elem
      !       case(3)
      !         self%data_3d => shared_3d_side
      !      end select
      !   endif

   end subroutine var_initialize

   !---------------------------------------------------------------------------
   subroutine var_assign_scribe(self)
      class(OutputVariable), intent(inout) :: self

      if (self%dimension == 2) then
         self%assigned_scribe = nscribes  ! Last scribe
      else
         self%assigned_scribe = mod(self%id, nscribes)
      endif

      ! Set MPI tag (must be unique)
      self%mpi_tag_base = 200 + self%id

   end subroutine var_assign_scribe

   !---------------------------------------------------------------------------
   subroutine var_receive_data(self)
      class(OutputVariable), intent(inout) :: self
      integer :: i, ierr, iproc, ipgb, iegb, isgb

      ! Only the assigned scribe receives data
      if (myrank_scribe /= self%assigned_scribe) return

      if (self%dimension == 2) then
         ! Receive 2D data from all compute ranks
         !  do iproc = 1, nproc_compute
         !     if (self%location == 1) then  ! node
         !        call mpi_irecv(var2dnode(:,:,iproc), np(iproc), MPI_REAL4, &
         !           iproc-1, self%mpi_tag_base, comm_schism, rrqst2(iproc), ierr)
         !     else if (self%location == 2) then  ! elem
         !        call mpi_irecv(var2delem(:,:,iproc), ne(iproc), MPI_REAL4, &
         !           iproc-1, self%mpi_tag_base, comm_schism, rrqst2(iproc), ierr)
         !     else  ! side
         !        call mpi_irecv(var2dside(:,:,iproc), ns(iproc), MPI_REAL4, &
         !           iproc-1, self%mpi_tag_base, comm_schism, rrqst2(iproc), ierr)
         !     endif
         !  enddo
         !  call mpi_waitall(nproc_compute, rrqst2, MPI_STATUSES_IGNORE, ierr)

         !  ! Assemble into global array using index mapping
         !  do iproc = 1, nproc_compute
         !     if (self%location == 1) then
         !        var2dnode_gb(iplg(1:np(iproc),iproc), :) = transpose(var2dnode(:,1:np(iproc),iproc))
         !     else if (self%location == 2) then
         !        var2delem_gb(iplg(1:np(iproc),iproc), :) = transpose(var2delem(:,1:np(iproc),iproc))
         !     else
         !        var2dside_gb(iplg(1:np(iproc),iproc), :) = transpose(var2dside(:,1:np(iproc),iproc))
         !     endif
         !  enddo

      else  ! 3D
         !  print *, 'Rank: ', myrank_schism, ' Receiving data for ', self%name
         !  print *, 'tag = ', self%mpi_tag_base
         ! Receive 3D data from all compute ranks
         do iproc = 1, nproc_compute
            if (self%location == 1) then  ! node
               call mpi_irecv(var3dnode(:,:,iproc), np(iproc)*nvrt, MPI_REAL4, &
                  iproc-1, self%mpi_tag_base, comm_schism, rrqst2(iproc), ierr)
            else if (self%location == 2) then  ! elem
               call mpi_irecv(var3delem(:,:,iproc), ne(iproc)*nvrt, MPI_REAL4, &
                  iproc-1, self%mpi_tag_base, comm_schism, rrqst2(iproc), ierr)
            else  ! side
               call mpi_irecv(var3dside(:,:,iproc), ns(iproc)*nvrt, MPI_REAL4, &
                  iproc-1, self%mpi_tag_base, comm_schism, rrqst2(iproc), ierr)
            endif
         enddo
         call mpi_waitall(nproc_compute, rrqst2, MPI_STATUSES_IGNORE, ierr)
         if(ierr/=MPI_SUCCESS) call parallel_abort('new_scribe: mpi_waitall receieve data',ierr)
         !  print *, 'Rank: ', myrank_schism, ' Got all data for ', self%name

         ! Assemble into global array
         do iproc = 1, nproc_compute
            if (self%location == 1) then
               var3dnode_gb(1:nvrt, iplg(1:np(iproc),iproc)) = var3dnode(1:nvrt, 1:np(iproc), iproc)
            else if (self%location == 2) then
               var3delem_gb(1:nvrt, ielg(1:ne(iproc),iproc)) = var3delem(1:nvrt, 1:ne(iproc), iproc)
            else
               var3dside_gb(1:nvrt, islg(1:ns(iproc),iproc)) = var3dside(1:nvrt, 1:ns(iproc), iproc)
            endif
         enddo

         ! Fill below-bottom values for node variables
         if (self%location == 1) then
            do i = 1, np_global
               var3dnode_gb(1:kbp00(i)-1, i) = NF90_FILL_FLOAT
            enddo
         endif
      endif

   end subroutine var_receive_data

   !---------------------------------------------------------------------------
   subroutine var_write_output(self, it)
      class(OutputVariable), intent(inout) :: self
      integer, intent(in) :: it

      ! Only the assigned scribe writes
      if (myrank_scribe /= self%assigned_scribe) return

      if (self%dimension == 2) then
         ! 2D variables written together to one file - handled separately
         ! This method not called for 2D vars in current design
         return
      else
         ! 3D variables get individual files
         call write_3d_variable(self, it)
      endif

   end subroutine var_write_output

   !---------------------------------------------------------------------------
   subroutine var_cleanup(self)
      class(OutputVariable), intent(inout) :: self
      integer :: iret

      if (self%file_open) then
         ! Close NetCDF file if open
         iret =  nf90_close(self%ncid)
      endif

   end subroutine var_cleanup

   !===============================================================================
   ! Main Public Interface
   !===============================================================================

   subroutine scribe_init(indir, iths, ntime)
      character(len=*), intent(in) :: indir
      integer, intent(out) :: iths, ntime

      ! [Existing initialization code - receive parameters, mesh data, etc.]
      ! ... (similar to original, but storing in module variables)

      character(len=1000) :: var_nm2
      integer :: i,j,m,rrqst,ierr,itmp,ipgb,iegb,isgb,ivar
      integer,allocatable :: iwork(:),iwork2(:,:),iwork3(:,:)
      real(rkind),allocatable :: work(:)

      nproc_compute=nproc_schism-nscribes
      allocate(np(nproc_compute),ne(nproc_compute),ns(nproc_compute),rrqst2(nproc_compute))

      out_dir=trim(adjustl(indir))//'/outputs/'
      if(myrank_scribe==0) then
         open(16,file=trim(adjustl(out_dir))//'mirror.out.scribe',status='replace')
      endif !myrank_scribe

      !Get basic info
      call mpi_recv(dt,1,rtype,0,100,comm_schism,rrqst,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('scribe_init: recv error')
      call mpi_recv(nspool,1,itype,0,101,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_2dnode,1,itype,0,102,comm_schism,rrqst,ierr)
      call mpi_recv(nc_out,1,itype,0,103,comm_schism,rrqst,ierr)
      call mpi_recv(nvrt,1,itype,0,104,comm_schism,rrqst,ierr)
      call mpi_recv(np_global,1,itype,0,105,comm_schism,rrqst,ierr)
      call mpi_recv(ne_global,1,itype,0,106,comm_schism,rrqst,ierr)
      call mpi_recv(ns_global,1,itype,0,107,comm_schism,rrqst,ierr)
      call mpi_recv(ihfskip,1,itype,0,108,comm_schism,rrqst,ierr)
      call mpi_recv(counter_out_name,1,itype,0,109,comm_schism,rrqst,ierr)
      if(counter_out_name>max_ncoutvar) call parallel_abort('scribe_init: increase out_name dim')
      call mpi_recv(iths,1,itype,0,110,comm_schism,rrqst,ierr)
      call mpi_recv(ntime,1,itype,0,111,comm_schism,rrqst,ierr)
      call mpi_recv(iof_hydro,40,itype,0,112,comm_schism,rrqst,ierr)
      call mpi_recv(iof_wwm,40,itype,0,113,comm_schism,rrqst,ierr)
      !Make sure char len is 20 in schism_init and nc_writeout2D()!
      call mpi_recv(out_name,counter_out_name*20,MPI_CHAR,0,114,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_2delem,1,itype,0,115,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_2dside,1,itype,0,116,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_3dnode,1,itype,0,117,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_3dside,1,itype,0,118,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_3delem,1,itype,0,119,comm_schism,rrqst,ierr)
      call mpi_recv(iout_23d,counter_out_name,itype,0,120,comm_schism,rrqst,ierr)
      call mpi_recv(h0,1,rtype,0,121,comm_schism,rrqst,ierr)
      call mpi_recv(ntrs,natrm,itype,0,122,comm_schism,rrqst,ierr)
      call mpi_recv(iof_cos,20,itype,0,124,comm_schism,rrqst,ierr)
      call mpi_recv(iof_fib,5,itype,0,125,comm_schism,rrqst,ierr)
      call mpi_recv(iof_sed2d,14,itype,0,126,comm_schism,rrqst,ierr)
      call mpi_recv(iof_ice,10,itype,0,127,comm_schism,rrqst,ierr)
      call mpi_recv(iof_ana,20,itype,0,128,comm_schism,rrqst,ierr)
      call mpi_recv(iof_marsh,2,itype,0,129,comm_schism,rrqst,ierr)
#ifdef USE_ICM
      call mpi_recv(nout_icm_3d,2,itype,0,142,comm_schism,rrqst,ierr)
      !call mpi_recv(nout_d3d,1,itype,0,143,comm_schism,rrqst,ierr)
      !allocate(iof_icm(nout_icm))
      !call mpi_recv(iof_icm,nout_icm,itype,0,144,comm_schism,rrqst,ierr)
#endif
      call mpi_recv(ics,1,itype,0,146,comm_schism,rrqst,ierr)
      call mpi_recv(iof_ugrid,1,itype,0,147,comm_schism,rrqst,ierr)

      if(myrank_scribe==0) then
         write(16,*)'Scribe ',myrank_scribe,myrank_schism,nproc_scribe,nproc_compute
         write(16,*)'Scribe, basic info:',dt,nspool,nvrt,np_global,ihfskip, &
         &iths,ntime,iof_hydro,ncount_2dnode,ncount_2delem,ncount_2dside,ncount_3dnode, &
         &ncount_3dside,ncount_3delem,ntrs
         write(16,*)'out_name and i23d:'
         do i=1,counter_out_name
            write(16,*)i,trim(adjustl(out_name(i))),iout_23d(i)
         enddo !i
      endif !myrank_scribe

      !Finish rest of recv for modules
      allocate(iof_gen(max(1,ntrs(3))),iof_age(max(1,ntrs(4))),iof_sed(3*ntrs(5)+20), &
      &iof_eco(max(1,ntrs(6))),iof_dvd(max(1,ntrs(12))))
      call mpi_recv(iof_gen,max(1,ntrs(3)),itype,0,130,comm_schism,rrqst,ierr)
      call mpi_recv(iof_age,max(1,ntrs(4)),itype,0,131,comm_schism,rrqst,ierr)
      call mpi_recv(iof_sed,3*ntrs(5)+20,itype,0,132,comm_schism,rrqst,ierr)
      call mpi_recv(iof_eco,max(1,ntrs(6)),itype,0,133,comm_schism,rrqst,ierr)
      call mpi_recv(iof_dvd,max(1,ntrs(12)),itype,0,134,comm_schism,rrqst,ierr)
      call mpi_recv(istart_sed_3dnode,1,itype,0,135,comm_schism,rrqst,ierr)
      call mpi_recv(start_year,1,itype,0,136,comm_schism,rrqst,ierr)
      call mpi_recv(start_month,1,itype,0,137,comm_schism,rrqst,ierr)
      call mpi_recv(start_day,1,itype,0,138,comm_schism,rrqst,ierr)
      call mpi_recv(start_hour,1,rtype,0,139,comm_schism,rrqst,ierr)
      call mpi_recv(utc_start,1,rtype,0,140,comm_schism,rrqst,ierr)

      !Write start time into a string for later write
      !> @todo fix fractional start_hour and utc_start
      write(start_time,'(i5,2(1x,i2),2(1x,f10.2))')start_year,start_month,start_day,start_hour,utc_start
      write(isotimestring,'(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A)') 'seconds since ', start_year, '-', start_month, &
         '-', start_day, 'T', int(start_hour), ':00:00'
      if (utc_start < 0) then
         write(isotimestring,'(A,I3.2)') trim(isotimestring),int(utc_start)
      else
         write(isotimestring,'(A,A,I2.2)') trim(isotimestring),'+', int(utc_start)
      endif


      iths0=iths !save to global var

      !Last scribe receives subdomain info and then bcast
      if(myrank_schism==nproc_schism-1) then
         !First 3 dimensions. I don't understand why I had to use irecv on some
         !systems
         do i=1,nproc_compute
            call mpi_irecv(itmp,1,itype,i-1,199,comm_schism,rrqst,ierr)
            call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
            np(i)=itmp
         enddo !i

         do i=1,nproc_compute
            call mpi_irecv(itmp,1,itype,i-1,198,comm_schism,rrqst,ierr)
            call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
            ne(i)=itmp
         enddo !i

         do i=1,nproc_compute
            call mpi_irecv(itmp,1,itype,i-1,197,comm_schism,rrqst,ierr)
            call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
            ns(i)=itmp
         enddo !i

         !        write(99,*)'np:',np
         !        write(99,*)'ne:',ne
         !        write(99,*)'ns:',ns
      endif !myrank_schism

      !Max dim
      call mpi_bcast(np,nproc_compute,itype,nscribes-1,comm_scribe,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('scribe_init: mpi_bcast')
      call mpi_bcast(ne,nproc_compute,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(ns,nproc_compute,itype,nscribes-1,comm_scribe,ierr)
      np_max=maxval(np)
      ne_max=maxval(ne)
      ns_max=maxval(ns)
      if(min(np_max,ne_max,ns_max)<1) call parallel_abort('scribe_init: dim<1')
      if(myrank_scribe==0) write(16,*)'max dim:',np_max,ne_max,ns_max,myrank_schism

      !Alloc
      allocate(iplg(np_max,nproc_schism),ielg(ne_max,nproc_schism),islg(ns_max,nproc_schism))
      allocate(iwork(max(np_max,ne_max)),iwork2(4,ne_max),iwork3(2,ns_max),work(np_max),xnd(np_global), &
      &ynd(np_global),dp(np_global),kbp00(np_global),i34(ne_global),elnode(4,ne_global),isidenode(2,ns_global))
      allocate(xsd(ns_global),ysd(ns_global),xel(ne_global),yel(ne_global))
      elnode=-1 !init
      if(myrank_schism==nproc_schism-1) then
         !Mapping index arrays first. Do not combine multiple into same loop
         do i=1,nproc_compute
            call mpi_irecv(iplg(1,i),np(i),itype,i-1,196,comm_schism,rrqst2(i),ierr)
         enddo !i
         call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)

         do i=1,nproc_compute
            call mpi_irecv(ielg(1,i),ne(i),itype,i-1,195,comm_schism,rrqst2(i),ierr)
         enddo !i
         call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)

         do i=1,nproc_compute
            call mpi_irecv(islg(1,i),ns(i),itype,i-1,194,comm_schism,rrqst2(i),ierr)
         enddo !i
         call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)

         !write(99,*)'ielg 36:',ne(36),ielg(1:ne(36),36)

         !check local and global dims
         do i=1,nproc_compute
            if(maxval(iplg(1:np(i),i))>np_global.or.minval(iplg(1:np(i),i))<1) &
            &call parallel_abort('scribe_init: overflow(1)')
            if(maxval(ielg(1:ne(i),i))>ne_global.or.minval(ielg(1:ne(i),i))<1) &
            &call parallel_abort('scribe_init: overflow(2)')
            if(maxval(islg(1:ns(i),i))>ns_global.or.minval(islg(1:ns(i),i))<1) &
            &call parallel_abort('scribe_init: overflow(3)')
         enddo !i

         !Other vars using index arrays
         do i=1,nproc_compute
            call mpi_irecv(work,np(i),rtype,i-1,193,comm_schism,rrqst,ierr)
            call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
            xnd(iplg(1:np(i),i))=work(1:np(i))
         enddo !i
         do i=1,nproc_compute
            call mpi_irecv(work,np(i),rtype,i-1,192,comm_schism,rrqst,ierr)
            call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
            ynd(iplg(1:np(i),i))=work(1:np(i))
         enddo !i
         do i=1,nproc_compute
            call mpi_irecv(work,np(i),rtype,i-1,191,comm_schism,rrqst,ierr)
            call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
            dp(iplg(1:np(i),i))=work(1:np(i))
         enddo !i
         do i=1,nproc_compute
            call mpi_irecv(iwork,np(i),itype,i-1,190,comm_schism,rrqst,ierr)
            call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
            kbp00(iplg(1:np(i),i))=iwork(1:np(i))
         enddo !i
         do i=1,nproc_compute
            call mpi_irecv(iwork,ne(i),itype,i-1,189,comm_schism,rrqst,ierr)
            call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
            i34(ielg(1:ne(i),i))=iwork(1:ne(i))
         enddo !i
         do i=1,nproc_compute
            call mpi_irecv(iwork2,4*ne(i),itype,i-1,188,comm_schism,rrqst,ierr)
            call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)

            do j=1,ne(i)
               iegb=ielg(j,i)
               do m=1,i34(iegb)
                  elnode(m,iegb)=iplg(iwork2(m,j),i)
               enddo !m
            enddo !j
         enddo !i

         do i=1,nproc_compute
            call mpi_irecv(iwork3,2*ns(i),itype,i-1,187,comm_schism,rrqst,ierr)
            call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
            do j=1,ns(i)
               isgb=islg(j,i)
               isidenode(1,isgb)=iplg(iwork3(1,j),i)
               isidenode(2,isgb)=iplg(iwork3(2,j),i)
            enddo !j
         enddo !i

         !        write(99,*)'x:',xnd
         !        write(99,*)'y:',ynd
         !        do i=1,ne_global
         !          write(99,*)i,i34(i),elnode(1:i34(i),i)
         !        enddo !i
         !        write(98,*)'kbp00:',kbp00
      endif !myrank_schism==nproc_schism-1

      call mpi_bcast(iplg,np_max*nproc_schism,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(ielg,ne_max*nproc_schism,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(islg,ns_max*nproc_schism,itype,nscribes-1,comm_scribe,ierr)

      call mpi_bcast(xnd,np_global,rtype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(ynd,np_global,rtype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(dp,np_global,rtype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(kbp00,np_global,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(i34,ne_global,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(elnode,4*ne_global,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(isidenode,2*ns_global,itype,nscribes-1,comm_scribe,ierr)

      deallocate(work,iwork,iwork2,iwork3)

      ! Calculate side centers as edge_x/y location. Error across lon jump
      do i=1,ns_global
         xsd(i) = sum(xnd(isidenode(1:2,i))) / 2.0
         ysd(i) = sum(ynd(isidenode(1:2,i))) / 2.0
      enddo

      ! Calculate elem centers as edge_x/y location. Error across lon jump
      do i=1,ne_global
         xel(i) = sum(xnd(elnode(1:i34(i),i)))/real(i34(i),rkind)
         yel(i) = sum(ynd(elnode(1:i34(i),i)))/real(i34(i),rkind)
      enddo


      ! Count total output variables
      call count_output_variables(num_output_vars, num_2d_vars, num_3d_vars)
      !   print *, 'num output vars = ', num_output_vars, ', num 2d = ', num_2d_vars, ', num 3d = ', num_3d_vars

      ! Allocate output variable array
      allocate(output_vars(num_output_vars))

      ! Build output variable list
      call build_output_variable_list()

      ! Assign scribes to variables
      do ivar = 1, num_output_vars
         call output_vars(ivar)%assign_scribe()
      enddo

      allocate(var2dnode(ncount_2dnode,np_max,nproc_compute),var2dnode_gb(np_global,ncount_2dnode), &
      &var2delem(ncount_2delem,ne_max,nproc_compute),var2delem_gb(ne_global,ncount_2delem), &
      &var2dside(ncount_2dside,ns_max,nproc_compute),var2dside_gb(ns_global,ncount_2dside), &
      &ivar_id2(ncount_2dnode+ncount_2delem+ncount_2dside))
      var2dside_gb(ns_global,ncount_2dside)=0. !touch mem

      allocate(var3dnode(nvrt,np_max,nproc_compute),var3dnode_gb(nvrt,np_global))
      var3dnode(nvrt,np_max,nproc_compute)=0.
      var3dnode_gb(nvrt,np_global)=0. !touch mem

      if(ncount_3dside>0) then
         allocate(var3dside(nvrt,ns_max,nproc_compute),var3dside_gb(nvrt,ns_global))
         var3dside(nvrt,ns_max,nproc_compute)=0.
      endif

      if(ncount_3delem>0) then
         allocate(var3delem(nvrt,ne_max,nproc_compute),var3delem_gb(nvrt,ne_global))
         var3delem(nvrt,ne_max,nproc_compute)=0.
      endif

      ! Allocate shared data arrays
      !   allocate(shared_2d_node(np_global, 1))
      !   allocate(shared_2d_elem(ne_global, 1))
      !   allocate(shared_2d_side(ns_global, 1))
      !   allocate(shared_3d_node(nvrt, np_global))
      !   allocate(shared_3d_elem(nvrt, ne_global))
      !   allocate(shared_3d_side(nvrt, ns_global))

      ! Allocate receive buffers
      !   allocate(recv_buffer_2d(1, max(np_max,ne_max,ns_max), nproc_compute))
      !   allocate(recv_buffer_3d(nvrt, max(np_max,ne_max,ns_max), nproc_compute))

      if(myrank_scribe==0) then
         write(16,*)'Scribe initialized with',num_output_vars,'output variables'
         write(16,*)'  2D variables:',num_2d_vars
         write(16,*)'  3D variables:',num_3d_vars
         write(16,*)'Variable assignments:'
         do ivar = 1, num_output_vars
            write(16,'(A,I4,A,A,A,I2)') '  Var ',ivar,': ',trim(output_vars(ivar)%name), &
               ' -> Scribe ',output_vars(ivar)%assigned_scribe
         enddo
      endif

   end subroutine scribe_init

   !---------------------------------------------------------------------------
   subroutine scribe_step(it)
      integer, intent(in) :: it
      integer :: ivar, i, num_2d
      logical :: new_stack

      ! Return if not output step
      if (nc_out == 0 .or. mod(it, nspool) /= 0) return

      new_stack = (mod(it-nspool, ihfskip) == 0)

      ! Handle 2D output (all 2D vars written together by last scribe)
      call handle_2d_output(it, new_stack)

      ! Handle 3D output (distributed across scribes)
      do ivar = 1, num_output_vars
         if (output_vars(ivar)%dimension == 3) then
            ! Each scribe processes only its assigned variables
            if (myrank_scribe == output_vars(ivar)%assigned_scribe) then
               call output_vars(ivar)%receive_data()
               call output_vars(ivar)%write_output(it)
            endif
         endif
      enddo
      !   print *, 'Rank: ', myrank_schism, ' finished scribe step'

   end subroutine scribe_step

   !---------------------------------------------------------------------------
   subroutine scribe_finalize()
      integer :: ivar

      ! Clean up all variables
      do ivar = 1, num_output_vars
         call output_vars(ivar)%cleanup()
      enddo

      ! Deallocate arrays
      if (allocated(output_vars)) deallocate(output_vars)
      !   if (allocated(shared_2d_node)) deallocate(shared_2d_node)
      !   if (allocated(shared_2d_elem)) deallocate(shared_2d_elem)
      !   if (allocated(shared_2d_side)) deallocate(shared_2d_side)
      !   if (allocated(shared_3d_node)) deallocate(shared_3d_node)
      !   if (allocated(shared_3d_elem)) deallocate(shared_3d_elem)
      !   if (allocated(shared_3d_side)) deallocate(shared_3d_side)
      !   if (allocated(recv_buffer_2d)) deallocate(recv_buffer_2d)
      !   if (allocated(recv_buffer_3d)) deallocate(recv_buffer_3d)

   end subroutine scribe_finalize

   !===============================================================================
   ! Helper Subroutines
   !===============================================================================

   subroutine count_output_variables(total, num_2d, num_3d)
      integer, intent(out) :: total, num_2d, num_3d

      num_2d = ncount_2dnode+ncount_2delem+ncount_2dside
      num_3d = ncount_3dnode+ncount_3delem+ncount_3dside
      total = ncount_3dnode+ncount_3delem+ncount_3dside
   end subroutine count_output_variables

   !---------------------------------------------------------------------------
   subroutine build_output_variable_list()
      integer :: ivar, j, icomp
      character(len=20) :: varname, comp_suffix(2)

      comp_suffix = (/'X', 'Y'/)
      ivar = 0
      ivar = ivar + 1

      ! Build 2D hydro variables
      !   if (iof_hydro(1) /= 0) then
      !      ivar = ivar + 1
      !      call output_vars(ivar)%initialize(out_name(ivar), ivar, 2, 1, 1, 1, 'hydro', 1)
      !   endif

      !   if (iof_hydro(2) /= 0) then
      !      ivar = ivar + 1
      !      call output_vars(ivar)%initialize(out_name(ivar), ivar, 2, 1, 1, 1, 'hydro', 2)
      !      ivar = ivar + 1
      !      call output_vars(ivar)%initialize(out_name(ivar), ivar, 2, 1, 1, 1, 'hydro', 2)
      !   endif

      !   ! ... more 2D variables

      ! Build 3D hydro variables
      do j = 17, 25
         if (iof_hydro(j) /= 0) then
            ivar = ivar + 1
            ! Get variable name from lookup or naming convention
            ! varname = get_hydro_varname(j)
            ! print *, 'Initializing variable name = ', out_name(num_2d_vars + ivar - 1), ' ivar = ', ivar
            call output_vars(ivar - 1)%initialize(out_name(num_2d_vars + ivar - 1), ivar, 3, 1, 1, 3, 'hydro', j)
         endif
      enddo

      ! 3D vectors
      if (iof_hydro(26) /= 0) then  ! horizontal velocity
         do icomp = 1, 2
            ivar = ivar + 1
            ! print *, 'Initializing variable name = ', out_name(num_2d_vars + ivar - 1), ' ivar = ', ivar
            call output_vars(ivar - 1)%initialize(out_name(num_2d_vars + ivar - 1), ivar, 3, 1, 1, 3, 'hydro', 26)
         enddo
      endif

      ! [Continue for all modules]
#ifdef USE_WWM
      ! WWM variables...
      do j=35,36
         if(iof_wwm(j)/=0) then
            ivar = ivar + 1
            ! print *, 'Initializing variable name = ', out_name(num_2d_vars + ivar - 1), ' ivar = ', ivar
            call output_vars(ivar - 1)%initialize(out_name(num_2d_vars + ivar - 1), ivar, 3, 1, 1, 3, 'wwm', j)
            !  if(iof_wwm(j)/=0) call scribe_recv_write(it,1,2,itotal,icount_out_name)
         endif
      enddo !j
#endif

#ifdef USE_SED
      ! Sediment variables...
#endif

      ! etc.
      do j=27,27
         !  if(iof_hydro(j)/=0) call scribe_recv_write(it,3,2,itotal,icount_out_name)
         if(iof_hydro(j)/=0) then
            do icomp = 1, 2
               ivar = ivar + 1
               !    print *, 'Initializing variable name = ', out_name(num_2d_vars + ivar - 1), ' ivar = ', ivar
               call output_vars(ivar - 1)%initialize(out_name(num_2d_vars + ivar - 1), ivar, 3, 1, 1, 3, 'hydro', j)
            enddo
         endif
      enddo

#ifdef USE_WWM
      if(iof_wwm(33)/=0) then
         ! call scribe_recv_write(it,3,1,itotal,icount_out_name)
         ivar = ivar + 1
         !  print *, 'Initializing variable name = ', out_name(num_2d_vars + ivar - 1), ' ivar = ', ivar
         call output_vars(ivar - 1)%initialize(out_name(num_2d_vars + ivar - 1), ivar, 3, 1, 1, 3, 'wwm', 33)
      endif

      !Vector
      !   if(iof_wwm(34)/=0) call scribe_recv_write(it,3,2,itotal,icount_out_name)
      if(iof_wwm(34)/=0) then
         ! call scribe_recv_write(it,3,1,itotal,icount_out_name)
         ivar = ivar + 1
         !  print *, 'Initializing variable name = ', out_name(num_2d_vars + ivar - 1), ' ivar = ', ivar
         call output_vars(ivar - 1)%initialize(out_name(num_2d_vars + ivar -1), ivar, 3, 1, 1, 3, 'wwm', 34)
         ivar = ivar + 1
         !  print *, 'Initializing variable name = ', out_name(num_2d_vars + ivar - 1), ' ivar = ', ivar
         call output_vars(ivar - 1)%initialize(out_name(num_2d_vars + ivar - 1), ivar, 3, 1, 1, 3, 'wwm', 34)
      endif


#endif /*USE_WWM*/
      do j=28,30
         if(iof_hydro(j)/=0) then
            ! call scribe_recv_write(it,2,1,itotal,icount_out_name)
            ivar = ivar + 1
            ! print *, 'Initializing variable name = ', out_name(num_2d_vars + ivar - 1), ' ivar = ', ivar
            call output_vars(ivar - 1)%initialize(out_name(num_2d_vars + ivar - 1), ivar, 3, 2, 1, 3, 'hydro', j)
         endif
      enddo !j
      !   print *, 'iof_hydro = ', iof_hydro
      !   print *, 'iof_wwm = ', iof_wwm

      !   print *, 'build_output_variable_list final ivar = ', ivar

   end subroutine build_output_variable_list

   !---------------------------------------------------------------------------
   subroutine handle_2d_output(it, new_stack)
      integer, intent(in) :: it
      logical, intent(in) :: new_stack

      integer :: i,j,k,m,rrqst,ierr,irank,itotal,icount_out_name,itmp5

!      write(*,*)'Ready for I/O...',myrank_schism,it
      itotal=0 !# of output/sends so far

      !2D: all 2D outputs share same scribe
      itotal=itotal+1 !used in tags and rank #
      irank=nproc_schism-itotal !last scribe
      if(myrank_schism/=irank) return
!------------------
      !   print *, 'Rank: ', myrank_schism, ' handle_2d_output'
      !2D node (modules already included inside the array)
      do i=1,nproc_compute
         call mpi_irecv(var2dnode(:,:,i),np(i)*ncount_2dnode,MPI_REAL4,i-1,200+itotal,comm_schism,rrqst2(i),ierr)
      enddo !i
      call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('new_scribe: mpi_waitall 2d nodes',ierr)

      !   print *, 'Rank: ', myrank_schism, ' got 2d node vars'

      do i=1,nproc_compute
         var2dnode_gb(iplg(1:np(i),i),:)=transpose(var2dnode(:,1:np(i),i)) !indiced reversed for write
!          write(99,*)'dry:',myrank_schism,it,i,var2dnode(1,1:np(i),i)
!          write(98,*)'elev:',myrank_schism,it,i,var2dnode(2,1:np(i),i)
!          write(97,*)'wind:',myrank_schism,it,i,var2dnode(3:4,1:np(i),i)
      enddo !i
!------------------
      !2D elem (modules already included inside the array)
      do i=1,nproc_compute
         call mpi_irecv(var2delem(:,:,i),ne(i)*ncount_2delem,MPI_REAL4,i-1,701,comm_schism,rrqst2(i),ierr)
      enddo !i
      call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('new_scribe: mpi_waitall 2d elems',ierr)
      do i=1,nproc_compute
         var2delem_gb(ielg(1:ne(i),i),:)=transpose(var2delem(:,1:ne(i),i)) !indiced reversed for write
!          write(99,*)'elem dry:',myrank_schism,it,i,var2delem(:,1:ne(i),i)
      enddo !i

!------------------
      !2D side (modules already included inside the array)
      do i=1,nproc_compute
         call mpi_irecv(var2dside(:,:,i),ns(i)*ncount_2dside,MPI_REAL4,i-1,702,comm_schism,rrqst2(i),ierr)
      enddo !i
      call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('new_scribe: mpi_waitall 2d side',ierr)
      do i=1,nproc_compute
         var2dside_gb(islg(1:ns(i),i),:)=transpose(var2dside(:,1:ns(i),i)) !indiced reversed for write
!          write(98,*)'side dry:',myrank_schism,it,i,var2dside(:,1:ns(i),i)
      enddo !i

!------------------
      !Output all 2D vars (modules included inside arrays)
      call nc_writeout2D(it,np_global,ne_global,ns_global,ncount_2dnode,ncount_2delem,ncount_2dside, &
      &var2dnode_gb,var2delem_gb,var2dside_gb,out_name(1:ncount_2dnode+ncount_2delem+ncount_2dside), &
      &iout_23d(1:ncount_2dnode+ncount_2delem+ncount_2dside))

   end subroutine handle_2d_output

   !---------------------------------------------------------------------------
   subroutine write_3d_variable(var, it)
      type(OutputVariable), intent(inout) :: var
      integer, intent(in) :: it

      integer :: irec, iret, ivar_id
      character(len=140) :: fname
      character(len=12) :: ifile_char
      logical :: new_stack
      real(rkind) :: a1d(1)

      new_stack = (mod(it-nspool, ihfskip) == 0)

      if (new_stack) then
         ifile = (it-1)/ihfskip + 1
         write(ifile_char,'(i12)') ifile
         fname = trim(adjustl(out_dir))//'/'//trim(adjustl(var%name))//'_' &
            //trim(adjustl(ifile_char))//'.nc'

         ! Create file and write header
         iret = nf90_create(trim(adjustl(fname)), OR(NF90_NETCDF4,NF90_CLOBBER), var%ncid)
         call write_3d_header(var, ifile_char)
         var%file_open = .true.
      endif

      ! Write data
      irec = (it-(ifile-1)*ihfskip)/nspool
      !   call write_3d_data_record(var, irec)

      if(irec<=0) call parallel_abort('nc_writeout3D: irec<=0')
      a1d(1)=dble(it*dt)
      data_start_1d(1)=irec; data_count_1d(1)=1
      iret=nf90_put_var(var%ncid,itime_id,a1d,data_start_1d,data_count_1d)
      if(iret/=NF90_NOERR) call parallel_abort('nc_writeout3D: put time')

      data_start_3d=(/1,1,irec/)
      data_count_3d=(/nvrt,var%horizontal_size,1/)
      select case (var%location)
       case (1)
         iret=nf90_put_var(var%ncid,var%varid,real(var3dnode_gb),data_start_3d,data_count_3d)
       case (2)
         iret=nf90_put_var(var%ncid,var%varid,real(var3delem_gb),data_start_3d,data_count_3d)
       case (3)
         iret=nf90_put_var(var%ncid,var%varid,real(var3dside_gb),data_start_3d,data_count_3d)
      end select
      if(iret/=NF90_NOERR) call parallel_abort('nc_writeout3D: put 3D var')

      ! Close file if end of stack
      if (mod(it, ihfskip) == 0) then
         iret = nf90_close(var%ncid)
         var%file_open = .false.
      endif

   end subroutine write_3d_variable

   subroutine write_3d_header(var, ifile_char)
      type(OutputVariable), intent(inout) :: var
      character(len=12), intent(in) :: ifile_char

      !   integer :: irec, iret
      !   character(len=140) :: fname
      !   character(len=12) :: ifile_char
      !   logical :: new_stack
      integer :: irec,iret, ivarid,chunks(3)

      !Define chunk size (contiguous block for access) for deflation: each chunk
      !must be < 4GB in size
      !TODO: add a scale for chunks(2) for large meshes
      chunks(1)=nvrt; chunks(3)=1
      if (var%location==1) then
         chunks(2) = np_global
      elseif(var%location==2) then
         chunks(2) = ne_global
      else
         chunks(2) = ns_global
      end if
      !Outputs are in 4 bytes; 4GB is in binary not decimal - not precise
      if(4.d0*chunks(1)*chunks(2)*chunks(3)>3.d9) call parallel_abort('nc_writeout3D: chunk size')

      call fill_header_static(iof_ugrid,var%ncid,itime_id,node_dim, &
      &nele_dim,nedge_dim,four_dim,nv_dim,one_dim,two_dim,time_dim)

      iret=nf90_redef(var%ncid)

      var3d_dims(1)=nv_dim
      var3d_dims(3)=time_dim
      if(var%location==1) then
         var3d_dims(2)=node_dim
      elseif(var%location==2) then
         var3d_dims(2)=nele_dim
      else !3
         var3d_dims(2)=nedge_dim
      endif
      iret=nf90_def_var(var%ncid,trim(adjustl(var%name)),NF90_FLOAT,var3d_dims,var%varid)
      if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout3D: var_dims')
      iret=nf90_put_att(var%ncid,var%varid,'i23d',var%i23d_flag) !set i23d flag
      iret=nf90_put_att(var%ncid,var%varid,'missing_value',NF90_FILL_FLOAT)

      iret=nf90_def_var_chunking(var%ncid,var%varid,NF90_CHUNKED,chunks)
      !function nf90_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level)
      !where deflate_level\in[0,9] with 9 being most compression
      iret=nf90_def_var_deflate(var%ncid,var%varid,0,1,4)


      call add_mesh_attributes(var%ncid, var%varid, iof_ugrid)
      if (iof_ugrid > 0) then
         ! We write CF metadata for iof_ugrid > 0.  In the case of iof_ugrid ==2, we do not
         ! duplicate the mesh information in the 3D output files
         call add_cf_variable_attributes(var%ncid, var%varid)
      endif

      if (iof_ugrid == 2) then
         ! The UGRID information is only available in out_2d files for iof_ugrid == 2, but it is not
         ! written to the 3D output files. There, it is referenced by the external_variables attribute.
         ! This global external_variables attribute is a blank-separated list of the names of variables
         ! which are named by attributes in the file but which are not present in the file.
         iret = nf90_put_att(var%ncid, NF90_GLOBAL, 'external_variables', 'SCHISM_hgrid_node_x SCHISM_hgrid_node_y')
         if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout3D: external_variables')
         iret = nf90_put_att(var%ncid, NF90_GLOBAL, 'external_variables_location', &
         & 'out2d_'//trim(adjustl(ifile_char))//'.nc')
         if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout3D: external_variables_location')
      endif !iof_ugrid == 2

      call add_user_attributes(var%ncid)

      iret=nf90_enddef(var%ncid)


   end subroutine write_3d_header

   ! [Other helper functions: write_2d_file, write_3d_header, write_3d_data_record, etc.]
   ! These would contain the actual NetCDF writing logic from the original module

   subroutine nc_writeout2D(it,np_gb,ne_gb,ns_gb,ncount_p,ncount_e,ncount_s, &
   &var2dnode_gb2,var2delem_gb2,var2dside_gb2,vname,i23da)
      implicit none

      integer, intent(in) :: it,np_gb,ne_gb,ns_gb,ncount_p,ncount_e,ncount_s,i23da(ncount_p+ncount_e+ncount_s)
      real(4), intent(in) :: var2dnode_gb2(np_gb,ncount_p),var2delem_gb2(ne_gb,ncount_e),var2dside_gb2(ns_gb,ncount_s)
      character(len=20), intent(in) :: vname(ncount_p+ncount_e+ncount_s)

      integer :: irec,iret,i,j,k,ih0_id2,ikbp_id2, ivarid
      character(len=140) :: fname
      character(len=12) :: ifile_char
      real(rkind) :: a1d(1)

      if(mod(it-nspool,ihfskip)==0) then !new stack
         ifile=(it-1)/ihfskip+1  !output stack #
         write(ifile_char,'(i12)') ifile
         fname=trim(adjustl(out_dir))//'/out2d_'//trim(adjustl(ifile_char))//'.nc'
         iret=nf90_create(trim(adjustl(fname)),OR(NF90_NETCDF4,NF90_CLOBBER),ncid_schism_2d)

         !Fill in header and static info
         call fill_header_static(1,ncid_schism_2d,itime_id2,node_dim2, &
         &nele_dim2,nedge_dim2,four_dim2,nv_dim2,one_dim2,two_dim2,time_dim2)

         iret=nf90_redef(ncid_schism_2d)
         !> Deal with all the variables with time/node dimension
         do i=1,ncount_p
            var2d_dims(1)=node_dim2
            var2d_dims(2)=time_dim2
            iret=nf90_def_var(ncid_schism_2d,trim(adjustl(vname(i))),NF90_FLOAT,var2d_dims,ivar_id2(i))
            if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: var_dims')
            iret=nf90_put_att(ncid_schism_2d,ivar_id2(i),'i23d',i23da(i)) !set i23d flag
            iret=nf90_def_var_deflate(ncid_schism_2d,ivar_id2(i),0,1,4)
            call add_mesh_attributes(ncid_schism_2d,ivar_id2(i), iof_ugrid)
            if(iof_ugrid/=0) call add_cf_variable_attributes(ncid_schism_2d,ivar_id2(i))
         enddo !i

         do i=1,ncount_e
            var2d_dims(1)=nele_dim2; var2d_dims(2)=time_dim2
            iret=nf90_def_var(ncid_schism_2d,trim(adjustl(vname(i+ncount_p))),NF90_FLOAT,var2d_dims,ivar_id2(i+ncount_p))
            if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: var_dims(2)')
            iret=nf90_put_att(ncid_schism_2d,ivar_id2(i+ncount_p),'i23d',i23da(i+ncount_p)) !set i23d flag
            iret=nf90_def_var_deflate(ncid_schism_2d,ivar_id2(i+ncount_p),0,1,4)
            call add_mesh_attributes(ncid_schism_2d,ivar_id2(i+ncount_p), iof_ugrid)
            if(iof_ugrid/=0) call add_cf_variable_attributes(ncid_schism_2d,ivar_id2(i+ncount_p))
         enddo !i

         do i=1,ncount_s
            var2d_dims(1)=nedge_dim2; var2d_dims(2)=time_dim2
            iret=nf90_def_var(ncid_schism_2d,trim(adjustl(vname(i+ncount_p+ncount_e))),NF90_FLOAT, &
            &var2d_dims,ivar_id2(i+ncount_p+ncount_e))
            if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: var_dims(3)')
            iret=nf90_put_att(ncid_schism_2d,ivar_id2(i+ncount_p+ncount_e),'i23d', &
            &i23da(i+ncount_p+ncount_e)) !set i23d flag
            iret=nf90_def_var_deflate(ncid_schism_2d,ivar_id2(i+ncount_p+ncount_e),0,1,4)
            call add_mesh_attributes(ncid_schism_2d,ivar_id2(i+ncount_p+ncount_e), iof_ugrid)
            if(iof_ugrid/=0) call add_cf_variable_attributes(ncid_schism_2d,ivar_id2(i+ncount_p+ncount_e))
         enddo !i

         call add_user_attributes(ncid_schism_2d)
         iret=nf90_enddef(ncid_schism_2d)
      endif !mod(it-

      !Output
      irec=(it-(ifile-1)*ihfskip)/nspool !time recod #
!      print*, 'inside nc_writeout2D:',irec,it
      if(irec<=0) call parallel_abort('nc_writeout2D: irec<=0')
      a1d(1)=dble(it*dt)
      data_start_1d(1)=irec; data_count_1d(1)=1
      iret=nf90_put_var(ncid_schism_2d,itime_id2,a1d,data_start_1d,data_count_1d)
      if(iret/=NF90_NOERR) call parallel_abort('nc_writeout2D: put time')

      data_start_2d=(/1,irec/)
      data_count_2d=(/np_global,1/)
      do i=1,ncount_p
         iret=nf90_put_var(ncid_schism_2d,ivar_id2(i),real(var2dnode_gb2(:,i)),data_start_2d,data_count_2d)
         if(iret/=NF90_NOERR) call parallel_abort('nc_writeout2D: put 2D var')
      enddo !i
      data_count_2d=(/ne_global,1/)
      do i=1,ncount_e
         iret=nf90_put_var(ncid_schism_2d,ivar_id2(i+ncount_p),real(var2delem_gb2(:,i)),data_start_2d,data_count_2d)
         if(iret/=NF90_NOERR) call parallel_abort('nc_writeout2D: put 2D elem var')
      enddo !i
      data_count_2d=(/ns_global,1/)
      do i=1,ncount_s
         iret=nf90_put_var(ncid_schism_2d,ivar_id2(i+ncount_p+ncount_e),real(var2dside_gb2(:,i)),data_start_2d,data_count_2d)
         if(iret/=NF90_NOERR) call parallel_abort('nc_writeout2D: put 2D elem var')
      enddo !i

      !Close the stack
      if(mod(it,ihfskip)==0) iret=nf90_close(ncid_schism_2d)

   end subroutine nc_writeout2D

   subroutine add_cf_variable_attributes(ncid, varid)
      implicit none
      integer, intent(in) :: ncid, varid

      integer :: iret, ndims, i
      character(len=255)   :: varname, long_name, standard_name, units, comment, cell_methods
      real :: valid_min, valid_max

      iret = nf90_inquire_variable(ncid, varid, name=varname, ndims=ndims)
      if(iret.ne.NF90_NOERR) call parallel_abort('add_cf_variable_attributes: inquire_variable')

      !write(16, '(A,A)') 'add_cf_variable_attributes: ', trim(varname)

      long_name = ''
      standard_name = ''
      comment = ''
      units = ''
      valid_min = NF90_FILL_FLOAT
      valid_max = NF90_FILL_FLOAT

      select case (trim(varname))
       case ('airPressure')
         long_name = 'Surface air pressure'
         standard_name = 'surface_air_pressure'
         units = 'Pa'
       case ('airTemperature')
         long_name = 'Surface air temperature'
         standard_name = 'air_temperature'
         units = 'degree_Celsius'
         ! units_metadata = 'on_scale'
       case ('depth')
         long_name = 'Depth below sea surface'
         standard_name = 'depth'
         units = 'm'
         !positive = 'down'
         comment = 'Positive downward'
         valid_min = -500000.0
         valid_max = 20000.0
       case ('depthAverageVelY')
         long_name = 'Depth-average y velocity'
         standard_name = 'sea_water_y_velocity'
         units = 'm s-1'
         valid_min = -1000.0
         valid_max = 1000.0
         cell_methods = 'nSCHISM_vgrid_layers: mean'
       case ('depthAverageVelX')
         long_name = 'Depth-average x velocity'
         standard_name = 'sea_water_x_velocity'
         units = 'm s-1'
         valid_min = -1000.0
         valid_max = 1000.0
         cell_methods = 'nSCHISM_vgrid_layers: mean'
       case ('zCoordinates')
         standard_name = 'depth'
         long_name = 'Layer interface depth'
         units = 'm'
       case ('waterDensity')
         standard_name = 'sea_water_density'
         units = 'kg m-3'
       case ('barotropicPresGradX')
         long_name = 'Eastward barotropic pressure gradient'
         units = 'Pa m-1'
       case ('barotropicPresGradY')
         long_name = 'Northward barotropic pressure gradient'
         units = 'Pa m-1'
       case ('temperatureAtElement')
         standard_name = 'sea_water_temperature'
         units = 'degree_C'
       case ('salinityAtElement')
         standard_name = 'sea_water_salinity'
         units = 'g kg-1'
       case ('verticalVelAtElement')
         standard_name = 'upward_sea_water_velocity'
         units = 'm s-1'
       case ('turbulentKineticEner')
         standard_name = 'specific_turbulent_kinetic_energy_of_sea_water '
         units = 'm2 s-2'
       case ('mixingLength')
         standard_name = 'turbulent_mixing_length_of_sea_water'
         units = 'm'
       case ('verticalVelocity')
         ! Check this one
         standard_name = 'upward_sea_water_velocity'
         units = 'm s-1'
       case ('temperature')
         standard_name = 'sea_water_temperature'
         units = 'degree_C'
       case ('salinity')
         standard_name = 'sea_water_salinity'
         units = 'g kg-1'
       case ('diffusivity')
         standard_name = 'ocean_vertical_diffusivity'
         units = 'm2 s-1'
       case ('viscosity')
         long_name = 'Viscosity'
         units = 'N s m-2'
       case ('elevation')
         standard_name = 'sea_surface_height_above_geoid'
         units = 'm'
       case ('specificHumidity')
         standard_name = 'specific_humidity'
         units = 'kg kg-1'
       case ('solarRadiation')
         standard_name = 'net_downward_shortwave_flux_at_sea_water_surface'
         units = 'W m-2'
       case ('latentHeat')
         standard_name = 'surface_downward_latent_heat_flux'
         units = 'W m-2'
       case ('sensibleHeat')
         standard_name = 'surface_downward_sensible_heat_flux'
         units = 'W m-2'
       case ('totalHeat')
         standard_name = 'surface_downward_heat_flux_in_sea_water'
         units = 'W m-2'
       case ('evaporationRate')
         standard_name = 'water_evapotranspiration_flux'
         units = 'kg m-2 s-1'
       case ('precipitationRate')
         standard_name = 'precipitation_flux'
         units = 'kg m-2 s-1'
       case ('upwardLongwave')
         standard_name = 'upwelling_longwave_flux_in_air'
         units = 'W m-2'
       case ('downwardLongwave')
         standard_name = 'downwelling_longwave_flux_in_air'
         units = 'W m-2'
       case ('bottomStressX')
         standard_name = 'sea_floor_horizontal_stress'
         long_name = 'Bottom stress in x direction'
         units = 'N m-2'
       case ('bottomStressY')
         standard_name = 'sea_floor_horizontal_stress'
         long_name = 'Bottom stress in y direction'
         units = 'N m-2'
       case ('horizontalSideVelX')
         standard_name = 'sea_water_x_velocity'
         units = 'm s-1'
       case ('horizontalSideVelY')
         standard_name = 'sea_water_y_velocity'
         units = 'm s-1'
       case ('windSpeedX')
         standard_name = 'eastward_wind'
         units = 'm s-1'
       case ('windSpeedY')
         standard_name = 'northward_wind'
         units = 'm s-1'
       case ('windStressX')
         standard_name = 'downward_x_stress_at_sea_water_surface'
         units = 'N m-2'
       case ('windStressY')
         standard_name = 'downward_y_stress_at_sea_water_surface'
         units = 'N m-2'
       case ('bottom_index_node')
         long_name = 'Bottom index at node'
         standard_name = 'model_level_number_at_sea_floor'
         units = '1'
         !positive = 'up'
         comment = 'Index of lowermost model level'
         valid_min = 1
       case default
      end select

      if (index(trim(varname), 'dry') > 0 .and. index(trim(varname), 'Flag') > 0) then
         standard_name = 'quality_flag'
         long_name = 'Dry flag'
         comment = '1: dry, 0: wet'
         units = '1'
         valid_min = 0
         valid_max = 1
      endif

      iret = NF90_NOERR
      if (trim(long_name) /= '')  iret=nf90_put_att(ncid,varid,'long_name',trim(long_name))
      if(iret.ne.NF90_NOERR) call parallel_abort(varname//' long_name')

      if (trim(standard_name) /= '')  iret=nf90_put_att(ncid,varid,'standard_name',trim(standard_name))
      if(iret.ne.NF90_NOERR) call parallel_abort(varname//' standard_name')

      if (trim(units) /= '')  iret=nf90_put_att(ncid,varid,'units',trim(units))
      if(iret.ne.NF90_NOERR) call parallel_abort(varname//' units')

      if (trim(comment) /= '')  iret=nf90_put_att(ncid,varid,'comment',trim(comment))
      if(iret.ne.NF90_NOERR) call parallel_abort(varname//' comment')

      if (valid_min /= NF90_FILL_FLOAT)  iret=nf90_put_att(ncid,varid,'valid_min',valid_min)
      if(iret.ne.NF90_NOERR) call parallel_abort(varname//' valid_min')

      if (valid_max /= NF90_FILL_FLOAT)  iret=nf90_put_att(ncid,varid,'valid_max',valid_max)
      if(iret.ne.NF90_NOERR) call parallel_abort(varname//' valid_max')

   end subroutine add_cf_variable_attributes

   subroutine add_user_attributes(ncid)
      implicit none
      integer, intent(in) :: ncid

      integer :: iret, varid, ncin, nvar, ivar, natt, dimid, iatt
      character(len=11), parameter :: infilename = 'metadata.nc'
      character(len=NF90_MAX_NAME) :: varname, attname
      logical :: file_exists

      inquire(file=infilename, exist=file_exists)
      if (.not. file_exists) return

      iret = nf90_open(trim(infilename), NF90_NOWRITE, ncin)
      if(iret.ne.NF90_NOERR) call parallel_abort('add_user_attributes: nf90_open '//trim(infilename))

      iret = nf90_inquire(ncid, nvariables=nvar)
      if (iret /= NF90_NOERR) call parallel_abort('add_user_attributes: nf90_inquire nvar')

      do varid = 1, nvar
         iret = nf90_inquire_variable(ncid, varid, name=varname)
         if (iret /= NF90_NOERR) call parallel_abort('add_user_attributes: nf90_inquire variable')

         iret = nf90_inq_varid(ncin, trim(varname), ivar)
         if (iret /= NF90_NOERR) cycle ! skipt to next

         call netcdf_copy_attributes(ncin, ivar, ncid, varid)
      end do

      ! The variable crs does not exist in SCHISM output if ics==1, so it needs
      ! to be copied unconditionally
      iret = nf90_inq_varid(ncid, 'crs', varid)
      if (iret /= NF90_NOERR) then
         iret = nf90_inq_varid(ncin, 'crs', ivar)
         if (iret == NF90_NOERR) then
            iret=nf90_inq_dimid(ncid, 'one', dimid)
            if(iret.ne.NF90_NOERR) call parallel_abort('add_user_attributes: inq_dimid one')
            iret=nf90_def_var(ncid,'crs',NF90_INT,(/dimid/),varid)
            if(iret.ne.NF90_NOERR) call parallel_abort('add_user_attributes: def_var crs')
            call netcdf_copy_attributes(ncin, ivar, ncid, varid)

            do varid =1, nvar
               iret = nf90_inquire_variable(ncid, varid, natts=natt)
               if (iret == NF90_NOERR .and. natt > 0) then
                  do iatt = 1, natt
                     iret = nf90_inq_attname(ncid, varid, iatt, attname)
                     if (iret == NF90_NOERR .and. trim(attname) == 'coordinates') then
                        iret = nf90_put_att(ncid, varid, 'grid_mapping', 'crs')
                        exit
                     end if
                  end do
               end if
            enddo
         end if
      end if

      call netcdf_copy_attributes(ncin, NF90_GLOBAL, ncid, NF90_GLOBAL)

   end subroutine add_user_attributes

   subroutine netcdf_copy_attributes(ncin, varid_in, ncout, varid_out)
      implicit none
      integer, intent(in) :: ncin, varid_in, ncout, varid_out

      integer :: iret, natt, attnum, ndim, nvar
      character(len=NF90_MAX_NAME) :: attname
      integer :: xtype, attlen
      integer :: attval_int
      character(len=:), allocatable :: attval_char
      real :: attval_real
      double precision :: attval_double

      if (varid_in == NF90_GLOBAL .and. varid_out == NF90_GLOBAL) then
         iret = nf90_inquire(ncin, ndim, nvar, natt)
         if (iret /= NF90_NOERR) call parallel_abort('netcdf_copy_attributes: inquire_global natts')
      elseif (varid_in /= NF90_GLOBAL .and. varid_out /= NF90_GLOBAL) then
         iret = nf90_inquire_variable(ncin, varid_in, natts=natt)
         if (iret /= NF90_NOERR) call parallel_abort('netcdf_copy_attributes: inquire_variable natts')
      else
         call parallel_abort('netcdf_copy_attributes: cannot mix global and variable attributes')
      end if
      !ierr = nf90_inquire(ncid_src, ndims, nvars, ngatts, unlimdimid)


      do attnum = 1, natt
         iret = nf90_inq_attname(ncin, varid_in, attnum, attname)
         if (iret /= NF90_NOERR) cycle

         iret = nf90_inquire_attribute(ncin, varid_in, attname, xtype=xtype, len=attlen)
         if (iret /= NF90_NOERR) cycle

         select case (xtype)
          case (NF90_CHAR)
            allocate(character(len=attlen) :: attval_char)
            iret = nf90_get_att(ncin, varid_in, attname, attval_char)
            if (iret == NF90_NOERR) iret = nf90_put_att(ncout, varid_out, attname, attval_char)
            deallocate(attval_char)
          case (NF90_INT)
            iret = nf90_get_att(ncin, varid_in, attname, attval_int)
            if (iret == NF90_NOERR) iret = nf90_put_att(ncout, varid_out, attname, attval_int)
          case (NF90_FLOAT)
            iret = nf90_get_att(ncin, varid_in, attname, attval_real)
            if (iret == NF90_NOERR) iret = nf90_put_att(ncout, varid_out, attname, attval_real)
          case (NF90_DOUBLE)
            iret = nf90_get_att(ncin, varid_in, attname, attval_double)
            if (iret == NF90_NOERR) iret = nf90_put_att(ncout, varid_out, attname, attval_double)
          case default
         end select
      end do

   end subroutine netcdf_copy_attributes

   subroutine add_mesh_attributes(ncid, varid, iheader)
      implicit none
      integer, intent(in) :: ncid, varid, iheader

      integer :: iret, ndims, i
      character(len=4)     :: location
      character(len=39)    :: coordinates
      character(len=255)   :: varname, dimname
      integer, allocatable :: dimids(:)

      iret = nf90_inquire_variable(ncid, varid, name=varname, ndims=ndims)
      if(iret.ne.NF90_NOERR) call parallel_abort('add_mesh_attributes: inquire_variable')

      allocate(dimids(ndims))
      iret = nf90_inquire_variable(ncid, varid, dimids=dimids)
      if(iret.ne.NF90_NOERR) call parallel_abort('add_mesh_attributes: inquire_variable')
      do i = 1, ndims
         iret = nf90_inquire_dimension(ncid, dimids(i), name=dimname)
         if (trim(dimname) == 'nSCHISM_hgrid_node') location = 'node'
         if (trim(dimname) == 'nSCHISM_hgrid_face') location = 'face'
         if (trim(dimname) == 'nSCHISM_hgrid_edge') location = 'edge'
      enddo
      deallocate(dimids)

      write(coordinates,'(6A)') 'SCHISM_hgrid_', location, '_x ', &
         'SCHISM_hgrid_', location, '_y'
      iret=nf90_put_att(ncid,varid,'coordinates',trim(coordinates))
      if(iret.ne.NF90_NOERR) call parallel_abort(varname)

      iret=nf90_put_att(ncid,varid,'location',trim(location))
      if(iret.ne.NF90_NOERR) call parallel_abort(varname//'(2)')
      if (ics > 1) then
         iret=nf90_put_att(ncid,varid,'grid_mapping','crs')
         if(iret.ne.NF90_NOERR) call parallel_abort(varname//'(3)')
      endif
      iret=nf90_put_att(ncid,varid,'mesh','SCHISM_hgrid')
      if(iret.ne.NF90_NOERR) call parallel_abort(varname//'(4)')
      ! todo add xtype-dependent fill value
      !iret=nf90_put_att(ncid,varid,'_FillValue',NF90_FILL_FLOAT)
      !if(iret.ne.NF90_NOERR) call parallel_abort(varname)

   end subroutine add_mesh_attributes

   subroutine fill_header_static(iheader,ncid_schism0,itime_id0,node_dim0, &
   &nele_dim0,nedge_dim0,four_dim0,nv_dim0,one_dim0,two_dim0,time_dim0)
      implicit none
      integer, intent(in) :: iheader,ncid_schism0 !header: more elaborate header per UGRID
      integer, intent(out) :: itime_id0,node_dim0,nele_dim0,nedge_dim0, &
      &four_dim0,nv_dim0,one_dim0,two_dim0,time_dim0

      integer :: irec,iret,i,j,k,ih0_id2,ikbp_id2, ivarid,time_dims(1)
      integer :: ix_id2,iy_id2,ih_id2,ixel_id2,iyel_id2,ixsd_id2,iysd_id2,elnode_id2,iside_id2

      !Header
      ! CF recommends that dimensions appear in the order T Z, X/Y, in the CDL description.
      iret=nf90_def_dim(ncid_schism0,'nSCHISM_hgrid_node',np_global,node_dim0)
      iret=nf90_def_dim(ncid_schism0,'nSCHISM_hgrid_face',ne_global,nele_dim0)
      iret=nf90_def_dim(ncid_schism0,'nSCHISM_hgrid_edge',ns_global,nedge_dim0)
      iret=nf90_def_dim(ncid_schism0,'nMaxSCHISM_hgrid_face_nodes',4, four_dim0)
      iret=nf90_def_dim(ncid_schism0,'nSCHISM_vgrid_layers',nvrt,nv_dim0)
      iret=nf90_def_dim(ncid_schism0,'one',1,one_dim0)
      iret=nf90_def_dim(ncid_schism0,'two',2,two_dim0)
      iret=nf90_def_dim(ncid_schism0,'time', NF90_UNLIMITED,time_dim0)

      ! Write the coordinate axis for the time dimension
      time_dims(1)=time_dim0
      iret=nf90_def_var(ncid_schism0,'time',NF90_DOUBLE,time_dims,itime_id0)
      if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: time dim')
      iret=nf90_put_att(ncid_schism0,itime_id0,'i23d',0) !set i23d flag
      iret=nf90_put_att(ncid_schism0,itime_id0,'base_date',start_time)
      iret=nf90_put_att(ncid_schism0,itime_id0,'units',trim(isotimestring))
      iret=nf90_put_att(ncid_schism0,itime_id0,'standard_name','time')
      iret=nf90_put_att(ncid_schism0,itime_id0,'axis','T')
      if (iheader > 0) then
         iret=nf90_put_att(ncid_schism0,itime_id0,'calendar','proleptic_gregorian')

         ! UGRID does not need to be specified as it is included in CF-1.12
         iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'Conventions', 'CF-1.12')
         iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'title', 'SCHISM unstructured grid output')
         iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'institution', 'anonymous')

         iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'source', 'SCHISM')
         iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'references', 'Zhang, Y., Ye, F., Stanev, E.V., Grashorn, S. (2016) Seamless cross-scale modeling with SCHISM, Ocean Modelling, 102, 64-81; http://schism.wiki/')
         iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'history', 'Created ' // trim(iso8601_now()))
         iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'creation_date', trim(iso8601_now()))
         iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'license', 'unspecified')
         iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'comment', 'SCHISM model output file')
         iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'originator', 'anonymous')

         if (ics > 1) then
            iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'crs', 'EPSG:4326') !WGS84
         endif

         ! For OpenDAP retrieval of time axis there should be StartTime/StopTime/StartLatitude etc.
         !iret = nf90_put_att(ncid_schism0, NF90_GLOBAL, 'StartTime', trim(isotimestring))

         ! Metadata that is dimension "one" should come here
         time_dims(1)=one_dim0
         iret=nf90_def_var(ncid_schism0,'minimum_depth',NF90_DOUBLE,time_dims,ih0_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: h0')
         iret=nf90_put_att(ncid_schism0,ih0_id2,'units','m')
         iret=nf90_put_att(ncid_schism0,ih0_id2,'long_name','Minimum depth at which water column is considered wet')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: h0')
      endif !iheader > 0

      if(iheader == 1) then
         ! Do not write the UGRID information for iheader == 2 or iheader == 0

         ! The CF convention requires for unstructured data a dimensionless
         ! field with the cf_role "mesh_topology", with pointers to the node/face/edge information
         iret=nf90_def_var(ncid_schism0,'SCHISM_hgrid',NF90_CHAR,time_dims,ivarid)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid')
         iret=nf90_put_att(ncid_schism0,ivarid,'long_name',"Topology data of 2d unstructured mesh")
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid(2)')
         iret=nf90_put_att(ncid_schism0,ivarid,'topology_dimension',2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid(3)')
         iret=nf90_put_att(ncid_schism0,ivarid,'cf_role',"mesh_topology")
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid(4)')
         iret=nf90_put_att(ncid_schism0,ivarid,'node_coordinates',"SCHISM_hgrid_node_x SCHISM_hgrid_node_y")
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid(5)')
         iret=nf90_put_att(ncid_schism0,ivarid,'edge_coordinates',"SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y")
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid(6)')
         iret=nf90_put_att(ncid_schism0,ivarid,'face_coordinates',"SCHISM_hgrid_face_x SCHISM_hgrid_face_y")
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid(7)')
         iret=nf90_put_att(ncid_schism0,ivarid,'edge_node_connectivity',"SCHISM_hgrid_edge_nodes")
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid(8)')
         iret=nf90_put_att(ncid_schism0,ivarid,'face_node_connectivity',"SCHISM_hgrid_face_nodes")
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid(9)')
         iret=nf90_put_att(ncid_schism0,ivarid,'long_name',"SCHISM unstructured grid topology")
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid(10)')
         iret=nf90_put_att(ncid_schism0,ivarid,'units',"1")
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: SCHISM_hgrid(11)')

         if (ics > 1) then
            iret=nf90_def_var(ncid_schism0,'crs',NF90_INT,time_dims,ivarid)
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: crs')
            iret=nf90_put_att(ncid_schism0,ivarid,'long_name',"Coordinate reference system (CRS) definition")
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: crs(2)')
            iret=nf90_put_att(ncid_schism0,ivarid,'grid_mapping_name',"latitude_longitude")
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: crs(3)')
            iret=nf90_put_att(ncid_schism0,ivarid,'longitude_of_prime_meridian',0.0)
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: crs(4)')
            iret=nf90_put_att(ncid_schism0,ivarid,'semi_major_axis',6378137.0)
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: crs(5)')
            iret=nf90_put_att(ncid_schism0,ivarid,'inverse_flattening',298.257223563)
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: crs(6)')
         endif

         time_dims(1)=node_dim0

         ! x and y coordinates
         iret=nf90_def_var(ncid_schism0,'SCHISM_hgrid_node_x',NF90_DOUBLE,time_dims,ix_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xnd')
         iret=nf90_put_att(ncid_schism0,ix_id2,'axis','X')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xnd(2)')
         iret=nf90_put_att(ncid_schism0,ix_id2,'location','node')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xnd(3)')
         iret=nf90_put_att(ncid_schism0,ix_id2,'mesh','SCHISM_hgrid')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xnd(4)')
         if (ics > 1) then
            iret=nf90_put_att(ncid_schism0,ix_id2,'units','degree_E')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xnd(5)')
            iret=nf90_put_att(ncid_schism0,ix_id2,'standard_name','longitude')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xnd(6)')
         else
            iret=nf90_put_att(ncid_schism0,ix_id2,'units','m')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xnd(7)')
            iret=nf90_put_att(ncid_schism0,ix_id2,'standard_name','projection_x_coordinate')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xnd(8)')
         endif

         iret=nf90_def_var(ncid_schism0,'SCHISM_hgrid_node_y',NF90_DOUBLE,time_dims,iy_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ynd')
         iret=nf90_put_att(ncid_schism0,iy_id2,'axis','Y')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ynd(2)')
         iret=nf90_put_att(ncid_schism0,iy_id2,'location','node')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ynd(3)')
         iret=nf90_put_att(ncid_schism0,iy_id2,'mesh','SCHISM_hgrid')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ynd(4)')
         if (ics > 1) then
            iret=nf90_put_att(ncid_schism0,iy_id2,'units','degree_N')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ynd(5)')
            iret=nf90_put_att(ncid_schism0,iy_id2,'standard_name','latitude')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ynd(6)')
         else
            iret=nf90_put_att(ncid_schism0,iy_id2,'units','m')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ynd(7)')
            iret=nf90_put_att(ncid_schism0,iy_id2,'standard_name','projection_y_coordinate')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ynd(8)')
         endif

         iret=nf90_def_var(ncid_schism0,'depth',NF90_FLOAT,time_dims,ih_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: dp')
         iret=nf90_put_att(ncid_schism0,ih_id2,'standard_name','depth')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: dp(1)')
         iret=nf90_put_att(ncid_schism0,ih_id2,'units','m')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: dp(2)')
         iret=nf90_put_att(ncid_schism0,ih_id2,'axis','Z')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: dp(3)')
         iret=nf90_put_att(ncid_schism0,ih_id2,'positive','down')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: dp(4)')
         call add_mesh_attributes(ncid_schism0,ih_id2, iheader)
         call add_cf_variable_attributes(ncid_schism0, ih_id2)

         iret=nf90_def_var(ncid_schism0,'bottom_index_node',NF90_INT,time_dims,ikbp_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: kbp')
         call add_mesh_attributes(ncid_schism0,ikbp_id2, iheader)
         call add_cf_variable_attributes(ncid_schism0, ikbp_id2)

         ! Switch dimension to elements
         time_dims(1)=nele_dim0
         iret=nf90_def_var(ncid_schism0,'SCHISM_hgrid_face_x',NF90_DOUBLE,time_dims,ixel_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xel')
         iret=nf90_put_att(ncid_schism0,ixel_id2,'axis','X')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xel(2)')
         iret=nf90_put_att(ncid_schism0,ixel_id2,'location','face')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xel(3)')
         iret=nf90_put_att(ncid_schism0,ixel_id2,'mesh','SCHISM_hgrid')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xel(4)')
         if (ics > 1) then
            iret=nf90_put_att(ncid_schism0,ixel_id2,'units','degree_E')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xel(5)')
            iret=nf90_put_att(ncid_schism0,ixel_id2,'standard_name','longitude')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xel(6)')
         else
            iret=nf90_put_att(ncid_schism0,ixel_id2,'units','m')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xel(7)')
            iret=nf90_put_att(ncid_schism0,ixel_id2,'standard_name','projection_x_coordinate')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xel(8)')
         endif

         iret=nf90_def_var(ncid_schism0,'SCHISM_hgrid_face_y',NF90_DOUBLE,time_dims,iyel_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: yel')
         iret=nf90_put_att(ncid_schism0,iyel_id2,'axis','Y')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: yel(2)')
         iret=nf90_put_att(ncid_schism0,iyel_id2,'location','face')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: yel(3)')
         iret=nf90_put_att(ncid_schism0,iyel_id2,'mesh','SCHISM_hgrid')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: yel(4)')
         if (ics > 1) then
            iret=nf90_put_att(ncid_schism0,iyel_id2,'units','degree_N')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: yel(5)')
            iret=nf90_put_att(ncid_schism0,iyel_id2,'standard_name','latitude')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: yel(6)')
         else
            iret=nf90_put_att(ncid_schism0,iyel_id2,'units','m')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: yel(7)')
            iret=nf90_put_att(ncid_schism0,iyel_id2,'standard_name','projection_y_coordinate')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: yel(8)')
         endif

         ! Switch dimension to sides
         time_dims(1)=nedge_dim0
         iret=nf90_def_var(ncid_schism0,'SCHISM_hgrid_edge_x',NF90_DOUBLE,time_dims,ixsd_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xsd')
         iret=nf90_put_att(ncid_schism0,ixsd_id2,'axis','X')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xsd(2)')
         iret=nf90_put_att(ncid_schism0,ixsd_id2,'location','edge')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xsd(3)')
         iret=nf90_put_att(ncid_schism0,ixsd_id2,'mesh','SCHISM_hgrid')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xsd(4)')
         if (ics > 1) then
            iret=nf90_put_att(ncid_schism0,ixsd_id2,'units','degree_E')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xsd(5)')
            iret=nf90_put_att(ncid_schism0,ixsd_id2,'standard_name','longitude')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xsd(6)')
         else
            iret=nf90_put_att(ncid_schism0,ixsd_id2,'units','m')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xsd(7)')
            iret=nf90_put_att(ncid_schism0,ixsd_id2,'standard_name','projection_x_coordinate')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: xsd(8)')
         endif

         iret=nf90_def_var(ncid_schism0,'SCHISM_hgrid_edge_y',NF90_DOUBLE,time_dims,iysd_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ysd')
         iret=nf90_put_att(ncid_schism0,iysd_id2,'axis','Y')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ysd(2)')
         iret=nf90_put_att(ncid_schism0,iysd_id2,'location','edge')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ysd(3)')
         iret=nf90_put_att(ncid_schism0,iysd_id2,'mesh','SCHISM_hgrid')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ysd(4)')
         if (ics > 1) then
            iret=nf90_put_att(ncid_schism0,iysd_id2,'units','degree_N')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ysd(5)')
            iret=nf90_put_att(ncid_schism0,iysd_id2,'standard_name','latitude')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ysd(6)')
         else
            iret=nf90_put_att(ncid_schism0,iysd_id2,'units','m')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ysd(7)')
            iret=nf90_put_att(ncid_schism0,iysd_id2,'standard_name','projection_y_coordinate')
            if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: ysd(8)')
         endif

         var2d_dims(1)=four_dim0
         var2d_dims(2)=nele_dim0
         iret=nf90_def_var(ncid_schism0,'SCHISM_hgrid_face_nodes',NF90_INT,var2d_dims,elnode_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: elnode')
         iret=nf90_put_att(ncid_schism0,elnode_id2,'start_index',1)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: elnode(2)')
         iret=nf90_put_att(ncid_schism0,elnode_id2,'_FillValue',-1)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: elnode(3)')
         iret=nf90_put_att(ncid_schism0,elnode_id2,'cf_role','face_node_connectivity')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: elnode(4)')
         iret=nf90_put_att(ncid_schism0,elnode_id2,'units','1')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: elnode(5)')
         iret=nf90_put_att(ncid_schism0,elnode_id2,'long_name','Face-node connectivity table')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: elnode(6)')

         var2d_dims(1)=two_dim0
         var2d_dims(2)=nedge_dim0
         iret=nf90_def_var(ncid_schism0,'SCHISM_hgrid_edge_nodes',NF90_INT,var2d_dims,iside_id2)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: iside')
         iret=nf90_put_att(ncid_schism0,iside_id2,'start_index',1)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: iside(2)')
         iret=nf90_put_att(ncid_schism0,iside_id2,'_FillValue',-1)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: iside(3)')
         iret=nf90_put_att(ncid_schism0,iside_id2,'cf_role','edge_node_connectivity')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: iside(4)')
         iret=nf90_put_att(ncid_schism0,iside_id2,'units','1')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: iside(5)')
         iret=nf90_put_att(ncid_schism0,iside_id2,'long_name','Edge-node connectivity table')
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: iside(6)')
      endif !iheader=1

      iret=nf90_enddef(ncid_schism0)

      !> @todo write these variables only for iheader == 1
      if(iheader == 1) then
         !Write static info (x,y...)
         iret=nf90_put_var(ncid_schism0,ih0_id2,h0)
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static:put h0')
         iret=nf90_put_var(ncid_schism0,ix_id2,xnd,(/1/),(/np_global/))
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: put node_x')
         iret=nf90_put_var(ncid_schism0,iy_id2,ynd,(/1/),(/np_global/))
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: put node_y')
         iret=nf90_put_var(ncid_schism0,ixel_id2,xel,(/1/),(/ne_global/))
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: put face_x')
         iret=nf90_put_var(ncid_schism0,iyel_id2,yel,(/1/),(/ne_global/))
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: put face_y')
         iret=nf90_put_var(ncid_schism0,ixsd_id2,xsd,(/1/),(/ns_global/))
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: put edge_x')
         iret=nf90_put_var(ncid_schism0,iysd_id2,ysd,(/1/),(/ns_global/))
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: put edge_y')
         iret=nf90_put_var(ncid_schism0,ih_id2,real(dp),(/1/),(/np_global/))
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: put depth')
         iret=nf90_put_var(ncid_schism0,ikbp_id2,kbp00,(/1/),(/np_global/))
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: put bottom_index')
         iret=nf90_put_var(ncid_schism0,elnode_id2,elnode,(/1,1/),(/4,ne_global/))
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: put elnode')
         iret=nf90_put_var(ncid_schism0,iside_id2,isidenode,(/1,1/),(/2,ns_global/))
         if(iret.ne.NF90_NOERR) call parallel_abort('fill_header_static: put sidenode')
      endif !iheader == 1

   end subroutine fill_header_static

   function iso8601_now() result(s)
      ! Returns local time in ISO 8601: YYYY-MM-DDThh:mm:ss±HH:MM

      implicit none

      character(len=32) :: s
      integer :: v(8)              ! (YYYY,MM,DD,zone?,hh,mm,ss,ms)
      character(len=5) :: zone_str ! minutes offset from UTC, signed
      integer :: zmin, zh, zm
      character(len=1) :: sgn

      call date_and_time(values=v, zone=zone_str)
      read(zone_str,'(I5)') zmin
      sgn = '+'
      if (zmin < 0) sgn = '-'

      zh = abs(zmin)/60
      zm = mod(abs(zmin),60)

      write(s,'(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,A1,I2.2,":",I2.2)') &
         v(1), v(2), v(3), v(5), v(6), v(7), sgn, zh, zm
   end function iso8601_now

end module scribe_io
