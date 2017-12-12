      module comm

      include 'mpif.h'

#ifdef DOUBLE_PREC
       integer, parameter,public :: mytype=8
       integer, parameter,public:: mpireal = MPI_DOUBLE_PRECISION
       integer,parameter,public:: mpicomplex = MPI_DOUBLE_COMPLEX
#else
       integer, parameter,public :: mytype=4
       integer, parameter,public:: mpireal = MPI_REAL
       integer,parameter,public:: mpicomplex = MPI_COMPLEX
#endif
       integer, parameter :: MAX_NS = 50
       integer, save :: ls
       type bound
          integer, dimension(:,:,:), pointer :: R
       end type bound
       type count
          integer, dimension(:,:), pointer :: RSndCnts
          integer, dimension(:,:), pointer :: RRcvCnts
          integer, dimension(:,:), pointer :: RSndStrt
          integer, dimension(:,:), pointer :: RRcvStrt
       end type count

! index bounds for data partition
       type(bound), dimension(MAX_NS),save :: st,sz,en

! message counts and displacements for mpi_alltoallv
       type(count), dimension(3,3,MAX_NS) :: Exch

      integer, public :: padd

      integer,save :: ns,GAR_DIM(3,MAX_NS),proc_grid(3,MAX_NS),gpid(2,MAX_NS)
      integer,save :: numtasks
!      integer,save :: ipid,jpid,taskid
!     integer,save :: iistart,iiend,iisize,jjstart,jjsize,jjend
!     integer,save ::jistart,kjstart,jisize,kjsize,jiend,kjend

      integer, parameter :: gblock = 1

! mpi process info
!
      logical,save :: mpi_set(MAX_NS)=.false.
      integer, save :: mpi_comm_row(MAX_NS), mpi_comm_col(MAX_NS)


      end
