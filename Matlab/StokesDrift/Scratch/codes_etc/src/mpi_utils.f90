subroutine start_mpi(myid,numprocs,logfile,comm) 
  implicit none
  integer             :: ierr,comm,myid,numprocs
  character(len=80)   :: logfile
  
  include 'mpif.h'
 
  comm=MPI_COMM_WORLD
  call mpi_init(ierr)
   if( ierr .ne. MPI_SUCCESS ) stop 'mpi_init error in start_mpi'
   
  call mpi_comm_rank(comm,myid,ierr)
   if( ierr .ne. MPI_SUCCESS ) stop 'mpi_comm_rank error in start_mpi'
   
  call mpi_comm_size(comm,numprocs,ierr)
   if( ierr .ne. MPI_SUCCESS ) stop 'mpi_comm_size error in start_mpi'
   
      
  if( myid .eq. 0 .and. ierr==MPI_SUCCESS ) then
    write(0,*)'----------------------------------------------------'
    if(numprocs==1) write(0,*) ' MPI initialized, np=1; ==> SINGLE PROCESSOR RUN'
    if(numprocs>1)  write(0,*) ' MPI initialized with ',numprocs,' live processors'
    write(0,*)'----------------------------------------------------'
    write(0,*) ' '
    open(1,file=logfile,position='rewind')
      write(1,*) ' '
      write(1,*) ' MPI initialized with ',numprocs,' live processors'
      write(1,*) '             comm =   ',comm
      write(1,*) ' MPI error flag   =   ',ierr
    close(1)
   endif 
   call mpi_barrier(comm,ierr)
 end subroutine start_mpi



 subroutine quit_mpi
   implicit none
   integer   :: jerr
   include 'mpif.h'
   call mpi_finalize(jerr)
   stop
  end subroutine quit_mpi

