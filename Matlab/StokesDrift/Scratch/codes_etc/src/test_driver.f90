program test_driver
use etc,         only: logfile,istep
use mpi_params
implicit none
integer, parameter  :: mmax=8,nmax=8,vlen=1024

integer             :: ncid,info
integer(kind=8)     :: forward_plan,inverse_plan 
integer             :: rank_trans,howmany_rank
integer             :: m,n,nb2,n_trans,n_loop
integer             :: fft_kind,myflag
integer             :: is_trans,is_loop
integer             :: os_trans,os_loop

integer             :: nx,ny,nz,nparticles
real(kind=8)        :: dt,t_start,t_end,Lx,Ly,Lz

integer             :: i,j,incx,lda,ipiv(nmax)
real(kind=8)        :: da,x(vlen)

real(kind=8), allocatable  :: in(:,:,:),out(:,:,:),A(:,:)
character(len=80)          :: runlabel

include 'fftw3.f'     
include 'netcdf.inc'
include 'mpif.h'

logfile='test_driver_log'
istep=0


if(myid==0) write(0,*) ' ................   ============================================ '
if(myid==0) write(0,*) ' ................   Saludos de la FlowSolve installacion prueba. '
if(myid==0) write(0,*) ' ................   ============================================ '
if(myid==0) write(0,*) ' ................   '

!!==============================================================================
!!  MPI
!!==============================================================================
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   Prueba de la biblioteca MPI .....'
 if(myid==0) write(0,*) ' ................   '
 comm=MPI_COMM_WORLD
 call mpi_init(ierr)
 call mpi_comm_rank(comm,myid,ierr)
 call mpi_comm_size(comm,numprocs,ierr)

 if( myid .eq. 0 .and. ierr==MPI_SUCCESS ) then
  write(0,*) ' ................        Soy MPI y estoy vivo! '
  write(0,*) ' ................        Tengo ',numprocs,' processero(s).'
 elseif( myid .eq. 0 .and. ierr .ne. MPI_SUCCESS ) then
  write(0,*) ' ................        Houston, tenemos un problema con MPI... '
  write(0,*) ' ................         ierr   '    ,ierr
  write(0,*) ' ................         myid   '    ,myid
  write(0,*) ' ................         numprocs   ',numprocs
 endif
 if(myid==0) write(0,*) ' ................   '  

!!==============================================================================
!!  NETCDF
!!==============================================================================
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   Prueba de la biblioteca NETCDF .....'
 if(myid==0) write(0,*) ' ................   '
 ierr=NF_CREATE('install_test.nc',NF_CLOBBER,ncid)
 if(ierr.eq.NF_NOERR) then
  if(myid==0) write(0,*) ' ................        voy a crear un archivo netcdf: install_test.nc'
 else
  if(myid==0) write(0,*) ' ................        hubo un problema creando el archivo netcdf: ',ierr
 endif

 ierr=NF_ENDDEF(ncid)
 if(ierr.eq.NF_NOERR) then
  if(myid==0) write(0,*) ' ................        salida del modo define bien'
 else
  if(myid==0) write(0,*) ' ................        hubo un error mientras saliendo el modo define'
 endif

 ierr=NF_CLOSE(ncid)
 if(ierr.ne.NF_NOERR .and. myid==0) &
  write(0,*) ' ................        hubo problemas cerrando el archivo netcdf'
 if(myid==0) write(0,*) ' ................   '
 
!!==============================================================================
!!  BLAS
!!==============================================================================
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   Prueba de  DBLAS .....'
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................        dscal 1024 identity vector'
 incx=1
 x(:)=1.d0
 da=4.d0*atan(1.d0)
 call dscal(VLEN,da,x,incx)
 if(myid==0) write(0,*) ' ................        el resultado puede ser pi: ',x(1)
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   '


!!==============================================================================
!!  LAPACK
!!==============================================================================
 lda=mmax
 m=lda
 n=nmax
 call RANDOM_SEED  !! f90 intrinsic
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   Prueba de  LAPACK .....'
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................        factor random dense matrix '
 if(myid==0) write(0,*) ' ................        ',lda,'   x',nmax

 allocate( A(lda,nmax) )
 do i=1,lda
  do j=1,nmax
   call RANDOM_NUMBER( A(i,j) )
  enddo
 enddo
 
 call dgetrf(m,n,A,lda,ipiv,info)
 if(info==0) then
  if(myid==0) write(0,*) ' ................        info=0 ==> Felicitaciones: la prueba fue un exito '
 elseif(info .gt. 0 ) then
  if(myid==0) write(0,*) ' ................        random matrix found to be singular: ',info
 else
  if(myid==0) write(0,*) ' ................                    m, n, lda: ',m,n,lda
  if(myid==0) write(0,*) ' ................                    ipiv: ',ipiv
  if(myid==0) write(0,*) ' ................                       A: ',A
 endif
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   '
 deallocate (A)

!!==============================================================================
!!  FFTW3
!!==============================================================================
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   Prueba de la biblioteca FFTW3 .....'
 if(myid==0) write(0,*) ' ................   '

 rank_trans=1          !! Transforms work on 1 dimension at a time.
 howmany_rank=1        !! Loop over the other 2 dimensions treated as 1d in memory
 n = 32                !! Data size in transform direction, measured in real words.
 nb2 = n/2
 n_trans=n             !! length of transform in real words     
 is_trans=1            !! stride of 1d transform --> in x
 os_trans=is_trans
 !! describe the loops over dims 2,3
 n_loop=n*n            !! 
 is_loop=n             !! stride between transforms
 os_loop=is_loop
 myflag=FFTW_PATIENT
 fft_kind=FFTW_R2HC
 allocate( in(n,n,n),out(n,n,n) )

 call dfftw_plan_guru_r2r(forward_plan,           &
                          rank_trans,n_trans,     &
                          is_trans,os_trans,      &
                          howmany_rank,n_loop,    &
                          is_loop,os_loop,        &
                          in,out,fft_kind,myflag)
 if(myid==0) write(0,*) ' ................        creando dfftw_plan_guru_r2r ok ..'
 in(:,:,:) = 0.0
 out(:,:,:)= 0.0
 call dfftw_execute_r2r(forward_plan,in,out)
 if(myid==0) write(0,*) ' ................        realizando dfftw_execute_r2r ok ..'
 if(myid==0) write(0,*) ' ................   '
 deallocate(in,out)
 
 
!!==============================================================================
!! Shut down MPI
!!==============================================================================
 call mpi_finalize(ierr)
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   Shutting down MPI ..... '
 if(myid==0) write(0,*) ' ................        MPI esta muerto ahora. ',ierr
 if(myid==0) write(0,*) ' ................   '

 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   '
 if(myid==0) write(0,*) ' ................   =========================================== '
 if(myid==0) write(0,*) ' ................   La installacion me parece bien. Adios. '
 if(myid==0) write(0,*) ' ................   =========================================== '
end program test_driver
