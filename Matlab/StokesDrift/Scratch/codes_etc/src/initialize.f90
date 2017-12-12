subroutine initialize 
 use mpi_params,              only: numprocs,myid,comm,ierr
 use decomposition_params,    only: np
 use etc,                     only: logfile,message
 use intermediate_variables,  only: tmpY
 implicit none 
 

!!===========================================================
!! Initialize mpi
!!===========================================================
 call start_mpi(myid,np,logfile,comm)
 numprocs=np

!!===========================================================
!! Initialize random number generator
!!===========================================================
 call random_seed  !! enable use of F90 intrinsic RV functions


 
!!===========================================================
!! Read in user data from input directory 
!!===========================================================
 if(myid==0) then
  message = '...  Calling ReadUserData...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif
 call ReadUserData
 


!!===========================================================
!! Set params for data decomposition
!!===========================================================
 if(myid==0) then
  message = '...  Calling Decomposition ...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif
 call Decomposition2D
 
 
 
!!===========================================================
!! Misc. preliminaries
!!===========================================================
 if(myid==0) then
  message =  '...  Calling PreliminaryTasks...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif
 call PreliminaryTasks



!!===========================================================
!! Setup spatial domain
!!===========================================================
 if(myid==0) then
  message = '...  Calling SetupDomain ...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif
 call SetupDomain
 
 
 
!!===========================================================
!! Setup differentiation matrices/fftw plans/wavenumbers etc.
!!===========================================================
 if(myid==0) then
  message = '...  Calling SetupDerivs ...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif
 call SetupDerivs


!!===========================================================
!! Allocate memory for dependent variables
!!===========================================================
 if(myid==0) then
  message = '...  Calling AllocateDependentVariables ...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif
 call AllocateDependentVariables 
 

!!===========================================================
!! Allocate memory for intermediate variables
!!=========================================================== 
 if(myid==0) then
  message = '...  Calling AllocateIntermediateVariables ...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif
 call AllocateIntermediateVariables 

!!===========================================================
!! Quick test of transform_xy
!!===========================================================
 !if(myid==0) then
  !message = '...  Calling test_transform_xy ...'
  !write(0,*) message
  !call LogMessage(message,logfile)
 !endif
 !call test_transform_xy(tmpY(1,1,1,1),tmpY(1,1,1,2))

 
!!===========================================================
!! Setup Diffusion operators 
!!===========================================================
 if(myid==0) then
  message = '...  Calling SetupDiffusion ...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif
 call SetupDiffusion 
 

 
!!===========================================================
!! Set the initial conditions
!!===========================================================
 if(myid==0) then
  message = '...  Calling InitialConditions ...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif 
 call InitialConditions 

 
!!==============================================================
!! Setup stuff for calls to cubic spline interpolation routines
!!==============================================================
 if(myid==0) then
  message = '...  Calling SetupInterpolation ...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif
 call SetupInterpolation
 
 


!!===========================================================
!! Set the immersed boundary
!!===========================================================
 if(myid==0) then
  message = '...  Calling SetupImmersedBoundary ...'
  write(0,*) message
  call LogMessage(message,logfile)
 endif 
 call SetupImmersedBoundary
 
 
 

                         
     
 if(myid==0) then
   message = ' -----> routine Initialize exiting normally <---------- '
   write(0,*) ' '
   write(0,*) '----------------------------------------------------'
   write(0,*) message
   write(0,*) '----------------------------------------------------'
   write(0,*) ' '
   call LogMessage(message,logfile)
  endif 



 call mpi_barrier(comm,ierr)

return
end subroutine initialize















 







