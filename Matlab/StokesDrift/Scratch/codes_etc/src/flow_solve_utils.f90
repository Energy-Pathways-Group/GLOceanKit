!!===========================================================
!!===========================================================
subroutine toggle
 use etc
 use mpi_params,             only: myid
 use dimensional_scales,     only: time_scale,        &
                                   velocity_scale,    &
                                   scalar_scale
 use independent_variables,  only: dt,t_secs
 use methods_params,         only: AB_order,AM_order
 use dependent_variables
 use timing
 implicit none
 
 integer                        :: itmp,iprint

 iprint=25
  
 if( myid .eq. 0 .and. mod(istep,iprint)==0 ) then
  write(0,*) '..........................................................'
  write(0,*) '............................     Toggling, istep  ',istep
  write(0,*) '..........................................................'
 endif
  
 t_secs = t_secs + dt*time_scale
 
!!Toggle 2-cycle pointer for AM solves (N,NM1).
 if( AM_order .eq. 4 ) then
  itmp=NM1
  NM1=N  
  N=itmp
  N_oldest=NM1
 elseif( AM_order .eq. 3 ) then
  N=1
  N_oldest=N
 endif
 
!!Toggle M cycle pointers for AB rhs fields
 if( AB_order .eq. 4 ) then
  itmp=MM3
  MM3=MM2
  MM2=MM1
  MM1=MM0
  MM0=itmp
  M_oldest=MM3
 elseif( AB_order .eq. 3 ) then
  itmp=MM2
  MM2=MM1
  MM1=MM0
  MM0=itmp
  M_oldest=MM2
 elseif( AB_order .eq. 2 ) then
  itmp=MM1
  MM1=MM0
  MM0=itmp
  M_oldest=MM1
 endif
 
!! Counters for Particle Trajectories
!! (always 4th order AB)
  itmp=PM3
  PM3=PM2
  PM2=PM1
  PM1=PM0
  PM0=itmp
  P_oldest=PM3

!! reset immersed boundary data if necessary 
 call SetupImmersedBoundary
   
 t_end_time_step = mpi_wtime()
 t_time_step = t_end_time_step - t_start_time_step
 t_total = t_total + t_time_step
 
 if( mod(istep,iprint)==0 ) call cfl
           
 istep=istep+1
 if(istep > iend) then
  if(myid==0) then
   write(0,*) ' time for last time step:   ',t_time_step
   write(0,*) ' total time (excluding io): ',t_total
  endif
  call finalize
 endif
end subroutine toggle
!!===========================================================
!!===========================================================


!!===========================================================
!!===========================================================
 subroutine LogMessage(msg,logfile)
  character(len=80)  :: msg,logfile
  open(1,file=logfile,position='append')
   write(1,*) msg
  close(1)
 end subroutine LogMessage
!!===========================================================
!!===========================================================



!!===========================================================
!!===========================================================
subroutine finalize 
 use etc,            only: istep,istart
 use mpi_params,     only: myid
 implicit none
 if(myid==0) &
 write(0,*) 'smooth exit after ',istep-istart-1,' time steps'
 call quit_mpi
 stop
end subroutine finalize
!!===========================================================
!!===========================================================


subroutine cfl
 use etc
 use mpi_params
 use dimensional_scales,     only: time_scale,        &
                                   velocity_scale,    &
                                   length_scale
 use independent_variables,  only: dt,x,y,z,nx,ny,nz
 use dependent_variables
 implicit none
 include 'mpif.h' 
 
 real(kind=8)      :: dx,dy,dz,dt_in_secs
 real(kind=8)      :: cfl_x,cfl_y,cfl_z
 real(kind=8)      :: umax,vmax,wmax,global_max
 integer           :: icount,root_pid
  
 icount=1         ! compute max of single value
 root_pid=0       ! return mpi_reduce results to pid=0
 dt_in_secs = dt*time_scale
 
 !------------------------------
 ! local value of max speeds
 !------------------------------
  umax = MAXVAL( ABS(u) )*velocity_scale
  vmax = MAXVAL( ABS(v) )*velocity_scale
  wmax = MAXVAL( ABS(w) )*velocity_scale
  
 !------------------------------
 ! local max CFL values
 !------------------------------
 if( nx > 1 ) then
  dx = (x(2)-x(1))*length_scale   ! [m]
  cfl_x = umax*dt_in_secs/dx
 else
  cfl_x=0.d0
 endif
 
 if( ny > 1 ) then
  dy = (y(2)-y(1))*length_scale   ! [m]
  cfl_y = vmax*dt_in_secs/dy
 else
  cfl_y = 0.d0             
 endif   
 
 dz = (z(2)-z(1))*length_scale    ! [m]
 cfl_z = wmax*dt_in_secs/dz
 
 !----------------------------------------------------
 !  find the maximum values across all the processors
 !----------------------------------------------------
 
 call MPI_REDUCE(cfl_x,global_max,icount,MPI_DOUBLE_PRECISION,MPI_MAX,root_pid,comm,ierr)
 if( myid==0 ) cfl_x = global_max
 
 call MPI_REDUCE(cfl_y,global_max,icount,MPI_DOUBLE_PRECISION,MPI_MAX,root_pid,comm,ierr)
 if( myid==0 ) cfl_y = global_max
 
 call MPI_REDUCE(cfl_z,global_max,icount,MPI_DOUBLE_PRECISION,MPI_MAX,root_pid,comm,ierr)
 if( myid==0 ) cfl_z = global_max
 
 if( myid==0 ) then
  write(0,*) '...............     ADVECTIVE CFL ratios '
  write(0,*) '...............              x :  ',cfl_x
  write(0,*) '...............              y :  ',cfl_y
  write(0,*) '...............              z :  ',cfl_z
  if( cfl_x > 1 ) stop 'CFL_x > 1 incipient instability '
  if( cfl_y > 1 ) stop 'CFL_y > 1 incipient instability '
  if( cfl_z > 1 ) stop 'CFL_z > 1 incipient instability '
  if( istep > 0 .and. cfl_x==0. .and. cfl_y==0. .and. cfl_z==0. )  &
   stop 'all CFL numbers zero: probably post-instability '
 endif

end subroutine cfl 




!!===========================================================
!!=========================================================== 
subroutine YBdotproduct(ans,sign,u1,u2,u3,v1,v2,v3)
!---------------------------------------------------
! Assume input and output arrays in YBLOCK format
!---------------------------------------------------
 use mpi_params,            only: myid
 use decomposition_params 
 implicit none
 integer                   :: i,j,k,n,sign
 real(kind=8)              :: ans(array_size(IDIM,YBLOCK,myid),     &
                                  array_size(JDIM,YBLOCK,myid),     &
                                  array_size(KDIM,YBLOCK,myid) )
                                  
 real(kind=8)              ::  u1(array_size(IDIM,YBLOCK,myid),     &
                                  array_size(JDIM,YBLOCK,myid),     &
                                  array_size(KDIM,YBLOCK,myid) )
                                  
 real(kind=8)              ::  u2(array_size(IDIM,YBLOCK,myid),     &
                                  array_size(JDIM,YBLOCK,myid),     &
                                  array_size(KDIM,YBLOCK,myid) )
                                  
 real(kind=8)              ::  u3(array_size(IDIM,YBLOCK,myid),     &
                                  array_size(JDIM,YBLOCK,myid),     &
                                  array_size(KDIM,YBLOCK,myid) )
                                  
 real(kind=8)              ::  v1(array_size(IDIM,YBLOCK,myid),     &
                                  array_size(JDIM,YBLOCK,myid),     &
                                  array_size(KDIM,YBLOCK,myid) )
                                  
 real(kind=8)              ::  v2(array_size(IDIM,YBLOCK,myid),     &
                                  array_size(JDIM,YBLOCK,myid),     &
                                  array_size(KDIM,YBLOCK,myid) )
                                  
 real(kind=8)              ::  v3(array_size(IDIM,YBLOCK,myid),     &
                                  array_size(JDIM,YBLOCK,myid),     &
                                  array_size(KDIM,YBLOCK,myid) )


 n = array_size(JDIM,YBLOCK,myid)  &
    *array_size(KDIM,YBLOCK,myid) 

if( sign == -1 ) then

 ans(:,:,:) = sign*(u1(:,:,:)*v1(:,:,:) + u2(:,:,:)*v2(:,:,:) + u3(:,:,:)*v3(:,:,:))

else

 ans(:,:,:) = (u1(:,:,:)*v1(:,:,:) + u2(:,:,:)*v2(:,:,:) + u3(:,:,:)*v3(:,:,:))

endif
 
 return
end subroutine YBdotproduct
!!===========================================================

!!===========================================================
!!=========================================================== 
subroutine axpy(alpha,x,y)
!---------------------------------------------------
! Assume input and output arrays in YBLOCK format
!---------------------------------------------------
 use mpi_params,            only: myid
 use decomposition_params 
 implicit none
 integer                   :: n,i,j,k
 integer,parameter         :: inc=1
 real(kind=8)              :: alpha
 real(kind=8)              :: x(array_size(IDIM,YBLOCK,myid),     &
                                array_size(JDIM,YBLOCK,myid),     &
                                array_size(KDIM,YBLOCK,myid) )
                                  
 real(kind=8)              :: y(array_size(IDIM,YBLOCK,myid),     &
                                array_size(JDIM,YBLOCK,myid),     &
                                array_size(KDIM,YBLOCK,myid) )
                                    
 n = array_size(JDIM,YBLOCK,myid)  &
    *array_size(KDIM,YBLOCK,myid)  
 
!--------------------------------
!    y <= alpha*x + y
!--------------------------------

 y(:,:,:) = alpha*x(:,:,:) + y(:,:,:)


end subroutine axpy
!!===========================================================
!!===========================================================



!!===========================================================
!!=========================================================== 
subroutine scale(alpha,y)
!---------------------------------------------------
! Assume input and output arrays in YBLOCK format
!---------------------------------------------------
 use mpi_params,            only: myid
 use decomposition_params
 implicit none
 integer                   :: n,i,j,k
 integer,parameter         :: inc=1
 real(kind=8)              :: alpha
 real(kind=8)              :: y(array_size(IDIM,YBLOCK,myid),     &
                                array_size(JDIM,YBLOCK,myid),     &
                                array_size(KDIM,YBLOCK,myid) )
 
 n = array_size(JDIM,YBLOCK,myid)  &
    *array_size(KDIM,YBLOCK,myid)
    
!--------------------------------
!    y <= alpha*y
!--------------------------------

 y(:,:,:) = alpha*y(:,:,:)


end subroutine scale
!!===========================================================
!!===========================================================


!!===========================================================
!!=========================================================== 
subroutine copy(x,y)
!---------------------------------------------------
! Assume input and output arrays in YBLOCK format
!---------------------------------------------------
 use mpi_params,            only: myid
 use decomposition_params 
 implicit none
 integer                   :: i,j,k,n
 integer,parameter         :: inc=1
 real(kind=8)              :: x(array_size(IDIM,YBLOCK,myid),     &
                                array_size(JDIM,YBLOCK,myid),     &
                                array_size(KDIM,YBLOCK,myid) )
                                  
 real(kind=8)              :: y(array_size(IDIM,YBLOCK,myid),     &
                                array_size(JDIM,YBLOCK,myid),     &
                                array_size(KDIM,YBLOCK,myid) )
                                    
 n = array_size(JDIM,YBLOCK,myid)  &
    *array_size(KDIM,YBLOCK,myid)  
 
!--------------------------------
!    y <= x
!--------------------------------

 y = x


end subroutine copy





function delta(x,eps)
 implicit none
 real(kind=8) :: delta,x,eps,recip
 
 !!=================================
 !! delta = (1/2eps) sech^2(x/eps)
 !!  (a) in lim eps-->0, delta-->inf
 !!  (b) int delta(x;eps) dx = 1
 !!=================================
 recip = 2.d0*eps*cosh(x/eps)**2
 delta = 1.d0/recip
                                 
 return
end function delta  



function step(x,x0,eps)
 implicit none
 real(kind=8) :: step,x,x0,eps
 
 !!=================================
 !! step = (1/2){ tanh((x-x0)/eps) + 1. }
 !!  [H] = 1. 
 !!  d/dx(H) = delta 
 !!=================================
 step = (0.5d0)*( 1.d0 + tanh((x-x0)/eps)  )
                                 
 return
end function step 



!*****************************************************************
function acosh2 ( x )
!! ACOSH2 returns the inverse hyperbolic cosine of a number.
!
!  Discussion:
!
!    Since a library function ACOSH may be available on some systems,
!    this routine is named ACOSH2 to avoid name conflicts.
!
!    One formula is:
!
!      ACOSH2 = LOG ( X + SQRT ( X**2 - 1.0 ) )
!
!    but this formula suffers from roundoff and overflow problems.
!    The formula used here was recommended by W Kahan, as discussed
!    by Moler.
!
!    Applying the inverse function
!
!      Y = ACOSH2(X)
!
!    implies that
!
!      X = COSH(Y) = 0.5 * ( EXP(Y) + EXP(-Y) ).
!
!    For every X greater than or equal to 1, there are two possible
!    choices Y such that X = COSH(Y), differing only in sign.  It
!    is usual to resolve this choice by taking the value of ACOSH2(X)
!    to be nonnegative.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cleve Moler,
!    Trigonometry is a Complex Subject,
!    MATLAB News and Notes,
!    Summer 1998.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose inverse hyperbolic 
!    cosine is desired.  X should be greater than or equal to 1.
!
!    Output, real ( kind = 8 ) ACOSH2, the inverse hyperbolic cosine of 
!    X.  The principal value (that is, the positive value of the two ) 
!    is returned.
!
  implicit none

  real ( kind = 8 ) acosh2
  real ( kind = 8 ) x

  if ( x < 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ACOSH2 - Fatal error!'
    write ( *, '(a)' ) '  Argument X must satisfy 1 <= X.'
    write ( *, '(a,g14.6)' ) '  The input X = ', x
    stop
  end if

  acosh2 = 2.0D+00 * log ( &
    sqrt ( 0.5D+00 * ( x + 1.0D+00 ) ) + sqrt ( 0.5D+00 * ( x - 1.0D+00 ) ) )

  return
end


integer function factorial(k)
integer::ans,n,k
 ans=1
 n=k
 do while (n.GT.1)
  ans=ans*n
  n=n-1
 enddo
 factorial=ans
 return
end function factorial
!===================================================    
!===================================================

subroutine dinit(n,alpha,y)
 implicit none
 integer       :: n,i
 real(kind=8)  :: y(n),alpha

!$omp parallel default(shared) private(i)
!$omp do 
 do i=1,n
  y(i)=alpha
 enddo
!$omp end do
!$omp end parallel

end subroutine dinit

subroutine get_2d_indices(k,nx,i,j)
!-------------------------------------------------------------
!
!  do j=1,ny             do k=1,nx*ny
!   do i=1,nx             call get_2d_indices(k,nx,i,j)
!    work(i,j)    ==>     work(i,j)
!   enddo                enddo
!  enddo
!
!  NB  ny not needed except for testing
!-------------------------------------------------------------

 implicit none
 integer        :: i,j,k,nx

  j = floor( (k-1.)/nx ) + 1
  i = k - (j-1)*nx 
 
end subroutine get_2d_indices


subroutine get_3d_indices(ii,nx,ny,i,j,k)
!-------------------------------------------------------------
! do k=1,nz
!  do j=1,ny             do ii=1,nx*ny*nz
!   do i=1,nx             call get_3d_indices(ii,nx,ny,i,j,k)
!    work(i,j,k)    ==>   work(i,j,k)
!   enddo                enddo
!  enddo
! enddo
!
!  NB  ny not needed except for testing
!-------------------------------------------------------------

 implicit none
 integer        :: ii,i,j,k,nx,ny,i2d

  k = floor( (ii-1.)/(nx*ny) ) + 1

  i2d = ii - (k-1)*(nx*ny)

  j = floor( (i2d-1.)/nx ) + 1
  i = i2d - (j-1)*nx 
 
end subroutine get_3d_indices


