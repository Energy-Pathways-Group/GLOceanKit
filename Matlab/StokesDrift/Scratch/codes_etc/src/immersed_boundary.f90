subroutine immersed_bdry
 use decomposition_params
 use mpi_params,              only: myid
 use methods_params
 use immersed_boundary 
 use dependent_variables
 use independent_variables,  only: z
 use dimensional_scales,     only: velocity_scale
 
 implicit none
 integer                         :: i,j,k,kg
 
 
 if( .not. do_immersed_boundary ) return
  
 !-------------------------------------------------------
 !  impose Dirichlet conditions at the immersed boundary
 !-------------------------------------------------------
  if( no_slip_no_flux_IB ) then
  
   call dirichlet_blend(u)           ! generic case: no slip IBs, u=v=w=0
   call dirichlet_blend(v)
   call dirichlet_blend(w)
   
  elseif( stokes2 ) then
  
   sigma_n = z(2)-z(1)            ! flat IB surface
   call dirichlet_blend(u)
   call dirichlet_blend(w)
   call stokes2_dirichlet_blend   !  u=w=0, v(t) specified on flat IB surface
   
  elseif( hc_IB_test ) then
  
   call dirichlet_blend(w)        ! w=0
   call hc_dirichlet_blend(s1)    ! s1(y)  specified on flat IB surface
   
  elseif( moving_cylinder ) then
  
   v = v - (U_cyl)/velocity_scale
   w = w - (W_cyl)/velocity_scale
   call dirichlet_blend(u)           ! generic case: no slip IBs, u=v=w=0
   call dirichlet_blend(v)
   call dirichlet_blend(w)
   v = v + (U_cyl)/velocity_scale
   w = w + (W_cyl)/velocity_scale
  endif

  
 !-------------------------------------------------------
 !  impose Neumann conditions at the immersed boundary
 !-------------------------------------------------------
  if( no_slip_no_flux_IB ) then
  
   !------------------------------------------------------------------
   !  add in the background gradient of s1 before imposing no flux bc
   !------------------------------------------------------------------
   do j=1,array_size(IDIM,YBLOCK,myid)
    do i=1,array_size(JDIM,YBLOCK,myid)
     do k=1,array_size(KDIM,YBLOCK,myid)
      kg = global_z_indices(START,YBLOCK,myid) + k - 1
      s1(j,i,k) = s1(j,i,k) + s1_bar(kg,1)
     enddo
    enddo
   enddo
   
   !------------------------------------------------------------------
   !  apply bc to s1_total
   !------------------------------------------------------------------
   call homogeneous_neumann_blend(s1)
   if( do_second_scalar ) then
    call homogeneous_neumann_blend(s2)
   endif
   
   !------------------------------------------------------------------
   !  subtract the background gradient of s1 back out
   !------------------------------------------------------------------
   do j=1,array_size(IDIM,YBLOCK,myid)
    do i=1,array_size(JDIM,YBLOCK,myid)
     do k=1,array_size(KDIM,YBLOCK,myid)
      kg = global_z_indices(START,YBLOCK,myid) + k - 1
      s1(j,i,k) = s1(j,i,k) - s1_bar(kg,1)
     enddo
    enddo
   enddo
  
  elseif( stokes2 ) then
  
   call homogeneous_neumann_blend(s1)
  
  elseif( hc_IB_test ) then
   
   call homogeneous_neumann_blend(u)    !   du/dz=0 on flat IB surface
   call homogeneous_neumann_blend(v)    !   dv/dz=0 on flat IB surface
   
  endif

return
end subroutine immersed_bdry



subroutine dirichlet_blend(phi)
 use mpi_params,            only: myid
 use immersed_boundary,     only: dist,sigma_n
 use decomposition_params
 implicit none
 real(kind=8)                  :: phi( array_size(IDIM,YBLOCK,myid),   &
                                      array_size(JDIM,YBLOCK,myid),    &
                                      array_size(KDIM,YBLOCK,myid)  )
 integer                       :: i,j,k
 real(kind=8)                  :: sigma,d,g,g1,g2
    
 sigma = 0.75*sigma_n   
  
  do k=1,array_size(KDIM,YBLOCK,myid)
   do i=1,array_size(JDIM,YBLOCK,myid)
    do j=1,array_size(IDIM,YBLOCK,myid)    
           
     !------------------------------------------------------------
     ! force odd symmetry and then smoothly clip inside the solid
     !------------------------------------------------------------
     d = dist(j,i,k)
     g1 = tanh( d/sigma )          ! old: g1 = 1.d0 - exp(-(d/sigma)**2)
     phi(j,i,k) = phi(j,i,k)*g1
     
     if( d < 0 ) then              ! taper result inside solid region
      g2 = exp(-(d/(1*sigma))**4)    
      phi(j,i,k) = phi(j,i,k)*g2
     endif
          
   
    enddo
   enddo
  enddo
    
return
end subroutine dirichlet_blend 




subroutine homogeneous_neumann_blend(phi)
 use immersed_boundary               
 use independent_variables,           only: x,y,z,Lx,Ly,Lz
 use dependent_variables,             only: u,v,w,s1,s2
 use mpi_params,                      only: myid,numprocs,comm,ierr
 use dimensional_scales,              only: length_scale
 use decomposition_params

                                      
 implicit none
 include 'mpif.h'
 real(kind=8)                  ::  phi( array_size(IDIM,YBLOCK,myid),   &
                                        array_size(JDIM,YBLOCK,myid),   &
                                        array_size(KDIM,YBLOCK,myid)  )
 
 integer                                 :: i,j,k,kg,N,m
 real(kind=8)                            :: sigma,d,xi,yi,zi,val
    
 !--------------------------------------------------
 ! spatial scale for taper function  in solid region
 !--------------------------------------------------
 sigma = 2*sigma_n   ! was 6

 !--------------------------------------------
 ! construct 3d interpolating spline for phi
 !--------------------------------------------
 call Spline_Interp_YB(phi)
 
 !----------------------------------------------------------
 ! Evaluate interpolant for all the "special" image points
 ! and communicate values globally ...
 !----------------------------------------------------------
 if( numprocs > 1 ) then
  xtmp(:) = 0.d0
  fval(:) = 0.d0
  N = total_num_remote_image_points
  !------------------------------------------
  ! For each special point, the image owner
  ! interpolates phi to the image point, storing
  ! the result in the array xtmp
  !------------------------------------------
  do m=1,N
   if( owner_image(m)==myid ) then     
    xi = xyz_image(m,1)
    yi = xyz_image(m,2)
    zi = xyz_image(m,3)
    call Spline_Eval_YB_ptwise( xi,yi,zi, xtmp(m) )  ! tmpY(:,:,:,1) = workspace 
   endif
  enddo  
  !-----------------------------------------------------
  ! Make all the results globally available in fval(:)
  !------------------------------------------------------
  call mpi_allreduce(xtmp,fval,N,mpi_double_precision,mpi_sum,comm,ierr)  ! this is blocking
 else
  N=0
 endif

!-------------------------------------------------------------------------------------
! Note: this logic depends on passing through the YBLOCK data in the same order,
! (i.e. k,i,j) as was done in SetupImmersedBoundary, with identical logic regarding
! TOL so that the special points are encountered in exactly the same order, otherwise
! the index in the lookup table will be incorrect.
!-------------------------------------------------------------------------------------
 m = my_start_index     ! starting location in N arrays for processor myid
 do k=1,array_size(KDIM,YBLOCK,myid)
  do i=1,array_size(JDIM,YBLOCK,myid)
   do j=1,array_size(IDIM,YBLOCK,myid)
    kg = global_z_indices(START,YBLOCK,myid) + k - 1
    
    if( dist(j,i,k) < 0  ) then   ! pt is outside fluid domain 
    
     d = abs( dist(j,i,k) )  ! work w/ pos value...
   
     if( d < TOL ) then  ! pt is close to IB    TOL = 25 sigma_n  
          
      !----------------------------------------------------
      !  find the image value inside the fluid domain
      !  2*d to get to the surface and equal dist beyond...
      !----------------------------------------------------
    
       xi = x(i)  + 2.d0*d*nhat(j,i,k,1)
       yi = y(j)  + 2.d0*d*nhat(j,i,k,2)
       zi = z(kg) + 2.d0*d*nhat(j,i,k,3)
       
    
      !---------------------------------------------------------------------
      ! evaluate interpolating functions to get phi at the image point
      !  if myid owns the image point, do the evaluation, otherwise, look up
      !  the function value from the precomputed values associated with
      !  all the "special" points
      !----------------------------------------------------------------------
      if( N > 0 ) then
       if( i==ijk_solid(m,1) .and. j==ijk_solid(m,2) .and. k==ijk_solid(m,3) ) then
        val = fval(m)   ! look up the interpolated value
        m = m + 1       ! increment position in table ...
       else
        call Spline_Eval_YB_ptwise( xi,yi,zi,val )  ! tmpY(:,:,:,1) = workspace
       endif
      else
       call Spline_Eval_YB_ptwise( xi,yi,zi,val )  ! tmpY(:,:,:,1) = workspace
      endif
    
      !----------------------------------------------
      ! make phi symmetric across immersed boundary
      !----------------------------------------------
    
       phi(j,i,k) = val 
    
      !----------------------------------------------------------------------
      ! force even extension to decay smoothly to 0 away from IB (into solid)
      !----------------------------------------------------------------------
     
       phi(j,i,k) = phi(j,i,k)*exp(-( d/sigma )**4)
       
     else  ! further from IB
      
       !---------------------------------------------------------------------
       ! Assume the smoothly decaying function is close enough to zero.
       ! This avoids interpolating when the answer won't have any effect.
       !---------------------------------------------------------------------
       phi(j,i,k) = 0.d0
       
     endif  ! end close/far from IB block
                                                
    endif  ! end outside fluid domain block
   
   enddo
  enddo
 enddo
 
return
end subroutine homogeneous_neumann_blend


subroutine hc_dirichlet_blend(phi)
 !--------------------------------------------------
 !--------------------------------------------------
 ! inhomogeneous dirichlet BC for s1 = rho at IB
 !--------------------------------------------------
 !--------------------------------------------------
 use immersed_boundary,      only: dist,nhat,TOL
 use independent_variables,  only: x,y,z,Ly
 use mpi_params,             only: myid,comm,ierr 
 use dimensional_scales,     only: density_scale,length_scale
 use decomposition_params  
 
 implicit none 
 include 'mpif.h'
 logical,save                   :: first_entry = .TRUE.
 logical,save                   :: IB_owner
 integer                        :: k0,k1
 integer                        :: i,j,k,kg
 integer,save                   :: Nvals
 real(kind=8),save              :: z0,z1,zs,pi
 real(kind=8),allocatable,save  :: fvals(:,:,:)
 real(kind=8)                   :: xs,ys,d
 real(kind=8)                   :: g1,g2,dz,sigma,f0
 real(kind=8)                   :: Rmax,ky,S
 real(kind=8)                   :: xi,yi,zi,val,phi_0
 real(kind=8)                   :: phi( array_size(IDIM,YBLOCK,myid),   &
                                        array_size(JDIM,YBLOCK,myid),   &
                                        array_size(KDIM,YBLOCK,myid)  )
 
 
 Rmax = 0.60874/100.d0              ! [kg/m3]  maximum density at IB surface
 Rmax = Rmax/density_scale          ! [1] dimensionless maximum density 
 ky = (  pi/(1.*Ly) )*length_scale  ! [1] dimensionless y wavenumber 
 dz = z(2)-z(1)                     ! [1] dless dz, we have a flat IB surface for horiz. conv problem 
 
 if( first_entry ) then
  pi=4.d0*datan(1.d0)
 
  !--------------------------------------------------------------
  ! slightly extended z domain for interpolation responsibilities
  !--------------------------------------------------------------  
  k0 = global_z_indices(START,YBLOCK,myid)
  k1 = global_z_indices(END,YBLOCK,myid)
  z0 = z(k0) - dz/2.d0   ! dless lower bound  of "interpolation responsibility" region
  z1 = z(k1) + dz/2.d0   ! dless upper bound  of "interpolation responsibility" region
  
  !--------------------------------
  ! locate zs for flat IB
  !--------------------------------
  d = dist(1,1,1)                            ! distance to IB from lowest z pt on each processor
  kg = global_z_indices(START,YBLOCK,myid)   ! lowest z pt on each processor
  if( d >= 0 ) then   ! in fluid
   zs = z(kg) + d
  else                ! in solid
   zs = z(kg) - d
  endif
  
  !--------------------------------
  ! find out if myid is "owner"
  !--------------------------------
  if( z0 <= zs .and. z1 > zs ) then
   IB_owner = .TRUE.
  else
   IB_owner = .FALSE.
  endif
  
  ! Allocate array(s) for interpolated function values/global reduction
  allocate( fvals(array_size(IDIM,YBLOCK,myid),array_size(JDIM,YBLOCK,myid),2) )    
  Nvals = array_size(IDIM,YBLOCK,myid)*array_size(JDIM,YBLOCK,myid)
  
  first_entry=.FALSE.
 endif
  
 
 !------------------------------------------
 ! construct 3d interpolating spline for s1
 !------------------------------------------ 
  call Spline_Interp_YB(phi)   
 
 !------------------------------------------
 ! get the function values at z=zs
 !------------------------------------------ 
 if( IB_owner ) then
  fvals(:,:,1)=0.d0
  do i=1,array_size(JDIM,YBLOCK,myid)
   do j=1,array_size(IDIM,YBLOCK,myid)
    call Spline_Eval_YB_ptwise( x(i),y(j),zs,fvals(j,i,2) )  ! tmpY(:,:,:,1) = workspace
   enddo
  enddo
 endif
 
 !----------------------------------------------------------
 ! Make the results globally available in fvals(:,:,1)
 !----------------------------------------------------------
  call mpi_allreduce(fvals(1,1,2),fvals(1,1,1),Nvals,mpi_double_precision,mpi_sum,comm,ierr) 
 
 
 do k=1,array_size(KDIM,YBLOCK,myid)
  do i=1,array_size(JDIM,YBLOCK,myid)
   do j=1,array_size(IDIM,YBLOCK,myid)
   
    d = dist(j,i,k)
    if( abs(d) < TOL  ) then   ! pt is near IB, either side
                 
       S  = Rmax * ( 1.d0 - cos(ky*y(j))**2  )    ! BC fn of y only
       f0 = fvals(j,i,1)
                  
       !----------------------------------------------------
       !  add localized adjustment to satisfy BC
       !----------------------------------------------------       
        sigma = dz
        if( d >= 0.d0 ) then   ! fluid
         phi(j,i,k) = phi(j,i,k) + (S-f0)*exp( -( d /sigma) )
        else                   ! solid
         phi(j,i,k) = phi(j,i,k) + (S-f0)*exp( -( d /sigma) - ( d /(4*sigma) )**4 )
        endif
                                  
    endif  ! end close/far from IB block
                                                
        
    
   enddo
  enddo
 enddo
 
  
return
end subroutine hc_dirichlet_blend




subroutine stokes2_dirichlet_blend
 !--------------------------------------------------
 !--------------------------------------------------
 ! inhomogeneous dirichlet BC for v at IB
 !   PRELIMINARY, THIS NEEDS TO BE REVISITED...
 !     NEEDS MULTIPROCESSOR LOGIC
 !--------------------------------------------------
 !--------------------------------------------------

 use immersed_boundary,      only: dist,nhat,sigma_d,TOL
 use dependent_variables,    only: u,v,w
 use independent_variables,  only: x,y,z,t_secs
 use mpi_params,             only: myid 
 use dimensional_scales,     only: velocity_scale,time_scale,length_scale
 use decomposition_params                                      
 implicit none 
 integer                        :: i,j,k,kg,N
 real(kind=8)                   :: xs,ys,zs,d,pi
 real(kind=8)                   :: g,sigma,f0,S0,S,urv,omega
 
 pi=4.d0*datan(1.d0)
   
 !------------------------------------------
 ! inhomogeneous dirichlet for v
 !  v = V0*cos(omega*t) at IB z=z0
    S0 = 0.10/velocity_scale               ! dless velocity scale 
    omega = (2.d0*pi/10.)                  ! frequency  [1/s] , use w/ time in seconds...
    
    S0 = S0*cos(omega*t_secs)              ! current time's BC for v, dless
    sigma= (z(2)-z(1))/2.d0                ! 1d problem, periodic Ly,dy meaningless...
 !------------------------------------------
 
 !------------------------------------------
 ! construct 3d interpolating spline for v
 !------------------------------------------
 call Spline_Interp_YB(v)
 
 do k=1,array_size(KDIM,YBLOCK,myid)
  kg = global_z_indices(START,YBLOCK,myid) + k - 1
  do i=1,array_size(JDIM,YBLOCK,myid)
   do j=1,array_size(IDIM,YBLOCK,myid)
   
    d = dist(j,i,k)
    if( abs(d) <= TOL ) then  ! near IB
    
    
      !----------------------------------------------------
      !  find the nearest IB surface point
      !----------------------------------------------------
    
       xs = x(i)  - d*nhat(j,i,k,1)
       ys = y(j)  - d*nhat(j,i,k,2)
       zs = z(kg) - d*nhat(j,i,k,3)
    
      !-----------------------------------------------------------------------
      ! evaluate interpolating functionto get v at the nearest surface point
      ! ALL VALUES ARE LOCAL to this processor in this example...
      !-----------------------------------------------------------------------
    
       call Spline_Eval_YB_ptwise( xs,ys,zs,f0 )  ! uses tmpY(:,:,:,1) for some workspace...
       
       call RANDOM_NUMBER( urv )        ! in [0,1]
       urv = 2.*(urv-.5)                ! in [-1 1]
       
       S = S0*(1. + urv/1000.)
       
       !----------------------------------------------------
       !  add "delta function" adjustment to satisfy BC
       !----------------------------------------------------
       v(j,i,k) = v(j,i,k) + (S-f0)*exp( -(d/sigma)**2 )
          
      
    endif
    
        
   enddo
  enddo
 enddo
 
  
return
end subroutine stokes2_dirichlet_blend



