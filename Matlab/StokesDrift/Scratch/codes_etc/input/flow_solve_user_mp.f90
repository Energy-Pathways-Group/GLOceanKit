subroutine eval_s1_bar(z,s1_bar,s1_bar_z,s1_bar_zz)  
use independent_variables, only: Lz
use dimensional_scales, only: g,rho_0
!use user_params,           only: N0,zm,B !! see user_params_module.f90 for definitions
!---------------------------------------------------------------- 
!  Code is written to solve the eqn of motion for s1'
!   where s1 = s1_bar(z) + s1'(x,y,z,t).
!   Supply s1_bar and its z derivatives here.
!   Setting s1_bar=0 is fine if there is no sensible 
!   ambient distribution of s1 in your problem.
!
! z:           vertical z spatial position  [m]
! output 
! s1_bar, s1_bar_z, s1_bar_zz
! scalar concentration and z derivs   [dimensional units]
!---------------------------------------------------------------- 
 implicit none 
 real(kind=8)             :: s1_bar,s1_bar_z,s1_bar_zz,z
 logical,save             :: first_entry=.TRUE.
 
 !------------------------------------------------------
 ! additional user declared parameters
   
 !------------------------------------------------------
 real(kind=8),save        :: density_gradient
 real(kind=8)             ::dRdz


  if (first_entry) then
     density_gradient = 2.9e-3   
     first_entry = .FALSE.  
  endif
!     zm  = Lz
!     dRdz = N0*N0*rho_0/g
!     s1_bar = rho_0 + (B/2.)*dRdz*(1.d0 -exp(-2*(zm-z)/B))
!     s1_bar_z = dRdz*exp(-2*(zm-z)/B)
!     s1_bar_zz = dRdz*(2/B)*exp(-2*(zm-z)/B)

    s1_bar = density_gradient*(Lz-z)  
    s1_bar_z = -density_gradient
    s1_bar_zz =0.d0


 return
end subroutine eval_s1_bar

 
 




subroutine eval_s2_bar(z,s2_bar,s2_bar_z,s2_bar_zz) 
use independent_variables,   only: Lz
use user_params,             only: delta,d,z0,zr  !! see user_params_module.f90 for definitions
!---------------------------------------------------------------- 
!  Code is written to solve the eqn of motion for s2'
!   where s2 = s2_bar(z) + s2'(x,y,z,t).
!   Supply s2_bar and its z derivatives here.
!   Setting s2_bar=0 is fine if there is no sensible 
!   ambient distribution of s2 in your problem.
!
! z:           vertical z spatial position  [m]
! output 
! s2_bar, s2_bar_z, s2_bar_zz
! scalar concentration and z derivs   [dimensional units]
!---------------------------------------------------------------- 
implicit none
 real(kind=8)             :: s2_bar,s2_bar_z,s2_bar_zz,z
 logical,save             :: first_entry=.FALSE.
 
 !------------------------------------------------------
 ! additional user declared parameters
 !------------------------------------------------------
 
 if( first_entry ) then 
   z0 = Lz/2.d0                   !! shear interface location [m]
   zr = z0 - d                    !! density interface position 
   first_entry=.FALSE.
 endif
 
 
 s2_bar = 0.d0
 s2_bar_z = 0.d0
 s2_bar_zz = 0.d0
      
 return
end subroutine eval_s2_bar

subroutine user_ics(x,y,z,F,id,nx,ny,nz)
!!-----------------------------------------------------------
!!  Inputs and outputs are all in dimensional units.
!!
!!  inputs:
!!   nx,ny,nz  array size parameters
!!   x,y,z     coordinate arrays (in meters) that define 
!!             the region over which forcing functions
!!             are to be evaluated
!! 
!!   id        field id [u,v,w,s1',s2'] <---[1,2,3,4,5]
!!             where s1 = s1bar(z) + s1'(y,x,z,t)
!!             and   s2 = s2bar(z) + s2'(y,x,z,t)
!!
!!  outputs:
!!          F(j,i,k)  i.e. F(y,x,z)
!!          the initial conditions for variable id
!!          F is to be returned in DIMENSIONAL units 
!!          (i.e. m/s or deg C etc )
!!
!!  NB STORAGE INDEX CONVENTION: F(y,x,z)
!!   
!!-----------------------------------------------------------
 use mpi_params,             only: myid,comm
 use independent_variables,  only: Lx,Ly,Lz
 use dimensional_scales,     only: coriolis=>f,dgrad,rho_0,g,N2
 use decomposition_params,   only: YBLOCK,proc_row,proc_col
 use user_params,            only: delta_rho,delta_U,h,delta,d,z0,zr,amp_pert  !! see user_params_module.f90 for definitions

 
 implicit none
 integer                  ::  id
 integer                  ::  i,j,k,ii
 integer                  ::  nx,ny,nz
 real(kind=8)             ::  pi
 real(kind=8)             ::  x(nx)
 real(kind=8)             ::  y(ny)
 real(kind=8)             ::  z(nz)
 real(kind=8)             ::  F(ny,nx,nz)

 
 !------------------------------------------------------------
 ! additional user's parameter declarations 
 !------------------------------------------------------------
 real(kind=8),allocatable     :: temp_u(:,:)
 real(kind=8),allocatable     :: temp_v(:,:)
 real(kind=8),allocatable     :: temp_w(:,:)
 real(kind=8),allocatable     :: temp_pd(:,:)
 real(kind=8)     :: fzramp(nz)
 real(kind=8)     :: urv,A,kx,ky,alphaz,amp_noise,strain_rate,k0
 real(kind=8)     :: dye_alpha,dye_beta,x00,z00,dye_radial_ext,dye_vert_ext
 integer          :: nbytes,len
 character(len=35) :: filename
 !------------------------------------------------------------
  pi=4.d0*atan(1.d0)
  kx=8.d0*(2.d0*pi/Lx)
!  kx=4.d0*(2.d0*pi/Lx)
  strain_rate  = 1.e-5
  ky=kx
  alphaz=1.d0*(2.d0*pi/Lz)
  A = 0.0
  amp_noise=0.e-1
  z0 = 2*Lz/3.d0                   !! shear interface location [m]
   zr = z0                     !! density interface position

allocate(temp_u(nx*nz,ny))
allocate(temp_v(nx*nz,ny))
! reading in ascii data from matlab
! construct filename  
write(0,*) 'myid,nz',myid,nz
write(filename,"('cyprus_u_',I2.2,'.txt')")myid
write(0,*)filename
 open(10,file=filename,status='old',form='formatted')
if (myid==15) then
 read(10,*) ((temp_u(i,j),i=1,nx*(nz-1)),j=1,ny)
 close(unit=10)
write(filename,"('cyprus_v_',I2.2,'.txt')")myid
 open(11,file=filename,status='old',form='formatted')
 read(11,*) ((temp_v(i,j),i=1,nx*(nz-1)),j=1,ny)
 close(unit=11)
else
 read(10,*) ((temp_u(i,j),i=1,nx*nz),j=1,ny)
 close(unit=10)
write(filename,"('cyprus_v_',I2.2,'.txt')")myid
 open(11,file=filename,status='old',form='formatted')
 read(11,*) ((temp_v(i,j),i=1,nx*nz),j=1,ny)
 close(unit=11)
endif
 select case(id) 
  case(1)        !! ics for u [m/s]
 do k=1,nz-1
    do j=1,ny
       do i=1,nx
      F(j,i,k) = temp_u(i*k,j)
     enddo
    enddo
  enddo
  case(2)        !! ics for v [m/s]
 do k=1,nz-1
    do j=1,ny
     do i=1,nx
      F(j,i,k) = temp_v(i*k,j)
     enddo
    enddo
  enddo
  case(3)        !! ics for w [m/s]
  
  do i=1,nx
    do j=1,ny
     do k=1,nz
      F(j,i,k)=0.d0
      enddo
    enddo
   enddo

  case(4)        !! ics for s1' 
  
   !--------------------------------------------------------------------------------
   !  add a perturbation to rho that is consistent with the boundary conditions...
   !--------------------------------------------------------------------------------
   do i=1,nx
    do j=1,ny
     do k=1,nz
        F(j,i,k)=0.d0 
     enddo
    enddo
   enddo
  case(5)        !! ics for s2' 
  
   x00 = Lx/2.
   z00 = Lz/2. 
   dye_radial_ext = 0.025*Lx
   dye_vert_ext = 0.05*Lz
   dye_alpha = 1.d0/(2*dye_radial_ext**2)
   dye_beta = 1.d0/(2*dye_vert_ext**2)
   k0 = 2*pi/Lx
   do i=1,nx
    do j=1,ny
     do k=1,nz
!        F(j,i,k)=1.d0*exp(-dye_alpha*(x(i)-x00)**2 -dye_beta*(z(k)-z00)**2)
!     F(j,i,k) = 1.d0*cos(k0*x(i))
     enddo
    enddo
   enddo
   case default
    stop 'error in calling user_ics: id out of range'
  end select

deallocate(temp_u)
deallocate(temp_v)
return
end subroutine user_ics







subroutine user_forcing(x,y,z,t,FF,     &
                        id,nx,ny,nz,    &
                        call_again)
!!  User defined subroutine to add time/space dependent
!!  forcing to the equations of motion. Inputs and outputs
!!  for this routine are all in dimensional units.
!!  
!!
!!  inputs:
!!   nx,ny,nz  array size parameters
!!   x,y,z     coordinate arrays [m] that define 
!!             the grid points at which the forcing
!!             functions are to be evaluated
!!   t         current time in seconds  
!!   id        field id [u,v,w,s1',s2'] <---[1,2,3,4,5]
!!
!!  outputs:
!!            FF(j,i,k) i.e. FF(y,x,z)
!!            the rhs term for equation id
!!            FF is to be returned in DIMENSIONAL units 
!!            (i.e. m/s2 or deg C/s etc )
!!   
!!            call_next_time
!!            set false if field id is not being forced
!!            e.g. suppose do_forcing=.true., but only u 
!!            is actually forced, set call_next_time=.FALSE.
!!            for v,w,s1 & s2 to prevent adding zeros to
!!            the respective rhs's at each time step.
!!
!!  NB STORAGE INDEX CONVENTION: FF(y,x,z)
!! 
!!===========================================================
 use dimensional_scales,     only: velocity_scale,scalar_scale,length_scale,coriolis=>f,dgrad,rho_0,g,N2 
 use independent_variables,  only: Lx,Ly,Lz,global_nz=>nz
 use dependent_variables,    only: u,v,w,s1,s2
 use mpi_params,             only: myid,comm,ierr
! use user_params,            only: nwaves 
 implicit none
 real(kind=8),dimension(:,:,:),pointer  :: phi
 integer                                :: i,j,k
 integer                                :: id,nx,ny,nz
 real(kind=8)                           :: x(nx),y(ny),z(nz),t
 real(kind=8)                           :: s1_bar_z(nz)
 real(kind=8)                           :: FF(ny,nx,nz) 
 real(kind=8)                           :: pi
 logical                                :: call_again
 real(kind=8),save                      :: scale
 logical,save                           :: first_entry=.TRUE.
  
 !------------------------------------------------------------------------------------
 ! user variables
 !------------------------------------------------------------------------------------
 !------------------------------------------------------------------------------------
 !-----------------------------------------------------------------------------------
  ! user defined forcing parameters and variables
    real(kind=8)                     :: kz,alphat
    real(kind=8)                     :: kx,u0,c_hat
    real(kind=8)                     :: omega,d_omega
    real(kind=8)                     :: ftramp,urv
    integer                          ::nn
  !-----------------------------------------------------------------------------------
  if(first_entry) then
    pi=4.*atan(1.)
    call_again = .FALSE.  ! switch to true below if desired, decide for each id

    first_entry=.TRUE.
  endif

    kz = 6.3e-3     ! corresponds to pi/Lz
    d_omega = coriolis(1) 
    alphat=coriolis(1)/6.28
 !------------------------------------------------------------------------------------
 !------------------------------------------------------------------------------------
 !------------------------------------------------------------------------------------
 !------------------------------------------------------------------------------------
 ftramp = tanh(alphat*t) 
 if( id==1 ) then
  phi=>u
  scale=velocity_scale
 elseif( id==2 ) then
  phi=>v
  scale=velocity_scale
 elseif( id==3 ) then
  phi=>w
  scale=velocity_scale
 elseif( id==4 ) then
  phi=>s1
  scale=scalar_scale(1)
 elseif( id==5 ) then
  phi=>s2
  scale=scalar_scale(2)
 endif

 FF(:,:,:) = 0.d0  
 if( id == 1) then       !  FF <===  specify desired forcing on rhs of u equation
  call_again = .FALSE.   !   keep calling routine for u

  
 elseif( id==2 ) then    !  FF <===  specify desired forcing on rhs of v equation
!  do i=1,nx
!      do j=1,ny
!          do k=1,nz
!          FF(j,i,k) = -(coriolis(1)/omega1)*u0*ftramp*sin(kx*x(i)+kz*z(k)-omega1*t)+(coriolis(1)/omega1)*u0*ftramp*sin(kx*x(i)+kz*z(k)+omega1*t)
!          
!          enddo
!      enddo
! enddo
  FF(:,:,:) = 0.d0
  call_again = .FALSE.

  
 elseif( id==3 ) then    !  FF <===  specify desired forcing on rhs of w equation
!  do i=1,nx
!     do j=1,ny
!        do k=1,nz
!        FF(j,i,k) = -(kx/kz)*u0*ftramp*cos(kx*x(i)+kz*z(k)-omega1*t)-(kx/kz)*u0*ftramp*cos(kx*x(i)+kz*z(k)+omega1*t)
!        enddo
!     enddo
!  enddo
  FF(:,:,:)=0.d0 
  call_again = .FALSE.


 elseif( id == 4 ) then  !  FF <===  specify desired forcing on rhs of s1 equation
!  do i=1,nx
!     do j=1,ny
!        do k=1,nz
!        FF(j,i,k) = -kx/(kz*omega1)*s1_bar_z(k)*u0*ftramp*sin(kx*x(i)+kz*z(k)-omega1*t)+kx/(kz*omega1)*s1_bar_z(k)*u0*ftramp*sin(kx*x(i)+kz*z(k)-omega1*t)
!        enddo
!     enddo
!  enddo
   FF(:,:,:)=0.d0
  call_again = .FALSE.


 elseif( id==5 ) then    !  FF <===  specify desired forcing on rhs of s2 equation
  FF(:,:,:) = 0.d0  
  call_again = .FALSE.
 endif 

end subroutine user_forcing


subroutine particle_positions(positions,npts,Lx,Ly,Lz)  
 use user_params,   only: h,delta,d,z0,zr     !! see user_params_module.f90 for definitions
 implicit none
 integer               :: i
 integer               :: npts                !! total number of particles
 real(kind=8)          :: positions(npts,3)   !! (1,2,3)->(x,y,z), in [m]
 real(kind=8)          :: Lx,Ly,Lz            !! domain size in [m]
 real(kind=8)          :: urv                 !! uniform rv in [0,1]
 
 !------------------------------------------------------------
 ! additional user's parameter declarations 
 !------------------------------------------------------------
 
 z0 = Lz/2.d0   ! shear center  [m]
 zr = z0 - d    ! strat center  [m]

!----------------------------------------------
! random_seed already called in initialize 
!----------------------------------------------
  do i=1,npts
         
   call RANDOM_NUMBER( urv )      ! in [0,1]
   positions(i,1) = Lx*urv        ! randomly in [0,lx]  
   
   call RANDOM_NUMBER( urv )      ! in [0,1]
   positions(i,2) = Ly*urv        ! randomly in [0,Ly]
   
   call RANDOM_NUMBER( urv )      ! in [0,1]
!   urv = 2.*(urv-0.5)             ! now in [-1 1]
!   positions(i,3) = zr  + 2*delta*urv  ! randomly near centerline
    positions(i,3) = Lz*urv        ! same as before
   
  enddo 

 return
end subroutine particle_positions


subroutine surface_height_derivs(x,y,h,hx,hy,hxx,hyy,hxy)
!----------------------------------------------------------------------
!  inputs are DIMENSIONAL x,y locations
!  output is the corresponding surface height h(x,y) in [m]
!  and various derivatives, also dimensional as necessary
!----------------------------------------------------------------------
 use independent_variables, only: Lx,Ly,Lz
 implicit none
 real(kind=8)      :: x,y,h,hx,hy,hxx,hyy,hxy
 real(kind=8)      :: pi
 
!----------------------------------------------------------------------
!    user defined immersed boundary parameters 
!----------------------------------------------------------------------
  pi = 4.d0*datan(1.d0)
!----------------------------------------------------------------------
   
  
!------------------------------------
! evaluate h(x,y) and derivs
!------------------------------------
 h = 0.d0                      ! [m]
 hx = 0.d0                     ! dimensionless
 hxx = 0.d0                    ! [1/m]
 hxy = 0.d0                    ! [1/m] 
 hy = 0.d0                     ! dimensionless
 hyy = 0.d0                    ! [1/m]
             
 return
end subroutine surface_height_derivs

