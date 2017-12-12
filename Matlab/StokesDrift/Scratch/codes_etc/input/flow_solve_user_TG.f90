subroutine eval_s1_bar(z,s1_bar,s1_bar_z,s1_bar_zz)  
use independent_variables, only: Lz
use user_params,           only: delta_rho,delta,d,z0,zr !! see user_params_module.f90 for definitions
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
 real(kind=8),save        :: sech2
 


  if (first_entry) then
   density_gradient = 2.6e-3   
   z0 = 2*Lz/3.d0                   !! shear interface location [m]
   zr = z0 - d                    !! density interface position   
   first_entry = .FALSE.
  endif
    sech2 = 1.d0 / cosh(2.*(z-zr)/(delta*Lz))**2


    s1_bar = 0.5*density_gradient*Lz*(1.-tanh(2*(z-zr)/(delta*Lz)))  
    s1_bar_z = -density_gradient*Lz*(1./(Lz*delta))*sech2
    s1_bar_zz =4*density_gradient*Lz/(1./((delta*Lz)**2))*sech2*tanh(2.*(z-zr)/(Lz*delta))

!    if(z .lt. (Lz-zr)) then
!    s1_bar=0.5d0*density_gradient*Lz*(1.-tanh(2*(z-zr)/(delta*Lz))) +density_gradient*((Lz-zr)-z)
!    s1_bar_z = -density_gradient*Lz*(1./(Lz*delta))*sech2 -density_gradient
!    s1_bar_zz = 4*density_gradient*Lz/(1./((delta*Lz)**2))*sech2*tanh(2.*(z-zr)/(Lz*delta))
!    endif

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
 logical,save             :: first_entry=.TRUE.
 
 !------------------------------------------------------
 ! additional user declared parameters
 !------------------------------------------------------
 real(kind=8),save        :: sech2
 
 if( first_entry ) then 
   z0 = Lz/2.d0                   !! shear interface location [m]
   zr = z0 - d                    !! density interface position 
   first_entry=.FALSE.
 endif
 
 sech2 = 1.d0 / cosh(2.*(z-zr)/delta)**2
 
 s2_bar = (1/2.)*(1.d0 - tanh(2.*(z-zr)/delta))
 s2_bar_z = -(1./delta)*sech2
 s2_bar_zz = 4*(1./(delta**2))*sech2*tanh(2.*(z-zr)/delta)
      
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
 use mpi_params,             only: myid,comm,ierr,numprocs
 use independent_variables,  only: Lx,Ly,Lz
 use dimensional_scales,     only: coriolis=>f,dgrad,rho_0,g,N2
 use decomposition_params,   only: YBLOCK,proc_row,proc_col
 use user_params,            only: delta_rho,delta_U,h,delta,d,z0,zr,amp_pert,    &
                                   up_restart,up_rs_basename, up_rs_islice,       &
                                   up_subtract_s1_bar,up_subtract_s2_bar  
!! see user_params_module.f90 for definitions

 
 implicit none
 integer                  ::  id
 integer                  ::  i,j,k
 integer                  ::  nx,ny,nz
 real(kind=8)             ::  pi
 real(kind=8)             ::  x(nx)
 real(kind=8)             ::  y(ny)
 real(kind=8)             ::  z(nz)
 real(kind=8)             ::  F(ny,nx,nz)
 real(kind=8),allocatable ::  Fscratch(:,:,:)
! real(kind=8)             :: Fscratch(nx,ny,nz)
 logical                  ::  restart 
 
 !------------------------------------------------------------
 ! additional user's parameter declarations 
 !------------------------------------------------------------
 integer                  :: ncid,varid
 character(len=7)         :: cid
 character(len=3)         :: c_row_id
 character(len=3)         :: c_col_id
 character(len=6)         :: cslice
 character(len=80)        :: ncfile
 integer                  :: start(4),count(4)  !! time,kdimension,jdimension,idimension stored 
 real(kind=8)             :: fzramp(nz)
 real(kind=8)             :: urv,A,kx,ky,alphaz,amp_noise
 real(kind=8)             :: s_bar,s_bar_z,s_bar_zz
 include 'netcdf.inc'
 restart = up_restart            ! [1]
 !------------------------------------------------------------
 if (.not. restart ) then
  pi=4.d0*atan(1.d0)
  kx=8.d0*(2.d0*pi/Lx)
!  kx=4.d0*(2.d0*pi/Lx)
  ky=kx
  alphaz=1.d0*(2.d0*pi/Lz)
  A = 0.2
  amp_noise=3.e-1
  z0 = 2*Lz/3.d0                   !! shear interface location [m]
   zr = z0                     !! density interface position


  do k=1,nz
!    fzramp(k) = exp(-alphaz*(Lz-z(k)))
    fzramp(k) = 1.d0
  enddo

  
 select case(id) 
  case(1)        !! ics for u [m/s]
  do i=1,nx
    do j=1,ny
     do k=1,nz
 
      call RANDOM_NUMBER(urv) ! random number between [0,1]
      urv= 2.*(urv-0.5)        ! rv between [-0.5,0.5]
      F(j,i,k) = A*fzramp(k)*(cos(kx*x(i))*sin(ky*y(j))+amp_noise*urv*sin(4*kx*x(i))*cos(2*ky*y(j)))
     enddo
    enddo
  enddo
 
  case(2)        !! ics for v [m/s]
  
   do i=1,nx
    do j=1,ny
     do k=1,nz
      call RANDOM_NUMBER(urv) ! random number between [0,1]
      urv= 2.*(urv-0.5)        ! rv between [-0.5,0.5]
      F(j,i,k) = -A*fzramp(k)*(kx/ky)*(sin(kx*x(i))*cos(ky*y(j))+2.d0*amp_noise*urv*cos(4*kx*x(i))*sin(2*ky*y(j)))
      enddo
    enddo
   enddo
  
  case(3)        !! ics for w [m/s]
  
  F(:,:,:)=0.d0 
  case(4)        !! ics for s1' 
  
   !--------------------------------------------------------------------------------
   !  add a perturbation to rho that is consistent with the boundary conditions...
   !--------------------------------------------------------------------------------
   do i=1,nx
    do j=1,ny
     do k=1,nz
      call RANDOM_NUMBER(urv) ! random number between [0,1]
      urv= 2.*(urv-0.5)        ! rv between [-0.5,0.5]
!      F(j,i,k) = -A*coriolis(1)*(rho_0/g)*(alphaz/kx)*fzramp(k)*(cos(kx*x(i))*cos(ky*y(j))+0.5d0*amp_noise*urv*sin(4*kx*x(i))*sin(2*ky*y(j)))
       F(j,i,k) = 0.d0
 
     enddo
    enddo
   enddo
  case(5)        !! ics for s2' 
  
   F(:,:,:)=0.d0
   case default
    stop 'error in calling user_ics: id out of range'
  end select

elseif (restart) then
allocate(Fscratch(nx,ny,nz))
!!=============================================================
!! construct the filename
!!=============================================================
  write(unit=c_row_id, fmt = 101) proc_row(YBLOCK,myid)
  write(unit=c_col_id, fmt = 101) proc_col(YBLOCK,myid)
  cid = c_row_id//'-'//c_col_id
101 format(I3.3)
  write(unit=cslice,fmt=102) up_rs_islice
102 format(I6.6)
    if (numprocs > 1) then
     ncfile=trim(up_rs_basename)//'_'//cid//'.nc'
    elseif (numprocs==1) then
     ncfile=trim(up_rs_basename)//'.nc'
    endif
    if (myid==0 .and. id==1) then
        write(0,*) '=============================================='
        write(0,*) 'READING ICS from restart files'
        write(0,*) '                         time slic...',up_rs_islice
        write(0,*) '=============================================='
    endif
!------------------------------------------------------------------
! ncdump -h reports e.g. u(timedimension,kdimension,jdimension,idimension)
! but the order of the indices needs to be reversed here for fortran
!-----------------------------------------------------------------
start=(/1,1,1,1/)    ! "rs_islice"-th slice saved time record
count = (/nx,ny,nz,1/)       ! 1 time slice, nx (global), ny (global),nz=nz(myid)==> (3D YBLOCK)
write(0,*) 'nx,ny,nz,ncfile',nx,ny,nz,ncfile
!!==========================================================
!! open the file, check for error
!!========================================================== 
ierr=NF_OPEN(ncfile,NF_NOWRITE,ncid)
if (ierr .ne. NF_NOERR) then
 write(0,*) '... ERROR OPENING NETCDF FILE: user_ics',trim(ncfile)
 write(0,*) '... mydi,ierr ',myid,ierr
 stop
endif

select case(id)
case(1)
!!==========================================================
!! extract variable id, check for error
!!==========================================================
ierr=NF_INQ_VARID(ncid,'u',varid)
if (ierr .ne. NF_NOERR) then
   write(0,*) myid,'NetCDF ERROR INQ_VARID -> u: ',ierr
stop
endif

!!==========================================================
!! read the corresponding variable (dimensional)
!!==========================================================
ierr=NF_GET_VARA_DOUBLE(ncid,varid,start,count,Fscratch)
if (ierr .ne. NF_NOERR) then
   write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> u: ',ierr
stop
endif
write(0,*)'have read u'
do k=1,nz
 do j=1,ny
  do i=1,nx
  F(j,i,k) = Fscratch(i,j,k)
  enddo
 enddo 
enddo
case(2)
!!==========================================================
!! extract variable id, check for error
!!==========================================================
ierr=NF_INQ_VARID(ncid,'v',varid)
if (ierr .ne. NF_NOERR) then
   write(0,*) myid,'NetCDF ERROR INQ_VARID -> v: ',ierr
stop
endif

!!==========================================================
!! read the corresponding variable (dimensional)
!!==========================================================
ierr=NF_GET_VARA_DOUBLE(ncid,varid,start,count,Fscratch)
if (ierr .ne. NF_NOERR) then
   write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> v: ',ierr
stop
endif
do k=1,nz
 do j=1,ny
  do i=1,nx
  F(j,i,k) = Fscratch(i,j,k)
  enddo
 enddo
enddo


case(3)
!!==========================================================
!! extract variable id, check for error
!!==========================================================
ierr=NF_INQ_VARID(ncid,'w',varid)
if (ierr .ne. NF_NOERR) then
   write(0,*) myid,'NetCDF ERROR INQ_VARID -> w: ',ierr
stop
endif

!!==========================================================
!! read the corresponding variable (dimensional)
!!==========================================================
ierr=NF_GET_VARA_DOUBLE(ncid,varid,start,count,Fscratch)
if (ierr .ne. NF_NOERR) then
   write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> w: ',ierr
stop
endif

do k=1,nz
 do j=1,ny
  do i=1,nx
  F(j,i,k) = Fscratch(i,j,k)
  enddo
 enddo
enddo

case(4)
!!==========================================================
!! extract variable id, check for error
!!==========================================================
ierr=NF_INQ_VARID(ncid,'s1',varid)
if (ierr .ne. NF_NOERR) then
   write(0,*) myid,'NetCDF ERROR INQ_VARID -> s1: ',ierr
stop
endif

!!==========================================================
!! read the corresponding variable (dimensional)
!!==========================================================
ierr=NF_GET_VARA_DOUBLE(ncid,varid,start,count,Fscratch)
if (ierr .ne. NF_NOERR) then
   write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> s1: ',ierr
stop
endif
do k=1,nz
 do j=1,ny
  do i=1,nx
  F(j,i,k) = Fscratch(i,j,k)
  enddo
 enddo
enddo
!!=========================================================
!! s1_bar+s1' in netcdf file ==> subtract s1_bar
!! whether this is the case depends on what was stored in
!! the original file
!!========================================================
if( up_subtract_s1_bar ) then
do k=1,nz
   call eval_s1_bar( z(k), s_bar, s_bar_z, s_bar_zz )
   F(:,:,k) = F(:,:,k) - s_bar
enddo
endif

case(5)
!!==========================================================
!! base run had no second scalar s2. s2 not computed or stored 
!! therefore set it equal to 0.
!!==========================================================
F(:,:,:) = 0.d0 

end select
deallocate(Fscratch)
endif  ! end of if restart logic...

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
 use dimensional_scales,     only: velocity_scale,scalar_scale,length_scale 
 use independent_variables,  only: Lx,Ly,Lz,global_nz=>nz
 use dependent_variables,    only: u,v,w,s1,s2
 use mpi_params,             only: myid,comm,ierr
  
 implicit none
 real(kind=8),dimension(:,:,:),pointer  :: phi
 integer                                :: i,j,k
 integer                                :: id,nx,ny,nz
 real(kind=8)                           :: x(nx),y(ny),z(nz),t
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
  ! user defined forcing parameters
    real(kind=8)                     :: wind_amp
    real(kind=8)                     :: ml_depth
    real(kind=8)                     :: alpha
  !-----------------------------------------------------------------------------------
 
 
 !------------------------------------------------------------------------------------
 !------------------------------------------------------------------------------------
 if(first_entry) then
  pi=4.*atan(1.)
  call_again = .FALSE.  ! switch to true below if desired, decide for each id
  
  first_entry=.TRUE.
 endif
 !------------------------------------------------------------------------------------
 !------------------------------------------------------------------------------------
 
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

   
 if( id == 1) then       !  FF <===  specify desired forcing on rhs of u equation
 wind_amp = 1.e-5
 ml_depth = 400.d0
 alpha = 1.d0/(Lz-ml_depth)
   if(t .gt. 2.e6) then
     do i=1,nx
        do j=1,ny
         do k=1,nz
         FF(j,i,k) = wind_amp*exp(-10*alpha*(Lz-z(k)))
         enddo
        enddo
       enddo
    else
         FF(:,:,:) = 0.d0
    endif

  call_again = .TRUE.   !   keep calling routine for u

  
 elseif( id==2 ) then    !  FF <===  specify desired forcing on rhs of v equation
  FF(:,:,:) = 0.d0  
  call_again = .FALSE.

  
 elseif( id==3 ) then    !  FF <===  specify desired forcing on rhs of w equation
  FF(:,:,:) = 0.d0   
  call_again = .FALSE.


 elseif( id == 4 ) then  !  FF <===  specify desired forcing on rhs of s1 equation
  FF(:,:,:) = 0.d0
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

