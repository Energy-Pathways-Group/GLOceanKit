subroutine eval_s1_bar(z,s1_bar,s1_bar_z,s1_bar_zz)  
use mpi_params,             only: myid,comm,ierr,numprocs
use dimensional_scales, 	only: rho_0
use decomposition_params,   only: global_z_indices
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
 integer(kind=8) nz,nbytes,len
 integer(kind=8),save    :: ncount
 
 !------------------------------------------------------------
 ! additional declarations for reading from netcdf files
 !------------------------------------------------------------
 include 'netcdf.inc'
 character(len=80)	      :: ncfile
 integer                  :: ncid,varid
 
 !------------------------------------------------------------
 ! provide netCDF filename with initial conditions.
 !------------------------------------------------------------
 ncfile = 'SaveIC_EarlyIWmodel_exp_strat_nonlin_correction.nc'
 
 !------------------------------------------------------------
 ! open the file, check for error
 !------------------------------------------------------------
 ierr=NF_OPEN(ncfile,NF_NOWRITE,ncid)
 if (ierr .ne. NF_NOERR) then
 	write(0,*) '... ERROR OPENING NETCDF FILE: eval_s1_bar ',trim(ncfile)
 	stop
 endif
 
 !------------------------------------------------------------
 ! keep track of vertical level for reading file
 !------------------------------------------------------------
 if (first_entry) then
 	ncount = 1  
 	first_entry = .FALSE.  
 endif
 
 !------------------------------------------------
 ! extract variable id, check for error
 !------------------------------------------------------------
 ierr=NF_INQ_VARID(ncid,'rhobar',varid)
 if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> s1_bar: ',ierr
 	stop
 endif
 !------------------------------------------------------------
 ! read the corresponding variable (dimensional)
 !------------------------------------------------------------
 ierr=NF_GET_VARA_DOUBLE(ncid,varid,ncount,1,s1_bar)
 if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> s1_bar: ',ierr	
 	stop
 endif

 !------------------------------------------------
 ! extract variable id, check for error
 !------------------------------------------------------------
 ierr=NF_INQ_VARID(ncid,'drhobar_dz',varid)
 if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> s1_bar_z: ',ierr
 	stop
 endif
 !------------------------------------------------------------
 ! read the corresponding variable (dimensional)
 !------------------------------------------------------------
 ierr=NF_GET_VARA_DOUBLE(ncid,varid,ncount,1,s1_bar_z)
 if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> s1_bar_z: ',ierr
 	stop
 endif
 
 !------------------------------------------------
 ! extract variable id, check for error
 !------------------------------------------------------------
 ierr=NF_INQ_VARID(ncid,'drhobar_dzdz',varid)
 if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> s1_bar_zz: ',ierr
 	stop
 endif
 !------------------------------------------------------------
 ! read the corresponding variable (dimensional)
 !------------------------------------------------------------
 ierr=NF_GET_VARA_DOUBLE(ncid,varid,ncount,1,s1_bar_zz)
 if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> s1_bar_zz: ',ierr
 	stop
 endif
 
 ncount = ncount+1
 ierr=NF_CLOSE(ncid)
 return
end subroutine eval_s1_bar

 
subroutine eval_s2_bar(z,s2_bar,s2_bar_z,s2_bar_zz) 
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
  
 s2_bar = 0.d0
 s2_bar_z = 0.d0
 s2_bar_zz = 0.d0
      
 return
end subroutine eval_s2_bar


subroutine user_ics(x,y,z,F,id,nx,ny,nz)
!!-----------------------------------------------------------
!!  Inputs and outputs are all i dimensional units.
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
 use decomposition_params,   only: global_z_indices

 
 implicit none
 integer                  ::  id
 integer                  ::  i,j,k,ii
 integer                  ::  nx,ny,nz
 real(kind=8)             ::  x(nx)
 real(kind=8)             ::  y(ny)
 real(kind=8)             ::  z(nz)
 real(kind=8)             ::  F(ny,nx,nz)
 
 !------------------------------------------------------------
 ! additional user's parameter declarations 
 !------------------------------------------------------------
 
 !------------------------------------------------------------
 ! additional declarations for reading from netcdf files
 !------------------------------------------------------------
 include 'netcdf.inc'
 real(kind=8),allocatable :: Fscratch(:,:,:)
 integer                  :: start_local(4),count(4)  !! time,kdimension,jdimension,idimension stored 
 character(len=80)	      :: ncfile
 integer                  :: ncid,varid

 !------------------------------------------------------------
 ! provide netCDF filename with initial conditions.
 !------------------------------------------------------------
 ncfile = 'SaveIC_EarlyIWmodel_exp_strat_nonlin_correction.nc'

 !------------------------------------------------------------
 ! open the file, check for error
 !------------------------------------------------------------
 ierr=NF_OPEN(ncfile,NF_NOWRITE,ncid)
 if (ierr .ne. NF_NOERR) then
	write(0,*) '... ERROR OPENING NETCDF FILE: user_ics ',trim(ncfile)
	write(0,*) '... myid,ierr ',myid,ierr
	stop
 endif

 !------------------------------------------------------------------
 ! ncdump -h reports e.g. u(timedimension,kdimension,jdimension,idimension)
 ! but the order of the indices needs to be reversed here for fortran.  Do 
 ! this for each variable (u,v,w,s1',s2') in order.
 !-----------------------------------------------------------------
 start_local=(/1,1,global_z_indices(1,1,myid),1/)    ! position to start reading netCDF file for this processor
 count = (/nx,ny,nz,1/)       ! 1 time slice, nx, ny, nz (local, per-processor array sizes)
 allocate(Fscratch(nx,ny,10*nz))		! pre-allocate larger space than necessary
 Fscratch(:,:,:) = 0.d0		! fill pre-allocated space

 select case(id)
 case(1)
 !------------------------------------------------------------
 ! extract variable id, check for error
 !------------------------------------------------------------
 ierr=NF_INQ_VARID(ncid,'u',varid)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR INQ_VARID -> u: ',ierr
 stop
 endif
 !------------------------------------------------------------
 ! read the corresponding variable (dimensional)
 !------------------------------------------------------------
 ierr=NF_GET_VARA_DOUBLE(ncid,varid,start_local,count,Fscratch)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> u: ',ierr
 stop
 endif
 do k=1,nz
  do j=1,ny
   do i=1,nx
   F(j,i,k) = Fscratch(i,j,k)
   enddo
  enddo 
 enddo
 write(0,*)'have read u'

 case(2)
 !------------------------------------------------------------
 ! extract variable id, check for error
 !------------------------------------------------------------
 ierr=NF_INQ_VARID(ncid,'v',varid)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR INQ_VARID -> v: ',ierr
 stop
 endif
 !------------------------------------------------------------
 ! read the corresponding variable (dimensional)
 !------------------------------------------------------------
 ierr=NF_GET_VARA_DOUBLE(ncid,varid,start_local,count,Fscratch)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> v: ',ierr
 stop
 endif
 do k=1,nz
  do j=1,ny
   do i=1,nx
   F(j,i,k) = Fscratch(i,j,k)
   enddo
  enddo 
 enddo
 write(0,*)'have read v'

 case(3)
 !------------------------------------------------------------
 ! extract variable id, check for error
 !------------------------------------------------------------
 ierr=NF_INQ_VARID(ncid,'w',varid)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR INQ_VARID -> w: ',ierr
 stop
 endif
 !------------------------------------------------------------
 ! read the corresponding variable (dimensional)
 !------------------------------------------------------------
 ierr=NF_GET_VARA_DOUBLE(ncid,varid,start_local,count,Fscratch)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> w: ',ierr
 stop
 endif
 do k=1,nz
  do j=1,ny
   do i=1,nx
   F(j,i,k) = Fscratch(i,j,k)
   enddo
  enddo 
 enddo
 write(0,*)'have read w'

 case(4)
 !------------------------------------------------------------
 ! extract variable id, check for error
 !------------------------------------------------------------
 ierr=NF_INQ_VARID(ncid,'s1',varid)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR INQ_VARID -> s1_prime: ',ierr
 stop
 endif
 !------------------------------------------------------------
 ! read the corresponding variable (dimensional)
 !------------------------------------------------------------
 ierr=NF_GET_VARA_DOUBLE(ncid,varid,start_local,count,Fscratch)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> s1_prime: ',ierr
 stop
 endif
 do k=1,nz
  do j=1,ny
   do i=1,nx
   F(j,i,k) = Fscratch(i,j,k)
   enddo
  enddo 
 enddo
 write(0,*)'have read s1_prime'

 case(5)
 !------------------------------------------------------------
 ! Currently not using second scalar s2', so just set to zero.
 !------------------------------------------------------------
 F(:,:,:)=0.d0
 write(0,*)'have read s2_prime'

 end select
 deallocate(Fscratch)

 return
end subroutine user_ics


subroutine user_forcing(x,y,z,t,FF,     &
                        id,nx,ny,nz,    &
                        call_again)
!!-----------------------------------------------------------
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
!!-----------------------------------------------------------
 use dimensional_scales,     only: velocity_scale,scalar_scale,length_scale,coriolis=>f,dgrad,rho_0,g,N2 
 use independent_variables,  only: Lx,Ly,Lz
 use decomposition_params,   only: global_z_indices
 use dependent_variables,    only: u,v,w,s1,s2
 use mpi_params,             only: myid,comm,ierr
 implicit none
 real(kind=8),dimension(:,:,:),pointer  :: phi
 integer                                :: i,j,k,m
 integer                                :: id,nx,ny,nz
 real(kind=8)                           :: x(nx),y(ny),z(nz),t
 real(kind=8)                           :: s1_bar_z(nz)
 real(kind=8)                           :: FF(ny,nx,nz) 
 logical                                :: call_again
 real(kind=8),save                      :: scale
 logical,save                           :: first_entry=.TRUE.
  
 !------------------------------------------------------------
 ! additional declarations for reading from netcdf files
 !------------------------------------------------------------
 include 'netcdf.inc'
 real(kind=8),allocatable,save	:: u_forcing(:,:),v_forcing(:,:),w_forcing(:,:),rho_forcing(:,:)
 real(kind=8),allocatable,save	:: omega_forcing(:),amp_forcing(:),phase_forcing(:),k_forcing(:),l_forcing(:)
 character(len=80)	      :: ncfile
 integer                  :: ncid,varid,forceDimID,forceLen,kDimID,kLen
  
 !------------------------------------------------------------------------------------
 ! user variables
 !------------------------------------------------------------------------------------
 ! user defined forcing parameters and variables
 !     real(kind=8)                     :: nio_amp,Twind
 !     real(kind=8)                     :: ml_depth
 !     real(kind=8)                     :: alpha,beta
 !     real(kind=8)                     :: rayleigh_z(nz),rayleigh_h(nx,ny)
 !     real(kind=8),save                      :: pi
 !-----------------------------------------------------------------------------------


 !------------------------------------------------------------
 ! read forcing variables from the netCDF file if this is the first loop.
 !------------------------------------------------------------
  if(first_entry) then
    call_again = .FALSE.  ! switch to true below if desired, decide for each id

  !------------------------------------------------------------
  ! provide netCDF filename with forcing frequencies and vertical structures.
  !------------------------------------------------------------
  ncfile = 'SaveIC_EarlyIWmodel_exp_strat_nonlin_correction.nc'
 
  !------------------------------------------------------------
  ! open the file, check for error, get some dimensions
  !------------------------------------------------------------
  ierr=NF_OPEN(ncfile,NF_NOWRITE,ncid)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) '... ERROR OPENING NETCDF FILE: user_forcing ',trim(ncfile)
 	write(0,*) '... myid,ierr ',myid,ierr
 	stop
  endif
  ! get the dimensions
  ierr = NF_INQ_DIMID(ncid, 'forcedimension', forceDimID)
  ierr = NF_INQ_DIMLEN(ncid, forceDimID, forceLen)
  ierr = NF_INQ_DIMID(ncid, 'kdimension', kDimID)
  ierr = NF_INQ_DIMLEN(ncid, kDimID, kLen) 
  ! pre-allocate space for forcing arrays
  allocate(u_forcing(kLen,forceLen))
  allocate(v_forcing(kLen,forceLen))
  allocate(w_forcing(kLen,forceLen))
  allocate(rho_forcing(kLen,forceLen))
  allocate(omega_forcing(forceLen))
  allocate(amp_forcing(forceLen))
  allocate(phase_forcing(forceLen))
  allocate(k_forcing(forceLen))
  allocate(l_forcing(forceLen))
  
  !------------------------------------------------------------
  ! read the forcing arrays
  !------------------------------------------------------------
  ! u_forcing vertical structure
  ierr=NF_INQ_VARID(ncid,'u_forcing',varid)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> u_forcing: ',ierr
  	stop
  endif
  ierr = NF_GET_VAR(ncid,varid,u_forcing)
  ierr=NF_GET_VARA_DOUBLE(ncid,varid,(/1,1/),(/kLen,forceLen/),u_forcing)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> u_forcing: ',ierr
  	stop
  endif
  ! v_forcing vertical structure
  ierr=NF_INQ_VARID(ncid,'v_forcing',varid)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> v_forcing: ',ierr
  	stop
  endif
  ierr=NF_GET_VARA_DOUBLE(ncid,varid,(/1,1/),(/kLen,forceLen/),v_forcing)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> v_forcing: ',ierr
  	stop
  endif
  ! w_forcing vertical structure
  ierr=NF_INQ_VARID(ncid,'w_forcing',varid)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> w_forcing: ',ierr
  	stop
  endif
  ierr=NF_GET_VARA_DOUBLE(ncid,varid,(/1,1/),(/kLen,forceLen/),w_forcing)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> w_forcing: ',ierr
  	stop
  endif
  ! rho_forcing vertical structure
  ierr=NF_INQ_VARID(ncid,'rho_forcing',varid)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> rho_forcing: ',ierr
  	stop
  endif
  ierr=NF_GET_VARA_DOUBLE(ncid,varid,(/1,1/),(/kLen,forceLen/),rho_forcing)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> rho_forcing: ',ierr
  	stop
  endif
 ! frequencies omega
  ierr=NF_INQ_VARID(ncid,'omega_forcing',varid)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> omega_forcing: ',ierr
  	stop
  endif
  ierr=NF_GET_VARA_DOUBLE(ncid,varid,1,forceLen,omega_forcing)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> omega_forcing: ',ierr
  	stop
  endif 
  ! amplitudes
  ierr=NF_INQ_VARID(ncid,'amp_forcing',varid)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> amp_forcing: ',ierr
  	stop
  endif
  ierr=NF_GET_VARA_DOUBLE(ncid,varid,1,forceLen,amp_forcing)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> amp_forcing: ',ierr
  	stop
  endif 
  ! random phases
  ierr=NF_INQ_VARID(ncid,'phase_forcing',varid)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> phase_forcing: ',ierr
  	stop
  endif
  ierr=NF_GET_VARA_DOUBLE(ncid,varid,1,forceLen,phase_forcing)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> phase_forcing: ',ierr
  	stop
  endif
 ! zonal wavenumber
  ierr=NF_INQ_VARID(ncid,'k_forcing',varid)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> k_forcing: ',ierr
  	stop
  endif
  ierr=NF_GET_VARA_DOUBLE(ncid,varid,1,forceLen,k_forcing)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> k_forcing: ',ierr
  	stop
  endif
 ! meridional wavenumber
  ierr=NF_INQ_VARID(ncid,'l_forcing',varid)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR INQ_VARID -> l_forcing: ',ierr
  	stop
  endif
  ierr=NF_GET_VARA_DOUBLE(ncid,varid,1,forceLen,l_forcing)
  if (ierr .ne. NF_NOERR) then
 	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> l_forcing: ',ierr
  	stop
  endif
  
 first_entry=.FALSE.
 endif

!  !------------------------------------------------------------
!  ! I'm not sure what, if anything, this does (cim).
!  !------------------------------------------------------------
!  if( id==1 ) then
!  phi=>u
!  scale=velocity_scale
!  elseif( id==2 ) then
!  phi=>v
!  scale=velocity_scale
!  elseif( id==3 ) then
!  phi=>w
!  scale=velocity_scale
!  elseif( id==4 ) then
!  phi=>s1
!  scale=scalar_scale(1)
!  elseif( id==5 ) then
!  phi=>s2
!  scale=scalar_scale(2)
!  endif

 !------------------------------------------------------------
 ! Apply forcing for each field id
 !------------------------------------------------------------    
 if( id == 1) then       !  FF <===  specify desired forcing on rhs of u equation
 	! calculate the forcing
    do i=1,nx
     do j=1,ny
      do k=1,nz
       FF(j,i,k) = 0.d0		! set FF to zero, then sum contribution from each forcing frequency.  There's gotta be a better way to do this.
       do m=1,size(omega_forcing)
        FF(j,i,k) = FF(j,i,k) + amp_forcing(m) * cos(k_forcing(m)*x(i)+l_forcing(m)*y(j)+omega_forcing(m)*t+phase_forcing(m)) * u_forcing(global_z_indices(1,1,myid)+k-1,m) 
       enddo
      enddo
     enddo 
    enddo
 call_again = .TRUE.
 
 elseif( id==2 ) then    !  FF <===  specify desired forcing on rhs of v equation
 	! calculate the forcing
    do i=1,nx
     do j=1,ny
      do k=1,nz
       FF(j,i,k) = 0.d0		! set FF to zero, then sum contribution from each forcing frequency.  There's gotta be a better way to do this.
       do m=1,size(omega_forcing)
        FF(j,i,k) = FF(j,i,k) + amp_forcing(m) * sin(k_forcing(m)*x(i)+l_forcing(m)*y(j)+omega_forcing(m)*t+phase_forcing(m)) * v_forcing(global_z_indices(1,1,myid)+k-1,m)
       enddo
      enddo
     enddo 
    enddo
 call_again = .TRUE. 
 
 elseif( id==3 ) then    !  FF <===  specify desired forcing on rhs of w equation
 	! calculate the forcing
    do i=1,nx
     do j=1,ny
      do k=1,nz
       FF(j,i,k) = 0.d0		! set FF to zero, then sum contribution from each forcing frequency.  There's gotta be a better way to do this.
       do m=1,size(omega_forcing)
        FF(j,i,k) = FF(j,i,k) + amp_forcing(m) * sin(k_forcing(m)*x(i)+l_forcing(m)*y(j)+omega_forcing(m)*t+phase_forcing(m)) * w_forcing(global_z_indices(1,1,myid)+k-1,m) 
       enddo
      enddo
     enddo 
    enddo
 call_again = .TRUE. 
 
 elseif( id == 4 ) then  !  FF <===  specify desired forcing on rhs of s1 equation
 	! calculate the forcing
    do i=1,nx
     do j=1,ny
      do k=1,nz
       FF(j,i,k) = 0.d0		! set FF to zero, then sum contribution from each forcing frequency.  There's gotta be a better way to do this.
       do m=1,size(omega_forcing)
        FF(j,i,k) = FF(j,i,k) + amp_forcing(m) * cos(k_forcing(m)*x(i)+l_forcing(m)*y(j)+omega_forcing(m)*t+phase_forcing(m)) * rho_forcing(global_z_indices(1,1,myid)+k-1,m) 
       enddo
      enddo
     enddo 
    enddo
 call_again = .TRUE. 
 
 elseif( id==5 ) then    !  FF <===  specify desired forcing on rhs of s2 equation
	! calculate the forcing
	FF(:,:,:) = 0.d0
 call_again = .FALSE.  
 
 endif

end subroutine user_forcing


subroutine particle_positions(positions,npts,Lx,Ly,Lz)  
 use user_params,   only: h,delta,d,z0,zr     !! see user_params_module.f90 for definitions
 use mpi_params,             only: myid,comm,ierr
 implicit none
 integer               :: i
 integer               :: npts                !! total number of particles
 real(kind=8)		   :: x_float(npts),y_float(npts),z_float(npts)
 real(kind=8)          :: positions(npts,3)   !! (1,2,3)->(x,y,z), in [m]
 real(kind=8)          :: Lx,Ly,Lz            !! domain size in [m]
 real(kind=8)          :: urv                 !! uniform rv in [0,1]

 !------------------------------------------------------------
 ! additional declarations for reading from netcdf files
 !------------------------------------------------------------
 include 'netcdf.inc'
 character(len=80)	      :: ncfile
 integer                  :: ncid,varid
 
 !------------------------------------------------------------
 ! provide netCDF filename with initial conditions.
 !------------------------------------------------------------
 ncfile = 'SaveIC_EarlyIWmodel_exp_strat_nonlin_correction.nc'

 !------------------------------------------------------------
 ! open the file, check for error
 !------------------------------------------------------------
 ierr=NF_OPEN(ncfile,NF_NOWRITE,ncid)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'... ERROR OPENING NETCDF FILE: particle_positions ',trim(ncfile)
	stop
 endif
 
 !------------------------------------------------------------
 ! extract variable id, check for error
 !------------------------------------------------------------
 ierr=NF_INQ_VARID(ncid,'x_float',varid)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR INQ_VARID -> x_float: ',ierr
 stop
 endif
 !------------------------------------------------------------
 ! read the corresponding variable (dimensional)
 !------------------------------------------------------------
 ierr=NF_GET_VARA_DOUBLE(ncid,varid,1,npts,x_float)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> x_float: ',ierr
 stop
 endif
 do i=1,npts
   positions(i,1) = x_float(i)
 enddo
 write(0,*)'have read x_float'
 
 !------------------------------------------------------------
 ! extract variable id, check for error
 !------------------------------------------------------------
 ierr=NF_INQ_VARID(ncid,'y_float',varid)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR INQ_VARID -> y_float: ',ierr
 stop
 endif
 !------------------------------------------------------------
 ! read the corresponding variable (dimensional)
 !------------------------------------------------------------
 ierr=NF_GET_VARA_DOUBLE(ncid,varid,1,npts,y_float)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> y_float: ',ierr
 stop
 endif
 do i=1,npts
   positions(i,2) = y_float(i)
 enddo
 write(0,*)'have read y_float'
 
 !------------------------------------------------------------
 ! extract variable id, check for error
 !------------------------------------------------------------
 ierr=NF_INQ_VARID(ncid,'z_float',varid)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR INQ_VARID -> z_float: ',ierr
 stop
 endif
 !------------------------------------------------------------
 ! read the corresponding variable (dimensional)
 !------------------------------------------------------------
 ierr=NF_GET_VARA_DOUBLE(ncid,varid,1,npts,z_float)
 if (ierr .ne. NF_NOERR) then
	write(0,*) myid,'NetCDF ERROR NF_GET_VARA_DOUBLE -> z_float: ',ierr
 stop
 endif
 do i=1,npts
   positions(i,3) = z_float(i)
 enddo
 write(0,*)'have read z_float'

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

