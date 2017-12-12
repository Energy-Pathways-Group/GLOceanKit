subroutine SetupImmersedBoundary
 use mpi_params
 use etc,                   only: logfile,runlabel
 use immersed_boundary
 use independent_variables, only: nx,ny,nz,Lx,Ly,Lz,x,y,z,t_secs                               
 use decomposition_params
 use dimensional_scales,    only: length_scale
 use methods_params
 
 
 implicit none 
 include 'mpif.h' 
 integer                  :: i,j,k,kg,n_remote,N
 integer                  :: i0,i1,k0,k1,m              ! edge indices in YBLOCK orientation
 integer                  :: locnx,locny,locnz          ! num pts in myid's YBLOCK
 real(kind=8)             :: pi
 real(kind=8)             :: dx,dy,dz                   ! dless grid spacings
 real(kind=8)             :: xi,yi,zi                   ! dless image point location
 real(kind=8)             :: x0,x1,z0,z1                ! limits on YBLOCK domain ownership
 real(kind=8)             :: d,R0,r,xc,yc,zc,W0,omega   ! simple test cases
 integer,allocatable,save :: tmp(:)
 logical,save             :: first_entry=.TRUE.
 logical                  :: smooth_surface=.TRUE.
 logical                  :: cylindrical_tank=.FALSE.
 logical                  :: flow_around_cylinder=.FALSE.
  
 !---------------------------------------------
 ! dist(y,x,z)
 ! closest distance to IB surface from the point
 ! ( x(i),y(j),z(kg) )
 ! taken positive for points in the fluid and
 ! negative for points in the neighboring solid
 !
 !   in this routine, work only w/ dless values
 !   get dimensional values from user routine
 !   and scale by length_scale
 !---------------------------------------------
 
 !---------------------------------------------
 ! nhat(y,x,z,1-3)
 ! unit normal at nearest surface point to
 !        ( x(i),y(j),z(kg) )
 ! (pointing INTO the fluid)
 !---------------------------------------------
 !     nhat(j,i,k,1) => x component
 !     nhat(j,i,k,2) => y component
 !     nhat(j,i,k,3) => z component
 !---------------------------------------------
  
 if( .not. do_immersed_boundary ) return
 if( stationary_IB .and. t_secs > 0 ) return
 pi = 4.d0*datan(1.d0)
 
 if(myid==0 .and. first_entry ) then
  write(0,*) ' ................'
  write(0,*) ' ................     hello world from SetupImmersedBoundary'
  open(1,file=logfile,position='append') 
  write(1,*) '  '
  write(1,*) '  '
  write(1,*) ' =========================================================== ' 		 
  write(1,*) ' =========================================================== ' 		 
  write(1,*) '                  SetupImmersedBoundary Report:' 		 
  write(1,*) ' =========================================================== ' 		 
  write(1,*) ' =========================================================== ' 		 
  write(1,*) '  '
 endif
 
 !-----------------------------------------------------------------------------
 !   check for special case of cylindrical tank & set params as appropriate
 !-----------------------------------------------------------------------------
  if( trim(runlabel)=='cylindrical_tank' ) then
   smooth_surface=.FALSE.           ! don't ask for user defined z_IB=h(x,y)
   cylindrical_tank=.TRUE.          ! define the sidewall surface here
   if( Lx .ne. Ly ) then
    if(myid==0) write(0,*) ' Lx should equal Ly for cylindrical tank geometry '
    stop
   endif
   R0 = .875*(Ly/2.)/length_scale    ! d'less tank radius
   if( nx > 1 ) then
    xc = (Lx/2.)/length_scale        ! dimensionless center point
   else
    xc = 0.d0
   endif
   yc = (Ly/2.)/length_scale        ! dimensionless center point
  endif
 !-----------------------------------------------------------------------------
 
 !-----------------------------------------------------------------------------
 !   check for special case of flow around cylinder & set params as appropriate
 !-----------------------------------------------------------------------------
  if( trim(runlabel)=='flow_around_cylinder' ) then
   no_slip_no_flux_IB = .FALSE.     ! disables generic IB treatment
   flow_around_cylinder = .TRUE.    ! invokes special geometrical logic here for dist,normals
   smooth_surface=.FALSE.           ! don't ask for user defined z_IB=h(x,y) here
   moving_cylinder=.TRUE.           ! flag used in the immersed boundary routine
   stationary_IB = .FALSE.          ! generic flag to see if IB needs to be recomputed each step
   
   R0 = (0.02225/2.)/length_scale   ! d'less cylinder radius
   d = 2.*R0*length_scale           ! [m]  DIMENSIONAL cylinder diameter
    !-----------------------------------------------------------------------------
    !  Browand & Winant, towed cylinder in closed tank
    !-----------------------------------------------------------------------------
     U_cyl = -.00204                  ! [m/s] speed of cylinder in y direction
     W0 = U_cyl*1.e-3                 ! [m/s] nominal vertical speed of cylinder
     omega = (U_cyl/d)*.2*2*pi        ! [1/s] frequency of cylinder oscillation
     W_cyl = W0*cos(omega*t_secs)     ! [m/s] speed of cylinder in z direction
     zc = (Lz/2.)/length_scale  &
        + ((W0/omega)/length_scale)*sin(omega*t_secs)! dimensionless cylinder center point
     yc = Ly/length_scale - 1*R0  &    ! 4*R0
        + (U_cyl*t_secs)/length_scale ! dimensionless cylinder center point
    !-----------------------------------------------------------------------------
    !  St Andrews cross -- vertically oscillating cylinder in closed tank
    !-----------------------------------------------------------------------------
    ! U_cyl = 0.d0                     ! [m/s] speed of cylinder in y direction
    ! W0 = 1.e-4                       ! [m/s] nominal vertical speed of cylinder
    ! omega = 0.07                     ! [1/s] frequency of cylinder oscillation
    ! W_cyl = W0*cos(omega*t_secs)     ! [m/s] speed of cylinder in z direction
    ! zc = (Lz/2.)/length_scale  &
    !    + ((W0/omega)/length_scale)*sin(omega*t_secs)! dimensionless cylinder center point
    ! yc = (Ly/2.)/length_scale        ! dimensionless cylinder center point
    !-----------------------------------------------------------------------------
  endif
 !-----------------------------------------------------------------------------
 
 
 dx=0.d0
 dy=0.d0
 dz=0.d0
 if(nx>1) dx = ( Lx/(nx-1.d0) )/length_scale  ! Lx/nx etc for Fourier
 if(ny>1) dy = ( Ly/(ny-1.d0) )/length_scale
 if(nz>1) dz = ( Lz/(nz-1.d0) )/length_scale
   
 sigma_d = 0.5*sqrt(dx**2 + dy**2 + dz**2)  ! Dirichlet  (No longer used anywhere...check)
 sigma_n = 1.0*sqrt(dx**2 + dy**2 + dz**2)  ! Neumann
 TOL = 25.*sigma_n                          ! don't bother w/ image vals further than TOL
 
 if( cylindrical_tank ) then         ! vertical sidewall at fixed radius
  sigma_d = 0.5*sqrt(dx**2 + dy**2)  ! Dirichlet
  sigma_n = 1.0*sqrt(dx**2 + dy**2)  ! Neumann
  TOL = 25.*sigma_n                  ! don't bother w/ image vals further than TOL
 elseif( flow_around_cylinder ) then
  sigma_d = 0.5*sqrt(dy**2 + dz**2)  ! Dirichlet
  sigma_n = 1.0*sqrt(dy**2 + dz**2)  ! Neumann
  TOL = 25.*sigma_n                  ! don't bother w/ image vals further than TOL 6
 endif
  
 !------------------------------------
 ! work here in YBLOCK orientation
 !------------------------------------
 locnx = array_size(JDIM,YBLOCK,myid)
 locny = array_size(IDIM,YBLOCK,myid)
 locnz = array_size(KDIM,YBLOCK,myid)
  
 
 !----------------------------------------------------------
 ! allocate some arrays 
 !----------------------------------------------------------
 if( first_entry ) then
  allocate( dist(locny,locnx,locnz) )        ! (y,x,z)                    
  allocate( nhat(locny,locnx,locnz,3) )      ! (y,x,z,1-3) 
                                             ! 1-3 => nx,ny,nz  (not ny,nx,nz)
  
  allocate( num_remote_image_points(numprocs) )
  num_remote_image_points(:) = 0
 endif
 
 allocate( tmp(numprocs) )
 tmp(:) = 0
  
 !--------------------------------------------------------------
 ! slightly extended domain for interpolation responsibilities
 ! only need x and z because all y vals contained in YBLOCKs
 !--------------------------------------------------------------
  i0 = global_x_indices(START,YBLOCK,myid)
  i1 = global_x_indices(END,YBLOCK,myid)
  x0 = x(i0) - dx/2.d0
  x1 = x(i1) + dx/2.d0
  
  k0 = global_z_indices(START,YBLOCK,myid)
  k1 = global_z_indices(END,YBLOCK,myid)
  z0 = z(k0) - dz/2.d0
  z1 = z(k1) + dz/2.d0
  
  

 !---------------------------------
 ! loop thru the entire YBLOCK...
 !---------------------------------
 do k=1,locnz
  kg = global_z_indices(START,YBLOCK,myid) + k - 1
  do i=1,locnx
   do j=1,locny
   

   !**********************************************************************************
   !   COMPUTE DISTANCE AND NORMAL DIRECTION...  CALL USER ROUTINE FOR IB DETAILS
   !**********************************************************************************
    
     !----------------------------------------------------
     ! use Newton minimization and user supplied h(x,y)
     ! and first and second derivatives
     ! to fill dist(:,:,:) and unit normal arrays
     !----------------------------------------------------
                                 !------------------------------------------
     if( smooth_surface ) then   ! this is the generic case z_IB = h(x,y)
                                 !------------------------------------------
     
      xpt = x(i)*length_scale    ! saved in immersed_boundary module...
      ypt = y(j)*length_scale    ! for passing around behind the scenes
      zpt = z(kg)*length_scale   ! N.B.  dimensional values
     
      call get_dist_nhat( dist(j,i,k),nhat(j,i,k,1),nhat(j,i,k,2),nhat(j,i,k,3) )                          
      dist(j,i,k)=dist(j,i,k)/length_scale  ! user routine returned dimensional values

                                 !------------------------------------------      
     elseif( cylindrical_tank ) then
                                 !------------------------------------------
                                 
      r = sqrt( (x(i)-xc)**2 + (y(j)-yc)**2 )   ! dimensionless radius
      dist(j,i,k) = R0 - r                      ! positive in fluid, neg in solid
      if( x(i)==xc .and. y(j)==yc ) then
       nhat(j,i,k,1:2) = 0.d0
      else
       nhat(j,i,k,1) = (xc-x(i))/r               ! points into interior
       nhat(j,i,k,2) = (yc-y(j))/r               ! points into interior
      endif
      nhat(j,i,k,3) = 0.d0
                                 !------------------------------------------      
     elseif( flow_around_cylinder ) then
                                 !------------------------------------------
                                 
      r = sqrt( (y(j)-yc)**2 + (z(kg)-zc)**2 )   ! dimensionless radius
      dist(j,i,k) = r - R0                       ! positive in fluid, neg in solid
      if( z(kg)==zc .and. y(j)==yc ) then
       nhat(j,i,k,2:3) = 0.d0                    ! nhat=0 at center
      else
       nhat(j,i,k,2) = (y(j)-yc)/r               ! points into fluid/outward from cylinder
       nhat(j,i,k,3) = (z(kg)-zc)/r              ! points into fluid/outward from cylinder       
      endif
      nhat(j,i,k,1) = 0.d0
                  
     endif
     
     !if( i==1 .and. k==1 ) write(0,*) y(j)*Lz,dist(j,i,k)*Lz,nhat(j,i,k,1),nhat(j,i,k,2)

!**********************************************************************************
!**********************************************************************************

     
     !---------------------------------------------------------------
     !  For all grid points in this YBLOCK outside fluid domain...
     !---------------------------------------------------------------
     if( dist(j,i,k) < 0 ) then   
     
      !----------------------------------------------------------
      !  Decide if corresponding image point is local or remote.
      !  Count the number of points on this YBLOCK
      !  where BOTH image point is nonlocal AND the
      !  distance to the IB is less than or equal to TOL
      !----------------------------------------------------------
      d = abs( dist(j,i,k) )
      xi = x(i)  + 2.d0*d*nhat(j,i,k,1)   ! x image location
      zi = z(kg) + 2.d0*d*nhat(j,i,k,3)   ! z image location
     
      if( nx > 1 .and. xi < x0 .and. abs(dist(j,i,k)) <= TOL ) then
       tmp(myid+1)=tmp(myid+1) + 1
       goto 997   ! location flagged as having nonlocal image location
      endif
     
      if( nx > 1 .and. xi >= x1 .and. abs(dist(j,i,k)) <= TOL ) then
       tmp(myid+1)=tmp(myid+1) + 1
       goto 997   ! location flagged as having nonlocal image location
      endif
     
      if( zi < z0 .and. abs(dist(j,i,k)) <= TOL ) then
       !write(0,*) 'image point below z0 ',zi,z0,z1
       tmp(myid+1)=tmp(myid+1) + 1
       goto 997   ! location flagged as having nonlocal image location
      endif
     
      if( zi >= z1 .and. abs(dist(j,i,k)) <= TOL ) then
       !write(0,*) 'image above z1 ',zi,z0,z1
       tmp(myid+1)=tmp(myid+1) + 1
       goto 997   ! location flagged as having nonlocal image location
      endif
     
997   continue  
     endif

     
   enddo
  enddo
 enddo
 

 !-------------------------------------------------------------------------- 
 !  do a global sum so all processors know the number of points
 !   => num_remote_image_points(pid+1) contains the number of points
 !      outside the fluid domain but close enough to the IB surface
 !      to matter on processor pid
 !--------------------------------------------------------------------------
 call mpi_allreduce(tmp,num_remote_image_points,numprocs,mpi_integer,mpi_sum,comm,ierr )
 deallocate( tmp )
 
 !--------------------------------------------------------------------------
 !   need the total number across all processors to allocate the rest 
 !   of the arrays we'll need...
 !--------------------------------------------------------------------------
 total_num_remote_image_points = SUM( num_remote_image_points(:) )
 N = total_num_remote_image_points
  
 
 if( N > 0 ) then
  if( allocated( xyz_solid ) ) then
   deallocate(xyz_solid,xyz_image,ijk_solid,image_vals,owner_solid,owner_image,xtmp,itmp,fval)
  endif
  
  allocate( xyz_solid(N,3) )    !  x,y,z for solid point near IB with remote image point
  allocate( xyz_image(N,3) )    !  x,y,z for corresponding image point
  allocate( ijk_solid(N,3) )    !  local i,j,k indices for processor that owns the solid point
  allocate( image_vals(N) )     !  array to store interpolated scalar values at image point
  allocate( owner_solid(N) )    !  pid of the owner of the solid point
  allocate( owner_image(N) )    !  pid of the owner of the image point
  allocate( fval(N) )           !  array for storing interpolated function values
  allocate( xtmp(N) )           !  tmp array for globally summing real(kind=8) arrays
  allocate( itmp(N) )           !  tmp array for globally summing integer arrays
  
  xyz_solid(:,:)=0.d0
  xyz_image(:,:)=0.d0
  ijk_solid(:,:)=0
  owner_solid(:)=0
  owner_image(:)=0
  
  !-----------------------------------------------------------
  ! Find out the starting index for myid in these arrays.
  ! It depends on how many special pts exist for pids < myid.
  !----------------------------------------------------------- 
   my_start_index=0
   do i=0,myid-1
    my_start_index = my_start_index + num_remote_image_points(i+1)
   enddo
   my_start_index = my_start_index + 1
  
  !-----------------------------------------------------------
  ! If there are any special points residing on processor myid,
  ! loop through entire myid YBLOCK and fill in known values.
  !-----------------------------------------------------------
  if( num_remote_image_points(myid) > 0 ) then
  
  m = my_start_index
  do k=1,locnz
   kg = global_z_indices(START,YBLOCK,myid) + k - 1
    do i=1,locnx
     do j=1,locny
     
     !--------------------------------------------------------------
     ! determine...(again) whether (x,y,z) is a "special" point...
     !--------------------------------------------------------------
     
     if( dist(j,i,k) < 0 ) then   ! pt is outside fluid domain
     
      xi = x(i)  + 2.d0*d*nhat(j,i,k,1)   ! x image location
      yi = y(j)  + 2.d0*d*nhat(j,i,k,2)   ! y image location
      zi = z(kg) + 2.d0*d*nhat(j,i,k,3)   ! z image location
     
      if( nx > 1 .and. xi < x0 .and. abs(dist(j,i,k)) <= TOL ) then
       goto 998   ! location flagged as having nonlocal image location
      endif
     
      if( nx > 1 .and. xi >= x1 .and. abs(dist(j,i,k)) <= TOL ) then
       goto 998   ! location flagged as having nonlocal image location
      endif
     
      if( zi < z0 .and. abs(dist(j,i,k)) <= TOL ) then
       goto 998   ! location flagged as having nonlocal image location
      endif
     
      if( zi >= z1 .and. abs(dist(j,i,k)) <= TOL ) then
       goto 998   ! location flagged as having nonlocal image location
      endif
      
      goto 999   ! solid but not a special point, don't need to save info...
     
998   continue    !  xyz point is in solid, close to bdry, owned by myid w/ remote image
      xyz_solid(m,1) = x(i)
      xyz_solid(m,2) = y(j)
      xyz_solid(m,3) = z(kg)
      
      ijk_solid(m,1) = i
      ijk_solid(m,2) = j
      ijk_solid(m,3) = k        ! local YBLOCK index
      owner_solid(m) = myid
      
      xyz_image(m,1) = xi
      xyz_image(m,2) = yi
      xyz_image(m,3) = zi
      
      m = m + 1                ! next row in arrays...

999   continue
     endif        ! end dist < 0 block,  go to next xyz point
          
     enddo
    enddo
   enddo
   
   endif   ! end of logic for processors with "special" points
   
   !--------------------------------------------------------------------
   ! Globally sum the results so all processors have the information.
   ! Global sums are written to tmp vectors, copy back where needed.
   !--------------------------------------------------------------------
   do j=1,3
    call mpi_allreduce(xyz_solid(:,j),xtmp,N,mpi_double_precision,mpi_sum,comm,ierr)
     xyz_solid(:,j)=xtmp 
    call mpi_allreduce(xyz_image(:,j),xtmp,N,mpi_double_precision,mpi_sum,comm,ierr)
     xyz_image(:,j)=xtmp
     
    call mpi_allreduce(ijk_solid(:,j),itmp,N,mpi_integer,mpi_sum,comm,ierr)
     ijk_solid(:,j)=xtmp
   enddo
   
   call mpi_allreduce(owner_solid(:),itmp,N,mpi_integer,mpi_sum,comm,ierr)
     owner_solid(:)=itmp
        
   
   !--------------------------------------------------------------------
   ! determine if myid owns any of the image points needed by other pids
   !--------------------------------------------------------------------
   do j=1,N
    xi = xyz_image(j,1)
    zi = xyz_image(j,3)
    !------------------------------------------
    ! if xi, zi are in myid's extended domain
    !------------------------------------------
    if( xi >= x0 .and. xi < x1 ) then
     if( zi >= z0 .and. zi < z1 ) then
      owner_image(j) = myid   ! myid owns the image point
     endif
    endif    
   enddo
   
   !--------------------------------------------------------------------
   ! make image point owners globally available
   !--------------------------------------------------------------------
   call mpi_allreduce(owner_image(:),itmp,N,mpi_integer,mpi_sum,comm,ierr )
     owner_image(:)=itmp
    
   !---------------------------------------------------------------------------------------------
   !  Now, rows [my_start_index] through [num_remote_image_points(myid+1)-1]
   !  contain solid points close to IB that live on my processor but have remote images
   
   !  For all rows j such that owner_image(j)==myid, processor myid has to do the interpolation
   
   !  Once the interpolations are done and globally summed, the values in image_vals
   !  in rows [my_start_index] through [num_remote_image_points(myid+1)-1] can be
   !  used by processor myid to do it's IB calculations
   
   !  So... each time step, each processor finds the s1 interpolating polynomial
   !  (it already does this...)
   !  and then goes through these arrays and does an eval when it owns the image location
   !  all the interpolated values at image pts are globally shared
   !---------------------------------------------------------------------------------------------
   
 endif   ! end N > 0 block
 
 
if(myid==0 .and. first_entry) then 
 write(0,*) ' ................ Total number of remote image points: ',N
 write(1,*) ' ................ Total number of remote image points: ',N
 write(1,*) ' -----> SetupImmersedBoundary routine exiting normally  <---------- ' 		 
 close(1)
endif    
 
 call mpi_barrier(comm,ierr)

first_entry = .FALSE.
return
end subroutine SetupImmersedBoundary



subroutine get_dist_nhat(d,nhat_x,nhat_y,nhat_z)
!--------------------------------------------------------------------------------
!   x0,y0,z0 = position in space,  dless on input , in immersed bdry module
!   we seek the nearest surface pt
!   work w/ DIMENSIONAL distance as we interface to a user routine
!   scale d value immediately after return to calling routine
!   NB x0,y0,z0 passed in via common for use w/  SNSQE solver
!--------------------------------------------------------------------------------
 use decomposition_params
 use mpi_params,               only: myid
 use independent_variables,    only: x,y,z
 use dimensional_scales,       only: length_scale
 use immersed_boundary,        only: xpt,ypt,zpt
 use methods_params                                ! special IB case flags
 
 implicit none
 real(kind=8)             :: d,nhat_x,nhat_y,nhat_z
 integer                  :: locnx,locny,i,j,imin,jmin,isign
 real(kind=8)             :: test,d2min,xx,yy,zz,x0,y0,z0
 logical,save             :: first_entry=.TRUE.
 
 real(kind=8)             :: XY(2),F,H(2,2),G(2)
 logical,save             :: PRNT=.FALSE.
 real(kind=8),save        :: R1,R2,GAMMA,BETA,EPS
 integer,save             :: NVARS,LIMIT,IFLAG,ITER
 real(kind=8)             :: hx,hy,hxx,hyy,hxy
 
 if( first_entry) then  ! set default params for newton solver
  NVARS = 2    ! x and y
  R1=0.00001
  R2=0.00005
  GAMMA=2.0
  BETA=0.5
  EPS=0.000001
  LIMIT=49
  if(myid==0) then
   if( hc_IB_test ) then
    write(0,*) '.............. Setting up IB for special case hc_IB_test'
   endif
  endif
  first_entry=.FALSE.
 endif
 
 !---------------------------------------------------------------
 !  just a change in notation, dimensional quantities...
 !  fixed position in space for which we want min dist and normal
 !---------------------------------------------------------------
  x0 = xpt
  y0 = ypt
  z0 = zpt    ! [m] fixed position in domain...
  
 !----------------------------------------------------
 ! surface height directly above or below (x0,y0,z0)
 ! to determine if point is in fluid domain
 !----------------------------------------------------
  call surface_height_derivs(x0,y0,zz,hx,hy,hxx,hyy,hxy)
  if( z0 > zz ) then
   isign = 1  ! (x0,y0,z0) in fluid region
  else
   isign = -1 ! (x0,y0,z0) not in fluid region
  endif

 
 
 !-----------------------------------------------------------
 !  find a good starting point...
 !-----------------------------------------------------------
  locnx = array_size(JDIM,YBLOCK,myid)   ! 
  locny = array_size(IDIM,YBLOCK,myid)   ! 
  d2min=1.e16
  do j=1,locny
   do i=1,locnx
   
    xx = x(i)*length_scale  ! [m]
    yy = y(j)*length_scale  ! [m]
    !call surface_height(xx,yy,zz)
    call surface_height_derivs(xx,yy,zz,hx,hy,hxx,hyy,hxy)
    test = (x0-xx)**2 + (y0-yy)**2 + (z0-zz)**2   ! dist**2
    if( test < d2min ) then
     d2min = test
     imin=i
     jmin=j
    endif
   
   enddo
  enddo  
  !----------------------------------------------------------
  !  x(imin) and y(jmin) should now be good starting guesses
  !  solve the minimization problem ....
  !----------------------------------------------------------  
  XY(1) = x(imin)*length_scale  ! [m]
  XY(2) = y(jmin)*length_scale  ! [m]
  call NEWTON(NVARS,XY,F,G,H,ITER,IFLAG,R1,R2,GAMMA,BETA,EPS,LIMIT,PRNT) 
  
  !----------------------------------------------------------
  !  x and y positions returned in XY(:)
  !  get corresponding z=h(x,y)
  !----------------------------------------------------------
  xx = XY(1)
  yy = XY(2)
  call surface_height_derivs(xx,yy,zz,hx,hy,hxx,hyy,hxy)
      
  
  !----------------------------------------------------------
  !  now compute signed distance and unit normal
  !  pointing into the fluid
  !----------------------------------------------------------
  d = sqrt( (x0-xx)**2 + (y0-yy)**2  + (z0-zz)**2 )
  if( isign > 0 ) then      ! (x0,y0,z0) in fluid region
   nhat_x = (x0-xx)/d
   nhat_y = (y0-yy)/d
   nhat_z = (z0-zz)/d
  elseif(isign < 0 ) then   ! (x0,y0,z0) in fluid region
   nhat_x = (xx-x0)/d
   nhat_y = (yy-y0)/d
   nhat_z = (zz-z0)/d
   d = isign*d      
  endif
  
  
  
  !------------------------------------------------------------
  !  this logic assumed that ABOVE h(x,y) was fluid and
  !  BELOW was solid... not always true...
  !------------------------------------------------------------
     if( hc_IB_test ) then  ! convention is reversed
      d = -d
      nhat_x = -nhat_x
      nhat_y = -nhat_y
      nhat_z = -nhat_z      
     endif
  !------------------------------------------------------------
  
 return
end subroutine get_dist_nhat



function FUNCT(n,xy)
use immersed_boundary,        only: xpt,ypt,zpt
implicit none
integer                          :: n
real(kind=8)                     :: h,xy(n),x,y,FUNCT
real(kind=8)                     :: hx,hy,hxx,hyy,hxy

 x=xy(1)
 y=xy(2)
 call surface_height_derivs(x,y,h,hx,hy,hxx,hyy,hxy)
 
 !-----------------------
 !  square of distance
 !-----------------------
 FUNCT = (x-xpt)**2 + (y-ypt)**2 + (h-zpt)**2

 return
end


subroutine DERIVS(n,xy,G,HH)
use immersed_boundary,        only: xpt,ypt,zpt
implicit none
integer                          :: n,i,j
real(kind=8)                     :: xy(n),G(n),HH(n,n)
real(kind=8)                     :: x,y,h,hx,hy,hxx,hyy,hxy

 x=xy(1)
 y=xy(2)
 call surface_height_derivs(x,y,h,hx,hy,hxx,hyy,hxy)
 
 ! f = distance^2
 ! function derivs G=[fx,fy]
 G(1) = 2.d0*(x-xpt) + 2.d0*(h-zpt)*hx
 G(2) = 2.d0*(y-ypt) + 2.d0*(h-zpt)*hy
 
 HH(1,1) = 2.d0 + 2.d0*(h-zpt)*hxx + 2.d0*hx**2    ! fxx
 HH(2,2) = 2.d0 + 2.d0*(h-zpt)*hyy + 2.d0*hy**2    ! fyy
 HH(1,2) = 2.d0*(h-zpt)*hxy + 2.d0*hx*hy           ! fxy
 HH(2,1) = HH(1,2)                                 ! fyx

 return
end











!subroutine unit_normal(x,n)
 ! input x=(x,y) are dimensional values
 ! eval_surface is a user routine that
 ! expects dimensional input/outputs
! implicit none
! real(kind=8) :: x(2),n(3),h,h_x,h_y
! real(kind=8) :: mag
 
!   call eval_surface(x,h,h_x,h_y)
!   mag = sqrt( h_x**2 + h_y**2 + 1.d0 )
!   n(1) = -h_x/mag
!   n(2) = -h_y/mag
!   n(3) = 1.d0/mag
   
 !return
!end subroutine unit_normal  ! end of function unit_normal



!subroutine unit_tangents(x,t1,t2)
 ! input x=(x,y) are dimensional values
 ! eval_surface is a user routine that
 ! expects dimensional input/outputs
! implicit none
! real(kind=8) :: x(2),t1(3),t2(3),h,h_x,h_y
! real(kind=8) :: n(3),mag,xx
 
!   call eval_surface(x,h,h_x,h_y)
!   mag = sqrt( h_x**2.d0 + h_y**2.d0 + 1.d0 )
!   n(1) = -h_x/mag
!   n(2) = -h_y/mag
!   n(3) = 1.d0/mag
   
   ! t1 is orthogonal to n but constrained
   ! so that the y component is equal to zero
!   xx = sqrt( n(1)**2 + n(3)**2 )
!   t1(1) = n(3)/xx
!   t1(2) = 0.d0
!   t1(3) = -n(1)/xx
   
   ! t2 = n x t1
!   t2(1) = -n(2)*n(1)/xx
!   t2(2) = ( n(3)**2 + n(1)**2 )/xx
!   t2(3) = -n(2)*n(3)/xx
   
   ! make sure t2 is exactly normalized
!   mag = sqrt( t2(1)**2 + t2(2)**2 + t2(3)**2 )
!   t2(:) = t2(:)/mag
      
! return
!end subroutine unit_tangents  ! end of function unit_tangents
