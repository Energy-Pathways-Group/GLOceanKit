subroutine SetupInterpolation
 use mpi_params
 use decomposition_params
 use interpolation
 use independent_variables,    only: x,y,z,          &
                                     xdim_periodic,  &
                                     ydim_periodic,  &
                                     zdim_periodic
 use etc,  only: logfile
 
 
 implicit none
 integer            :: i0,nwrk
 integer            :: nx,ny,nz

 if(myid==0) then
  write(0,*) ' ................'
  write(0,*) ' ................     hello world from SetupInterpolation'
  open(1,file=logfile,position='append') 
  write(1,*) '  '
  write(1,*) '  '
  write(1,*) ' =========================================================== ' 		 
  write(1,*) ' =========================================================== ' 		 
  write(1,*) '                  SetupInterpolation Report:' 		 
  write(1,*) ' =========================================================== ' 		 
  write(1,*) ' =========================================================== ' 		 
  write(1,*) '  '
 endif
 
 !------------------------------------------------------------
 ! XBLOCK
 !------------------------------------------------------------
   nx = array_size(IDIM,XBLOCK,myid)   ! x is dim 1 in XBLOCK
   ny = array_size(JDIM,XBLOCK,myid)   ! y is dim 2 in XBLOCK
   nz = array_size(KDIM,XBLOCK,myid)   ! z is dim 3 in XBLOCK
  
   ! method order in (x,y,z) <-> (1,2,3)  directions
   interp_order_XBLOCK(1:3)=3+1   !! cubic spline
 
   if( interp_order_XBLOCK(1) .ge. nx ) interp_order_XBLOCK(1)=nx-1
   if( interp_order_XBLOCK(2) .ge. ny ) interp_order_XBLOCK(2)=ny-1
   if( interp_order_XBLOCK(3) .ge. nz ) interp_order_XBLOCK(3)=nz-1
  
   nwrk = nx*ny*nz + 2*max(interp_order_XBLOCK(1)*(nx+1),   &
                           interp_order_XBLOCK(2)*(ny+1),   &
                           interp_order_XBLOCK(3)*(nz+1) )

   allocate( wrk_XBLOCK(nwrk) )
   allocate( tx_XBLOCK(nx+interp_order_XBLOCK(1))  )
   allocate( ty_XBLOCK(ny+interp_order_XBLOCK(2))  )
   allocate( tz_XBLOCK(nz+interp_order_XBLOCK(3))  )
   
   !-----------------------------------------
   !  independent variable is x, use generic
   !  knot strategy for global x dimension, 
   !  unless x is periodic
   !-----------------------------------------
   if( nx > 1 ) then
    if( xdim_periodic ) then
     call dbknot_extrap_both(x,nx,interp_order_XBLOCK(1),tx_XBLOCK)
    else
     call dbknot(x,nx,interp_order_XBLOCK(1),tx_XBLOCK)
    endif     
   endif
   
   !-----------------------------------------
   !  use exterior knots in y to extrapolate
   !  when pt is slightly outside myid's range
   !-----------------------------------------
   if( ny > 1 ) then
    i0=global_y_indices(START,XBLOCK,myid)
    if( p1 > 1 .or. ydim_periodic ) then
     call dbknot_extrap_both(y(i0),ny,interp_order_XBLOCK(2),ty_XBLOCK)
    else
     call dbknot(y(i0),ny,interp_order_XBLOCK(2),ty_XBLOCK)
    endif
   endif
   
   !-----------------------------------------
   !  use exterior knots in z to extrapolate
   !  when pt is slightly outside myid's range
   !-----------------------------------------
   if( nz > 1 ) then
    i0=global_z_indices(START,XBLOCK,myid)
    if( p2 > 1 ) then
     call dbknot_extrap_both(z(i0),nz,interp_order_XBLOCK(3),tz_XBLOCK)
    else
     call dbknot(z(i0),nz,interp_order_XBLOCK(3),tz_XBLOCK)
    endif
   endif
   
   
 !------------------------------------------------------------
 ! YBLOCK
 !------------------------------------------------------------
   nx = array_size(JDIM,YBLOCK,myid)   ! x is dim 2 in YBLOCK
   ny = array_size(IDIM,YBLOCK,myid)   ! y is dim 1 in YBLOCK
   nz = array_size(KDIM,YBLOCK,myid)   ! z is dim 3 in YBLOCK
  
   ! method order in (x,y,z) <-> (1,2,3)  directions
   interp_order_YBLOCK(1:3)=3+1   !! cubic spline
 
   if( interp_order_YBLOCK(1) .ge. nx ) interp_order_YBLOCK(1)=nx-1
   if( interp_order_YBLOCK(2) .ge. ny ) interp_order_YBLOCK(2)=ny-1
   if( interp_order_YBLOCK(3) .ge. nz ) interp_order_YBLOCK(3)=nz-1
  
   nwrk = nx*ny*nz + 2*max(interp_order_YBLOCK(1)*(nx+1),   &
                           interp_order_YBLOCK(2)*(ny+1),   &
                           interp_order_YBLOCK(3)*(nz+1) )

   allocate( wrk_YBLOCK(nwrk) )
   allocate( tx_YBLOCK(nx+interp_order_YBLOCK(1))  )
   allocate( ty_YBLOCK(ny+interp_order_YBLOCK(2))  )
   allocate( tz_YBLOCK(nz+interp_order_YBLOCK(3))  )
   
   !-----------------------------------------
   !  independent variable is y, use generic
   !  knot strategy for global y dimension, 
   !  unless y is periodic
   !-----------------------------------------
   if( ny > 1 ) then
    if( ydim_periodic ) then
     call dbknot_extrap_both(y,ny,interp_order_YBLOCK(2),ty_YBLOCK)
    else
     call dbknot(y,ny,interp_order_YBLOCK(2),ty_YBLOCK)
    endif
   endif
   
   !-----------------------------------------
   !  use exterior knots in x to extrapolate
   !  when pt is slightly outside myid's range
   !-----------------------------------------
   if( nx > 1 ) then
    i0=global_x_indices(START,YBLOCK,myid)
    if( p1 > 1 .or. xdim_periodic ) then
     call dbknot_extrap_both(x(i0),nx,interp_order_YBLOCK(1),tx_YBLOCK)
    else
     call dbknot(x(i0),nx,interp_order_YBLOCK(1),tx_YBLOCK)
    endif
   endif
   
   !-----------------------------------------
   !  use exterior knots in z to extrapolate
   !  when pt is slightly outside myid's range
   !-----------------------------------------
   if( nz > 1 ) then
    i0=global_z_indices(START,YBLOCK,myid)
    if( p2 > 1 ) then
     call dbknot_extrap_both(z(i0),nz,interp_order_YBLOCK(3),tz_YBLOCK)
    else
     call dbknot(z(i0),nz,interp_order_YBLOCK(3),tz_YBLOCK)
    endif
   endif
 
 
 !------------------------------------------------------------
 ! ZBLOCK
 !------------------------------------------------------------
   nx = array_size(JDIM,ZBLOCK,myid)   ! x is dim 2 in ZBLOCK
   ny = array_size(KDIM,ZBLOCK,myid)   ! y is dim 3 in ZBLOCK
   nz = array_size(IDIM,ZBLOCK,myid)   ! z is dim 1 in ZBLOCK
  
   ! method order in (x,y,z) <-> (1,2,3)  directions
   interp_order_ZBLOCK(1:3)=3+1   !! cubic spline
 
   if( interp_order_ZBLOCK(1) .ge. nx ) interp_order_ZBLOCK(1)=nx-1
   if( interp_order_ZBLOCK(2) .ge. ny ) interp_order_ZBLOCK(2)=ny-1
   if( interp_order_ZBLOCK(3) .ge. nz ) interp_order_ZBLOCK(3)=nz-1
  
   nwrk = nx*ny*nz + 2*max(interp_order_ZBLOCK(1)*(nx+1),   &
                           interp_order_ZBLOCK(2)*(ny+1),   &
                           interp_order_ZBLOCK(3)*(nz+1) )

   allocate( wrk_ZBLOCK(nwrk) )
   allocate( tx_ZBLOCK(nx+interp_order_ZBLOCK(1))  )
   allocate( ty_ZBLOCK(ny+interp_order_ZBLOCK(2))  )
   allocate( tz_ZBLOCK(nz+interp_order_ZBLOCK(3))  )
   
   !-----------------------------------------
   !  independent variable is z, use generic
   !  knot strategy for global z dimension, 
   !  unless z is periodic
   !-----------------------------------------
   if( nz > 1 ) then
    if( zdim_periodic ) then
     call dbknot_extrap_both(z,nz,interp_order_ZBLOCK(3),tz_ZBLOCK)
    else
     call dbknot(z,nz,interp_order_ZBLOCK(3),tz_ZBLOCK)
    endif
   endif
   
   !-----------------------------------------
   !  use exterior knots in x to extrapolate
   !  when pt is slightly outside myid's range
   !-----------------------------------------
   if( nx > 1 ) then
    i0=global_x_indices(START,ZBLOCK,myid)
    if( p1 > 1 .or. xdim_periodic ) then
     call dbknot_extrap_both(x(i0),nx,interp_order_ZBLOCK(1),tx_ZBLOCK)
    else
     call dbknot(x(i0),nx,interp_order_ZBLOCK(1),tx_ZBLOCK)
    endif
   endif
   
   !-----------------------------------------
   !  use exterior knots in y to extrapolate
   !  when pt is slightly outside myid's range
   !-----------------------------------------
   if( ny > 1 ) then
    i0=global_y_indices(START,ZBLOCK,myid)
    if( p2 > 1 .or. ydim_periodic ) then
     call dbknot_extrap_both(y(i0),ny,interp_order_ZBLOCK(2),ty_ZBLOCK)
    else
     call dbknot(y(i0),ny,interp_order_ZBLOCK(2),ty_ZBLOCK)
    endif
   endif
 
 if( myid==0 ) then
  write(1,*) ' -----> SetupInterpolation routine exiting normally  <---------- ' 		 
  close(1)
 endif
 
end subroutine SetupInterpolation


subroutine Spline_Interp_XB(F)
 use decomposition_params 
 use mpi_params
 use interpolation,               order=>interp_order_XBLOCK,  &
                                  tx => tx_XBLOCK,             &
                                  ty => ty_XBLOCK,             &
                                  tz => tz_XBLOCK,             &
                                  wrk => wrk_XBLOCK
                                 
 use independent_variables, only: x,y,z
  
 use intermediate_variables,only: e3 => tmpX
                            
 implicit none
 real(kind=8)                                :: F(*)
 integer                                     :: iflag 
 integer,save                                :: istart,jstart,kstart
 integer,save                                :: nx,ny,nz
 logical,save                                :: first_entry=.TRUE.
 
 if( first_entry ) then
  istart = global_x_indices(START,XBLOCK,myid)
  jstart = global_y_indices(START,XBLOCK,myid)
  kstart = global_z_indices(START,XBLOCK,myid)
  nx = array_size(IDIM,XBLOCK,myid)   ! x is dim 1 in XBLOCK
  ny = array_size(JDIM,XBLOCK,myid)   ! y is dim 2 in XBLOCK
  nz = array_size(KDIM,XBLOCK,myid)   ! z is dim 3 in XBLOCK
  first_entry=.FALSE.
 endif

 !!========================================================
 !! Compute coefficients of the 3D tensor basis functions
 !! using deBoor routine from nist/gams/cmlib
 !!========================================================  
 
  
 !-------------------------------
 !  3d  nx,ny > 1
 !-------------------------------
 if( nx > 1 .and. ny > 1 ) then
 
  iflag = 1   ! knots precomputed
  call db3ink(x(istart),nx,                   &   ! XBLOCK format
              y(jstart),ny,                   &
              z(kstart),nz,                   &
              F,nx,ny,                        &
              order(1),order(2),order(3),     &
              tx(1),ty(1),tz(1),              &
              e3,wrk,iflag)
              
   if(iflag > 1 ) then
     write(0,*) ' myid ',myid
    write(0,*) 'error condition detected by db3ink:',iflag
    write(0,*) 'nx,ny,nz,order ',nx,ny,nz,order(1),order(2),order(3)
    write(0,*) ' x  '
    write(0,*) x 
    write(0,*) ' y  '
    write(0,*) y 
    write(0,*) ' z  '
    write(0,*) z 
    write(0,*) ' tx '
    write(0,*) tx
    write(0,*) ' ty '
    write(0,*) ty
    write(0,*) ' tz '
    write(0,*) tz
    stop
   endif
  
  !-------------------------------
  !  2d  ny=1
  !-------------------------------
  elseif( nx > 1 .and. ny == 1 ) then
    
   iflag = 1   ! knots precomputed
   call db2ink(x (istart),nx,              &
               z (kstart),nz,              &
               F,nx,order(1),order(3),     &
               tx(1),tz(1),e3,wrk,iflag)
               
   if(iflag > 1 ) then
     write(0,*) ' myid ',myid
    write(0,*) 'error condition detected by db2ink:',iflag
    write(0,*) 'nx,nz,order ',nx,nz,order(1),order(3)
    write(0,*) ' x  '
    write(0,*) x 
    write(0,*) ' z  '
    write(0,*) z 
    write(0,*) ' tx '
    write(0,*) tx
    write(0,*) ' tz '
    write(0,*) tz
    stop
   endif
  
  !-------------------------------
  !  2d  nx=1
  !-------------------------------
  elseif( nx == 1 .and. ny > 1 ) then
  
   iflag = 1   ! knots precomputed
   call db2ink(y (jstart),ny,              &
               z (kstart),nz,              &
               F,ny,order(2),order(3),     &
               ty(1),tz(1),e3,wrk,iflag)
               
   if(iflag > 1 ) then
    write(0,*) ' myid ',myid
    write(0,*) 'error condition detected by db2ink:',iflag
    write(0,*) 'ny,nz,order ',ny,nz,order(2),order(3)
    write(0,*) ' y  '
    write(0,*) y 
    write(0,*) ' z  '
    write(0,*) z 
    write(0,*) ' ty '
    write(0,*) ty
    write(0,*) ' tz '
    write(0,*) tz
    stop
   endif
  
  endif
end subroutine Spline_Interp_XB


subroutine Spline_Interp_YB(F)
 use decomposition_params 
 use mpi_params
 use interpolation,               order=>interp_order_YBLOCK,  &
                                  tx => tx_YBLOCK,             &
                                  ty => ty_YBLOCK,             &
                                  tz => tz_YBLOCK,             &
                                  wrk => wrk_YBLOCK
                                 
 use independent_variables, only: x,y,z

 use intermediate_variables,only: e3 => tmpY
                            
 implicit none
 real(kind=8)                                :: F(*)
 integer                                     :: iflag
 integer,save                                :: istart,jstart,kstart
 integer,save                                :: nx,ny,nz
 logical,save                                :: first_entry=.TRUE.
 
 if( first_entry ) then
  istart = global_x_indices(START,YBLOCK,myid)
  jstart = global_y_indices(START,YBLOCK,myid)
  kstart = global_z_indices(START,YBLOCK,myid)
  nx = array_size(JDIM,YBLOCK,myid)   ! x is dim 2 in YBLOCK
  ny = array_size(IDIM,YBLOCK,myid)   ! y is dim 1 in YBLOCK
  nz = array_size(KDIM,YBLOCK,myid)   ! z is dim 3 in YBLOCK
  first_entry=.FALSE.
 endif

 !!========================================================
 !! Compute coefficients of the 3D tensor basis functions
 !! using deBoor routine from nist/gams/cmlib
 !!========================================================  
 
  
 !-------------------------------
 !  3d  nx,ny > 1
 !-------------------------------
 if( nx > 1 .and. ny > 1 ) then
 
  iflag = 1   ! knots precomputed
  call db3ink(y (jstart),ny,                  &   ! YBLOCK format
              x (istart),nx,                  &
              z (kstart),nz,                  &
              F,ny,nx,                        &
              order(2),order(1),order(3),     &
              ty(1),tx(1),tz(1),              &
              e3,wrk,iflag)
              
   if(iflag > 1 ) then
     write(0,*) ' myid ',myid
    write(0,*) 'error condition detected by db3ink:',iflag
    write(0,*) 'nx,ny,nz,order ',nx,ny,nz,order(1),order(2),order(3)
    write(0,*) ' x  '
    write(0,*) x 
    write(0,*) ' y  '
    write(0,*) y 
    write(0,*) ' z  '
    write(0,*) z 
    write(0,*) ' tx '
    write(0,*) tx
    write(0,*) ' ty '
    write(0,*) ty
    write(0,*) ' tz '
    write(0,*) tz
    stop
   endif
  
  !-------------------------------
  !  2d  ny=1
  !-------------------------------
  elseif( nx > 1 .and. ny == 1 ) then
    
   iflag = 1   ! knots precomputed
   call db2ink(x(istart),nx,               &
               z(kstart),nz,               &
               F,nx,order(1),order(3),     &
               tx(1),tz(1),e3,wrk,iflag)
               
   if(iflag > 1 ) then
     write(0,*) ' myid ',myid
    write(0,*) 'error condition detected by db2ink:',iflag
    write(0,*) 'nx,nz,order ',nx,nz,order(1),order(3)
    write(0,*) ' x '
    write(0,*) x
    write(0,*) ' z '
    write(0,*) z
    write(0,*) ' tx '
    write(0,*) tx
    write(0,*) ' tz '
    write(0,*) tz
    stop
   endif
  
  !-------------------------------
  !  2d  nx=1
  !-------------------------------
  elseif( nx == 1 .and. ny > 1 ) then
  
   iflag = 1   ! knots precomputed
   call db2ink(y(jstart),ny,               &
               z(kstart),nz,               &
               F,ny,order(2),order(3),     &
               ty(1),tz(1),e3,wrk,iflag)
               
   if(iflag > 1 ) then
    write(0,*) ' myid ',myid
    write(0,*) 'error condition detected by db2ink:',iflag
    write(0,*) 'ny,nz,order ',ny,nz,order(2),order(3)
    write(0,*) ' y '
    write(0,*) y
    write(0,*) ' z '
    write(0,*) z
    write(0,*) ' ty '
    write(0,*) ty
    write(0,*) ' tz '
    write(0,*) tz
    stop
   endif
  
  endif
end subroutine Spline_Interp_YB


subroutine Spline_Interp_ZB(F)
 use decomposition_params 
 use mpi_params
 use interpolation,               order=>interp_order_ZBLOCK,  &
                                  tx => tx_ZBLOCK,             &
                                  ty => ty_ZBLOCK,             &
                                  tz => tz_ZBLOCK,             &
                                  wrk => wrk_ZBLOCK
                                 
 use independent_variables, only: x,y,z
 
 use intermediate_variables,only: e3 => tmpZ
                            
 implicit none
 real(kind=8)                                :: F(*)
 integer                                     :: iflag 
 integer,save                                :: istart,jstart,kstart
 integer,save                                :: nx,ny,nz
 logical,save                                :: first_entry=.TRUE.
 
 if( first_entry ) then
  istart = global_x_indices(START,ZBLOCK,myid)
  jstart = global_y_indices(START,ZBLOCK,myid)
  kstart = global_z_indices(START,ZBLOCK,myid)
  nx = array_size(JDIM,ZBLOCK,myid)   ! x is dim 2 in ZBLOCK
  ny = array_size(KDIM,ZBLOCK,myid)   ! y is dim 3 in ZBLOCK
  nz = array_size(IDIM,ZBLOCK,myid)   ! z is dim 1 in ZBLOCK
  first_entry=.FALSE.
 endif

 !!========================================================
 !! Compute coefficients of the 3D tensor basis functions
 !! using deBoor routine from nist/gams/cmlib
 !!========================================================  
 
  
 !-------------------------------
 !  3d  nx,ny > 1
 !-------------------------------
 if( nx > 1 .and. ny > 1 ) then
 
  iflag = 1   ! knots precomputed
  call db3ink(z(kstart),nz,                   &   ! ZBLOCK format
              x(istart),nx,                   &
              y(jstart),ny,                   &
              F,nz,nx,                        &
              order(3),order(1),order(2),     &
              tz(1),tx(1),ty(1),              &
              e3,wrk,iflag)
              
   if(iflag > 1 ) then
     write(0,*) ' myid ',myid
    write(0,*) 'error condition detected by db3ink:',iflag
    write(0,*) 'nx,ny,nz,order ',nx,ny,nz,order(1),order(2),order(3)
    write(0,*) ' x '
    write(0,*) x
    write(0,*) ' y '
    write(0,*) y
    write(0,*) ' z '
    write(0,*) z
    write(0,*) ' tx '
    write(0,*) tx
    write(0,*) ' ty '
    write(0,*) ty
    write(0,*) ' tz '
    write(0,*) tz
    stop
   endif
  
  !-------------------------------
  !  2d  ny=1
  !-------------------------------
  elseif( nx > 1 .and. ny == 1 ) then
    
   iflag = 1   ! knots precomputed
   call db2ink(z(kstart),nz,               &
               x(istart),nx,               &
               F,nz,order(3),order(1),     &
               tz(1),tx(1),e3,wrk,iflag)
               
   if(iflag > 1 ) then
     write(0,*) ' myid ',myid
    write(0,*) 'error condition detected by db2ink:',iflag
    write(0,*) 'nx,nz,order ',nx,nz,order(1),order(3)
    write(0,*) ' x '
    write(0,*) x
    write(0,*) ' z '
    write(0,*) z
    write(0,*) ' tx '
    write(0,*) tx
    write(0,*) ' tz '
    write(0,*) tz
    stop
   endif
  
  !-------------------------------
  !  2d  nx=1
  !-------------------------------
  elseif( nx == 1 .and. ny > 1 ) then
  
   iflag = 1   ! knots precomputed
   call db2ink(z(jstart),nz,               &
               y(kstart),ny,               &
               F,nz,order(3),order(2),     &
               tz(1),ty(1),e3,wrk,iflag)
               
   if(iflag > 1 ) then
    write(0,*) ' myid ',myid
    write(0,*) 'error condition detected by db2ink:',iflag
    write(0,*) 'ny,nz,order ',ny,nz,order(2),order(3)
    write(0,*) ' y '
    write(0,*) y
    write(0,*) ' z '
    write(0,*) z
    write(0,*) ' ty '
    write(0,*) ty
    write(0,*) ' tz '
    write(0,*) tz
    stop
   endif
  
  endif
end subroutine Spline_Interp_ZB




subroutine Spline_Eval_XB(g)
use decomposition_params
use mpi_params
use particles,                  only: my_labels,              &
                                      npts => nparticles,     &
                                      positions
                                      
use interpolation,               order=>interp_order_XBLOCK,  &
                                 tx => tx_XBLOCK,             &
                                 ty => ty_XBLOCK,             &
                                 tz => tz_XBLOCK,             &
                                 wrk => wrk_XBLOCK
                                 
use intermediate_variables,   only: e3 => tmpX

implicit none
integer                          :: m(3)=0   ! eval function only
integer                          :: i
real(kind=8)                     :: g(npts)
integer,save                     :: nx,ny,nz
logical,save                     :: first_entry=.TRUE.
real(kind=8),external            :: db2val,db3val
 
 if( first_entry ) then
  nx = array_size(IDIM,XBLOCK,myid)   ! x is dim 1 in XBLOCK
  ny = array_size(JDIM,XBLOCK,myid)   ! y is dim 2 in XBLOCK
  nz = array_size(KDIM,XBLOCK,myid)   ! z is dim 3 in XBLOCK
  first_entry=.FALSE.
 endif

!!===========================================================
!! Loop through each position & interpolate the function
!! note, positions detected as outside the appropriate 
!!         range will trigger a fatal error
!!===========================================================
 
  if( nx > 1 .and. ny > 1) then
   do i=1,npts
    
    if( my_labels(i)==1 ) then
     g(i) =  db3val(positions(i,1),                 &
                    positions(i,2),                 &
                    positions(i,3),                 &
                    m(1),m(2),m(3),                 &
                    tx,ty,tz,                       &
                    nx,ny,nz,                       &
                    order(1),order(2),order(3),     &
                    e3,wrk)
    else
     g(i)=0.d0
    endif
   enddo
   
   
   
  elseif( nx > 1 .and. ny==1 ) then
   do i=1,npts

    if( my_labels(i)==1 ) then
     g(i) =  db2val(positions(i,1),                 &
                    positions(i,3),                 &
                    m(1),m(3),                      &
                    tx,tz,                          &
                    nx,nz,                          &
                    order(1),order(3),              &
                    e3,wrk)                   
    else
     g(i)=0.d0
    endif
   enddo
   
   elseif( nx == 1 .and. ny>1 ) then
   do i=1,npts

    if( my_labels(i)==1 ) then
     g(i) =  db2val(positions(i,2),                 &
                    positions(i,3),                 &
                    m(2),m(3),                      &
                    ty,tz,                          &
                    ny,nz,                          &
                    order(2),order(3),              &
                    e3,wrk)                   
    else
     g(i)=0.d0
    endif
   enddo
   
  endif

end subroutine Spline_Eval_XB






subroutine Spline_Eval_YB(g)
use decomposition_params
use mpi_params
use particles,                  only: my_labels,              &
                                      npts => nparticles,     &
                                      positions
                                      
use interpolation,               order=>interp_order_YBLOCK,  &
                                 tx => tx_YBLOCK,             &
                                 ty => ty_YBLOCK,             &
                                 tz => tz_YBLOCK,             &
                                 wrk => wrk_YBLOCK
                                 
use intermediate_variables,   only: e3 => tmpY

implicit none
integer                          :: m(3)=0   ! eval function only
integer                          :: i
real(kind=8)                     :: g(npts)
integer,save                     :: nx,ny,nz
logical,save                     :: first_entry=.TRUE.
real(kind=8),external            :: db2val,db3val
 
 if( first_entry ) then
  nx = array_size(JDIM,YBLOCK,myid)   ! x is dim 2 in YBLOCK
  ny = array_size(IDIM,YBLOCK,myid)   ! y is dim 1 in YBLOCK
  nz = array_size(KDIM,YBLOCK,myid)   ! z is dim 3 in YBLOCK
  first_entry=.FALSE.
 endif

!!===========================================================
!! Loop through each position & interpolate the function
!! note, positions detected as outside the appropriate 
!!         range will trigger a fatal error
!!===========================================================
 
  if( nx > 1 .and. ny > 1) then
   do i=1,npts
    
    if( my_labels(i)==1 ) then
     g(i) =  db3val(positions(i,2),                 &
                    positions(i,1),                 &
                    positions(i,3),                 &
                    m(2),m(1),m(3),                 &
                    ty,tx,tz,                       &
                    ny,nx,nz,                       &
                    order(2),order(1),order(3),     &
                    e3,wrk)
    else
     g(i)=0.d0
    endif
   enddo
   
   
   
  elseif( nx > 1 .and. ny==1 ) then
   do i=1,npts

    if( my_labels(i)==1 ) then
     g(i) =  db2val(positions(i,1),                 &
                    positions(i,3),                 &
                    m(1),m(3),                      &
                    tx,tz,                          &
                    nx,nz,                          &
                    order(1),order(3),              &
                    e3,wrk)                   
    else
     g(i)=0.d0
    endif
   enddo
   
   elseif( nx == 1 .and. ny>1 ) then
   do i=1,npts

    if( my_labels(i)==1 ) then
     g(i) =  db2val(positions(i,2),                 &
                    positions(i,3),                 &
                    m(2),m(3),                      &
                    ty,tz,                          &
                    ny,nz,                          &
                    order(2),order(3),              &
                    e3,wrk)                   
    else
     g(i)=0.d0
    endif
   enddo
   
  endif

end subroutine Spline_Eval_YB



subroutine Spline_Eval_ZB(g)
use decomposition_params
use mpi_params
use particles,                  only: my_labels,              &
                                      npts => nparticles,     &
                                      positions
                                      
use interpolation,               order=>interp_order_ZBLOCK,  &
                                 tx => tx_ZBLOCK,             &
                                 ty => ty_ZBLOCK,             &
                                 tz => tz_ZBLOCK,             &
                                 wrk => wrk_ZBLOCK
                                 
use intermediate_variables,   only: e3 => tmpZ

implicit none
integer                          :: m(3)=0   ! eval function only
integer                          :: i
real(kind=8)                     :: g(npts)
integer,save                     :: nx,ny,nz
logical,save                     :: first_entry=.TRUE.
real(kind=8),external            :: db2val,db3val
 
 if( first_entry ) then
  nx = array_size(JDIM,ZBLOCK,myid)   ! x is dim 2 in ZBLOCK
  ny = array_size(KDIM,ZBLOCK,myid)   ! y is dim 3 in ZBLOCK
  nz = array_size(IDIM,ZBLOCK,myid)   ! z is dim 1 in ZBLOCK
  first_entry=.FALSE.
 endif

!!===========================================================
!! Loop through each position & interpolate the function
!! note, positions detected as outside the appropriate 
!!         range will trigger a fatal error
!!===========================================================
 
  if( nx > 1 .and. ny > 1) then
   do i=1,npts
    
    if( my_labels(i)==1 ) then
     g(i) =  db3val(positions(i,3),                 &
                    positions(i,1),                 &
                    positions(i,2),                 &
                    m(3),m(1),m(2),                 &
                    tz,tx,ty,                       &
                    nz,nx,ny,                       &
                    order(3),order(1),order(2),     &
                    e3,wrk)
    else
     g(i)=0.d0
    endif
   enddo
   
   
   
  elseif( nx > 1 .and. ny==1 ) then
   do i=1,npts

    if( my_labels(i)==1 ) then
     g(i) =  db2val(positions(i,3),                 &
                    positions(i,1),                 &
                    m(3),m(1),                      &
                    tz,tx,                          &
                    nz,nx,                          &
                    order(3),order(1),              &
                    e3,wrk)                   
    else
     g(i)=0.d0
    endif
   enddo
   
   elseif( nx == 1 .and. ny>1 ) then
   do i=1,npts

    if( my_labels(i)==1 ) then
     g(i) =  db2val(positions(i,3),                 &
                    positions(i,2),                 &
                    m(3),m(2),                      &
                    tz,ty,                          &
                    nz,ny,                          &
                    order(3),order(2),              &
                    e3,wrk)                   
    else
     g(i)=0.d0
    endif
   enddo
   
  endif

end subroutine Spline_Eval_ZB

subroutine DBKNOT_EXTRAP_BOTH(X,N,K,T)
!  --------------------------------------------------------------------
!  DBKNOT CHOOSES A KNOT SEQUENCE FOR INTERPOLATION OF ORDER K AT THE
!  DATA POINTS X(I), I=1,..,N.  THE N+K KNOTS ARE PLACED IN THE ARRAY
!  T.  K KNOTS ARE PLACED AT EACH ENDPOINT AND NOT-A-KNOT END
!  CONDITIONS ARE USED.  THE REMAINING KNOTS ARE PLACED AT DATA POINTS
!  IF N IS EVEN AND BETWEEN DATA POINTS IF N IS ODD.  THE RIGHTMOST
!  KNOT IS SHIFTED SLIGHTLY TO THE RIGHT TO INSURE PROPER INTERPOLATION
!  AT X(N) (SEE PAGE 350 OF THE REFERENCE).
!  DOUBLE PRECISION VERSION OF BKNOT.
!
!   KW  I pushed the ENDVALS out 1/2 gridpoint so that DBVALU
!       could be used for slight extrapolation as well
!  --------------------------------------------------------------------
!
!  ------------
!  DECLARATIONS
!  ------------
!
!  PARAMETERS
!
      integer       :: N,K
      real(kind=8)  :: X(N), T(*)
!
      INTEGER       :: I, J, IPJ, NPJ, IP1
      real(kind=8)  :: RNOT,LNOT

!  ----------------------------
!  PUT K KNOTS AT EACH ENDPOINT
!  ----------------------------
!     (SHIFT BOTH ENPOINTS 1/2 GRDPT)
      LNOT = X(1) - 0.750D0*( X(2)-X(1) )       !.5
      RNOT = X(N) + 0.750D0*( X(N)-X(N-1) )     !.5
      DO 110 J=1,K
         T(J) = LNOT
         NPJ = N + J
         T(NPJ) = RNOT
  110 CONTINUE
!  --------------------------
!  DISTRIBUTE REMAINING KNOTS
!  --------------------------
      IF (MOD(K,2) .EQ. 1)  GO TO 150
!
!     CASE OF EVEN K --  KNOTS AT DATA POINTS
!
      I = (K/2) - K
      JSTRT = K+1
      DO 120 J=JSTRT,N
         IPJ = I + J
         T(J) = X(IPJ)
  120 CONTINUE
      GO TO 200
!
!     CASE OF ODD K --  KNOTS BETWEEN DATA POINTS
!
  150 CONTINUE
      I = (K-1)/2 - K
      IP1 = I + 1
      JSTRT = K + 1
      DO 160 J=JSTRT,N
         IPJ = I + J
         T(J) = 0.50D0*( X(IPJ) + X(IPJ+1) )
  160 CONTINUE
  200 CONTINUE
!
      RETURN
END subroutine DBKNOT_EXTRAP_BOTH

 
subroutine DBKNOT_EXTRAP_LEFT(X,N,K,T)
!  --------------------------------------------------------------------
!  DBKNOT CHOOSES A KNOT SEQUENCE FOR INTERPOLATION OF ORDER K AT THE
!  DATA POINTS X(I), I=1,..,N.  THE N+K KNOTS ARE PLACED IN THE ARRAY
!  T.  K KNOTS ARE PLACED AT EACH ENDPOINT AND NOT-A-KNOT END
!  CONDITIONS ARE USED.  THE REMAINING KNOTS ARE PLACED AT DATA POINTS
!  IF N IS EVEN AND BETWEEN DATA POINTS IF N IS ODD.  THE RIGHTMOST
!  KNOT IS SHIFTED SLIGHTLY TO THE RIGHT TO INSURE PROPER INTERPOLATION
!  AT X(N) (SEE PAGE 350 OF THE REFERENCE).
!  DOUBLE PRECISION VERSION OF BKNOT.
!
!   KW  I pushed the ENDVALS out 1/2 gridpoint so that DBVALU
!       could be used for slight extrapolation as well
!  --------------------------------------------------------------------
!
!  ------------
!  DECLARATIONS
!  ------------
!
!  PARAMETERS
!
      integer       :: N,K
      real(kind=8)  :: X(N), T(*)
!
      INTEGER       :: I, J, IPJ, NPJ, IP1
      real(kind=8)  :: RNOT,LNOT

!  ----------------------------
!  PUT K KNOTS AT EACH ENDPOINT
!  ----------------------------
!     (SHIFT LEFT ENPOINT 1/2 GRDPT)
      LNOT = X(1) - 0.50D0*( X(2)-X(1) )
      RNOT = X(N) + 0.10D0*( X(N)-X(N-1) )  ! default config already does this slight shift, keep
      DO 110 J=1,K
         T(J) = LNOT
         NPJ = N + J
         T(NPJ) = RNOT
  110 CONTINUE
!  --------------------------
!  DISTRIBUTE REMAINING KNOTS
!  --------------------------
      IF (MOD(K,2) .EQ. 1)  GO TO 150
!
!     CASE OF EVEN K --  KNOTS AT DATA POINTS
!
      I = (K/2) - K
      JSTRT = K+1
      DO 120 J=JSTRT,N
         IPJ = I + J
         T(J) = X(IPJ)
  120 CONTINUE
      GO TO 200
!
!     CASE OF ODD K --  KNOTS BETWEEN DATA POINTS
!
  150 CONTINUE
      I = (K-1)/2 - K
      IP1 = I + 1
      JSTRT = K + 1
      DO 160 J=JSTRT,N
         IPJ = I + J
         T(J) = 0.50D0*( X(IPJ) + X(IPJ+1) )
  160 CONTINUE
  200 CONTINUE
!
      RETURN
END subroutine DBKNOT_EXTRAP_LEFT

subroutine DBKNOT_EXTRAP_RIGHT(X,N,K,T)
!  --------------------------------------------------------------------
!  DBKNOT CHOOSES A KNOT SEQUENCE FOR INTERPOLATION OF ORDER K AT THE
!  DATA POINTS X(I), I=1,..,N.  THE N+K KNOTS ARE PLACED IN THE ARRAY
!  T.  K KNOTS ARE PLACED AT EACH ENDPOINT AND NOT-A-KNOT END
!  CONDITIONS ARE USED.  THE REMAINING KNOTS ARE PLACED AT DATA POINTS
!  IF N IS EVEN AND BETWEEN DATA POINTS IF N IS ODD.  THE RIGHTMOST
!  KNOT IS SHIFTED SLIGHTLY TO THE RIGHT TO INSURE PROPER INTERPOLATION
!  AT X(N) (SEE PAGE 350 OF THE REFERENCE).
!  DOUBLE PRECISION VERSION OF BKNOT.
!
!   KW  I pushed the ENDVALS out 1/2 gridpoint so that DBVALU
!       could be used for slight extrapolation as well
!  --------------------------------------------------------------------
!
!  ------------
!  DECLARATIONS
!  ------------
!
!  PARAMETERS
!
      integer       :: N,K
      real(kind=8)  :: X(N), T(*)
!
      INTEGER       :: I, J, IPJ, NPJ, IP1
      real(kind=8)  :: RNOT,LNOT

!  ----------------------------
!  PUT K KNOTS AT EACH ENDPOINT
!  ----------------------------
!     (SHIFT RIGHT ENPOINT 1/2 GRDPT)
      LNOT = X(1)
      RNOT = X(N) + 0.50D0*( X(N)-X(N-1) )
      DO 110 J=1,K
         T(J) = LNOT
         NPJ = N + J
         T(NPJ) = RNOT
  110 CONTINUE
!  --------------------------
!  DISTRIBUTE REMAINING KNOTS
!  --------------------------
      IF (MOD(K,2) .EQ. 1)  GO TO 150
!
!     CASE OF EVEN K --  KNOTS AT DATA POINTS
!
      I = (K/2) - K
      JSTRT = K+1
      DO 120 J=JSTRT,N
         IPJ = I + J
         T(J) = X(IPJ)
  120 CONTINUE
      GO TO 200
!
!     CASE OF ODD K --  KNOTS BETWEEN DATA POINTS
!
  150 CONTINUE
      I = (K-1)/2 - K
      IP1 = I + 1
      JSTRT = K + 1
      DO 160 J=JSTRT,N
         IPJ = I + J
         T(J) = 0.50D0*( X(IPJ) + X(IPJ+1) )
  160 CONTINUE
  200 CONTINUE
!
      RETURN
END subroutine DBKNOT_EXTRAP_RIGHT






subroutine Spline_Eval_YB_ptwise(xval,yval,zval,ans)
use decomposition_params                                      
use interpolation,               order=>interp_order_YBLOCK,  &
                                 tx => tx_YBLOCK,             &
                                 ty => ty_YBLOCK,             &
                                 tz => tz_YBLOCK,             &
                                 wrk => wrk_YBLOCK
                                 
use intermediate_variables,   only: e3 => tmpY
use mpi_params,               only: myid

implicit none
real(kind=8)                     :: xval,yval,zval,ans
integer                          :: m(3)=0   ! eval function only
integer,save                     :: nx,ny,nz,npts
logical,save                     :: first_entry=.TRUE.
real(kind=8),external            :: db2val,db3val

 
 if( first_entry ) then
  nx = array_size(JDIM,YBLOCK,myid)   ! x is dim 2 in YBLOCK
  ny = array_size(IDIM,YBLOCK,myid)   ! y is dim 1 in YBLOCK
  nz = array_size(KDIM,YBLOCK,myid)   ! z is dim 3 in YBLOCK
  npts = 1                       
    
  first_entry=.FALSE.
 endif


  if( nx > 1 .and. ny > 1) then
    
      ans     =  db3val(yval,                           &
                        xval,                           &
                        zval,                           &
                        m(2),m(1),m(3),                 &
                        ty,tx,tz,                       &
                        ny,nx,nz,                       &
                        order(2),order(1),order(3),     &
                        e3,wrk)
   
   
   
  elseif( nx > 1 .and. ny==1 ) then

     ans      = db2val(xval,                           &
                       zval,                           &
                       m(1),m(3),                      &
                       tx,tz,                          &
                       nx,nz,                          &
                       order(1),order(3),              &
                       e3,wrk) 
   
  elseif( nx == 1 .and. ny>1 ) then

      ans     = db2val(yval,                           &
                       zval,                           &
                       m(2),m(3),                      &
                       ty,tz,                          &
                       ny,nz,                          &
                       order(2),order(3),              &
                       e3,wrk) 
   
  endif

return
end subroutine Spline_Eval_YB_ptwise


subroutine Spline_Eval_YB_vec(xval,yval,zval,ans,N)
use decomposition_params                                      
use interpolation,               order=>interp_order_YBLOCK,  &
                                 tx => tx_YBLOCK,             &
                                 ty => ty_YBLOCK,             &
                                 tz => tz_YBLOCK,             &
                                 wrk => wrk_YBLOCK
                                 
use intermediate_variables,   only: e3 => tmpY
use mpi_params,               only: myid

implicit none
integer                          :: i,N
real(kind=8)                     :: xval(N),yval(N),zval(N),ans(N)
integer                          :: m(3)=0   ! eval function only
integer,save                     :: nx,ny,nz
logical,save                     :: first_entry=.TRUE.
real(kind=8),external            :: db2val,db3val

 
 if( first_entry ) then
  nx = array_size(JDIM,YBLOCK,myid)   ! x is dim 2 in YBLOCK
  ny = array_size(IDIM,YBLOCK,myid)   ! y is dim 1 in YBLOCK
  nz = array_size(KDIM,YBLOCK,myid)   ! z is dim 3 in YBLOCK
    
  first_entry=.FALSE.
 endif


  if( nx > 1 .and. ny > 1) then
    
     do i=1,N
      ans(i)     =  db3val(yval(i),                           &
                           xval(i),                           &
                           zval(i),                           &
                           m(2),m(1),m(3),                    &
                           ty,tx,tz,                          &
                           ny,nx,nz,                          &
                           order(2),order(1),order(3),        &
                           e3,wrk)
     enddo
   
   
  elseif( nx > 1 .and. ny==1 ) then

    do i=1,N
     ans(i)    = db2val(xval(i),                        &
                        zval(i),                        &
                        m(1),m(3),                      &
                        tx,tz,                          &
                        nx,nz,                          &
                        order(1),order(3),              &
                        e3,wrk) 
    enddo
   
  elseif( nx == 1 .and. ny>1 ) then

     do i=1,N
      ans(i)  = db2val(yval(i),                        &
                       zval(i),                        &
                       m(2),m(3),                      &
                       ty,tz,                          &
                       ny,nz,                          &
                       order(2),order(3),              &
                       e3,wrk)
     enddo
   
  endif

return
end subroutine Spline_Eval_YB_vec







subroutine Spline_Eval_ZB_ptwise(xval,yval,zval,ans)
use decomposition_params
use mpi_params
use interpolation,               order=>interp_order_ZBLOCK,  &
                                 tx => tx_ZBLOCK,             &
                                 ty => ty_ZBLOCK,             &
                                 tz => tz_ZBLOCK,             &
                                 wrk => wrk_ZBLOCK
                                 
use intermediate_variables,   only: e3 => tmpZ

implicit none
real(kind=8)                     :: xval,yval,zval,ans
integer                          :: m(3)=0   ! eval function only
integer,save                     :: nx,ny,nz,npts
logical,save                     :: first_entry=.TRUE.
real(kind=8),external            :: db2val,db3val
 
 if( first_entry ) then
  nx = array_size(JDIM,ZBLOCK,myid)   ! x is dim 2 in ZBLOCK
  ny = array_size(KDIM,ZBLOCK,myid)   ! y is dim 3 in ZBLOCK
  nz = array_size(IDIM,ZBLOCK,myid)   ! z is dim 1 in ZBLOCK
  npts=1
  first_entry=.FALSE.
 endif

!!===========================================================
!! Loop through each position & interpolate the function
!! note, positions detected as outside the appropriate 
!!         range will trigger a fatal error
!!===========================================================
 
  if( nx > 1 .and. ny > 1) then
    
     ans  =  db3val(zval,                           &
                    xval,                           &
                    yval,                           &
                    m(3),m(1),m(2),                 &
                    tz,tx,ty,                       &
                    nz,nx,ny,                       &
                    order(3),order(1),order(2),     &
                    e3,wrk)
   
   
   
  elseif( nx > 1 .and. ny==1 ) then

     ans =   db2val(zval,                           &
                    xval,                           &
                    m(3),m(1),                      &
                    tz,tx,                          &
                    nz,nx,                          &
                    order(3),order(1),              &
                    e3,wrk)                   
   
   elseif( nx == 1 .and. ny>1 ) then

     ans =   db2val(zval,                           &
                    yval,                           &
                    m(3),m(2),                      &
                    tz,ty,                          &
                    nz,ny,                          &
                    order(3),order(2),              &
                    e3,wrk)                   
   
  endif

end subroutine Spline_Eval_ZB_ptwise
      
