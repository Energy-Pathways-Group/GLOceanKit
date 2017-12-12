subroutine SetupDomain
 use etc
 use mpi_params
 use dimensional_scales
 use independent_variables
 use methods_params
 use decomposition_params
 implicit none
 
 real(kind=8)               :: pi
 integer                    :: i,k
 logical                    :: closed_interval(3)=.TRUE.

 if(myid==0) then
  write(0,*) '................'
  write(0,*) '................      hello world from SetupDomain'
  open(1,file=logfile,position='append') 
  write(1,*) '  '
  write(1,*) '  '
  write(1,*) ' =========================================================== '
  write(1,*) ' =========================================================== '
  write(1,*) '                  SetupDomain Report:'
  write(1,*) ' =========================================================== '
  write(1,*) ' =========================================================== '
  write(1,*) '  '
 endif
  
 allocate ( x(nx),x_s(nx),x_ss(nx),s_x(nx) )  
 allocate ( y(ny),y_s(ny),y_ss(ny),s_y(ny) )
 allocate ( z(nz),z_s(nz),z_ss(nz),s_z(nz) )
 
 !! dimension with 1 xtra location, needed for
 !! interpolation in open intervals, otherwise
 !! extra value not actually used
 allocate ( s_of_x(nx+1) )
 allocate ( s_of_y(ny+1) )
 allocate ( s_of_z(nz+1) )
 
 pi=4.d0*datan(1.d0)
 
 !! periodic coordinates have open intervals
 if( xdim_periodic ) closed_interval(1)=.FALSE.
 if( ydim_periodic ) closed_interval(2)=.FALSE.
 if( zdim_periodic ) closed_interval(3)=.FALSE.
 
 !!==================================================
 !! computational coordinate s 
 !! [0,1]  for compact, sin, cos
 !! [0,1)  for periodic dimensions
 !! [0,pi] for cheby
 !!================================================== 
 if( closed_interval(IDIM) ) then  
  if( trim(deriv_type(1,IDIM,1)) == 'cheby' ) then
   ds(IDIM) = pi/(dfloat(nx)-1.d0)
  else
   ds(IDIM) = 1.d0/(dfloat(nx)-1.d0)
  endif
 else
  ds(IDIM) = 1.d0/dfloat(nx)
 endif 
 do i=1,nx+1
  s_of_x(i) = (dfloat(i)-1.d0)*ds(IDIM)
 enddo
 
 if( closed_interval(JDIM) ) then
  if( trim(deriv_type(1,JDIM,1)) == 'cheby' ) then
   ds(JDIM) = pi/(dfloat(ny)-1.d0)
  else
   ds(JDIM) = 1.d0/(dfloat(ny)-1.d0)
  endif
 else
  ds(JDIM) = 1.d0/dfloat(ny)
 endif
 do i=1,ny+1
  s_of_y(i) = (dfloat(i)-1.d0)*ds(JDIM)
 enddo
 
 if( closed_interval(KDIM) ) then
  if( trim(deriv_type(1,KDIM,1)) == 'cheby' ) then
   ds(KDIM) = pi/(dfloat(nz)-1.d0)
  else
   ds(KDIM) = 1.d0/(dfloat(nz)-1.d0)
  endif
 else
  ds(KDIM) = 1.d0/dfloat(nz)
 endif 
 do i=1,nz+1
  s_of_z(i) = (dfloat(i)-1.d0)*ds(KDIM)
 enddo

 !!=============================================================
 !! call user routine to get the mapping parameters 
 !!=============================================================
 
 ! no longer ask user... set for NO stretching...
 !call coordinate_map_params( cmap_params )
 cmap_params(:,:)=1.d0
 cmap_params(2,IDIM)=Lx
 cmap_params(2,JDIM)=Ly
 cmap_params(2,KDIM)=Lz
 stretch(:)=.FALSE.      !! always for this version...

 
 !!=============================================================
 ! for cheby dimensions, set indicator in case user forgot
 !!=============================================================
 do i=1,3
  if( trim(deriv_type(1,i,1)) == 'cheby' ) then
   write(0,*) 'CHEBY OBSOLETE...setup_domain.f90'
   STOP
   cmap_params(1,i)=-999
   stretch(i)=.FALSE.    ! chain rule built into cheby recursion
  endif
 enddo
 
 
 !!=============================================================================
 !! define x coordinate mapping 
 !!=============================================================================
 if( cmap_params(1,IDIM) == -999 .and. deriv_type(1,IDIM,1) .ne. 'cheby' ) &
  stop 'gauss lobatto pts specified for non cheby differentiation in x'
 call evaluate_mapping(nx,ds(IDIM),Lx,x,x_s,x_ss,cmap_params(1,IDIM))
 
  
 !!=============================================================================
 !! define y coordinate mapping 
 !!=============================================================================
 if( cmap_params(1,JDIM) == -999 .and. deriv_type(1,JDIM,1) .ne. 'cheby' ) &
  stop 'gauss lobatto pts specified for non cheby differentiation in y'
 call evaluate_mapping(ny,ds(JDIM),Ly,y,y_s,y_ss,cmap_params(1,JDIM))
 
 
 !!=============================================================================
 !! define z coordinate mapping 
 !!=============================================================================
 if( cmap_params(1,KDIM) == -999 .and. deriv_type(1,KDIM,1) .ne. 'cheby' ) &
  stop 'gauss lobatto pts specified for non cheby differentiation in z'
 call evaluate_mapping(nz,ds(KDIM),Lz,z,z_s,z_ss,cmap_params(1,KDIM))
 
 
 
 !!================================================================
 !! Scale the results using "length_scale" [m], s is dimensionless
 !!================================================================
 x(:)=x(:)/length_scale
 x_s(:)=x_s(:)/length_scale
 x_ss(:)=x_ss(:)/length_scale 
 do k=1,nx
  s_x(k) = 1.d0/x_s(k)
 enddo
 
 y(:)=y(:)/length_scale
 y_s(:)=y_s(:)/length_scale
 y_ss(:)=y_ss(:)/length_scale
 do k=1,ny
  s_y(k) = 1.d0/y_s(k)
 enddo
 
 z(:)=z(:)/length_scale
 z_s(:)=z_s(:)/length_scale
 z_ss(:)=z_ss(:)/length_scale
 do k=1,nz
  s_z(k) = 1.d0/z_s(k)
 enddo
  
 if(myid==0) then
  write(1,*) ' -----> SetupDomain routine exiting normally  <---------- '
  close(1)
 endif

	 
 return
end subroutine SetupDomain


 subroutine evaluate_mapping(n,ds,Lx,x,dxds,d2xds2,cmap)       
  implicit none
  integer            :: n,i
  real(kind=8)       :: ds
  real(kind=8)       :: Lx
  real(kind=8)       :: x(*)
  real(kind=8)       :: dxds(*)
  real(kind=8)       :: d2xds2(*)
  real(kind=8)       :: cmap(5)
  real(kind=8)       :: b
  real(kind=8)       :: p,beta,s0,s
  real(kind=8)       :: alpha,gamma
  real(kind=8)       :: num,denom
  real(kind=8)       :: pi,theta

  
  if(n==1) then 
   x(1)=0.d0
   dxds(1)=Lx  !! I divide by dxds to get dsdx
   d2xds2(1)=0.d0
   return
  endif
 
  b = Lx 
  p=cmap(1)
  beta=cmap(2)
  s0=cmap(3)

 if( p .ne. -999 ) then        !! s in  [0,1], x_s generally not constant
  !! some derived quantities    
   num=(2*p+1.d0)*(b-beta);
   denom=b*( ((1.d0-s0)/b)**(2*p+1.d0) + (s0/b)**(2*p+1.d0) )
   alpha = num/denom;
   gamma = (alpha*b/(2*p+1.d0) )*(s0/b)**(2*p+1.d0)
   cmap(4) = alpha
   cmap(5) = gamma
    
  !! evaluate the mapping at gridpoints
   do i=1,n
    s=(dfloat(i)-1.d0)*ds
    x(i)=(alpha*b/(2.d0*p+1.d0))*((s-s0)/b)**(2.d0*p+1.d0) + beta*s + gamma
    dxds(i)= alpha*((s-s0)/b)**(2*p) + beta
    d2xds2(i)=(2.d0*p*alpha/b)*((s-s0)/b)**(2.d0*p-1.d0)
   enddo
  elseif( p == -999 ) then  !! cheby method: s=theta in [0,pi], x in [0,Lx]
   pi=4.d0*datan(1.d0)
   do i=1,n
    theta=pi*(dfloat(i)-1.d0)/(dfloat(n)-1.d0)
    x(i)=-(Lx/2.d0)*( dcos(theta) -1.d0 )
    dxds(i)=(Lx/2.d0)*dsin(theta)
    d2xds2(i)=(Lx/2.d0)*dcos(theta)
   enddo
  endif
                 
 end subroutine evaluate_mapping


