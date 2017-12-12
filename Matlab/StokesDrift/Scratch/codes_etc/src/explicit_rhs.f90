subroutine explicit_rhs_1
 use decomposition_params
 use mpi_params,             only: myid 
 use etc,                    only: MM0,N
 use pde_params,             only: Rot,f_plane_key 
 use methods_params,         only: deriv_type,do_nonlinear
 use independent_variables,  only: nx,ny,nz,Lx,Ly,Lz,x,y,z
 use dimensional_scales,     only: length_scale,y_pivot
 use dependent_variables,    only: u,v,w      
 use intermediate_variables, only: tmpY,          &
                                   explicit_rhs,  &
                                   implicit_rhs
                                    
 implicit none  !
 real(kind=8)                            :: y0
 integer                                 :: i,j,k,kg
 integer                                 :: order
 integer,parameter                       :: id=1
 integer,parameter                       :: sign=-1
 integer,parameter                       :: inc=1
 character(len=80),dimension(3)          :: method
 
 real(kind=8),dimension(:,:,:),pointer   :: rhs      ! -u dot grad u + fv 
                                                     ! (+ beta*(y-y0)*v)
                                                     ! (- f~ * w )
 real(kind=8),dimension(:,:,:),pointer   :: diff     ! generalized {D[u]}
 
 !-----------------------------------------------------------
 !-----  params for forcing f~,g -> 0 at z=0,Lz and
 !                          f,beta -> 0 at sidewalls
 !    not always needed...
 !-----------------------------------------------------------
  real(kind=8),save            :: dx,sigma_x,xl,xr
  real(kind=8),save            :: dy,sigma_y,yl,yr
  real(kind=8),save            :: dz,sigma_z,zb,zt,xx
  logical,save                 :: first_entry=.TRUE.
  
  if( first_entry ) then
   if( nx > 1 ) then
    dx = ( Lx/(nx-1.d0) )/length_scale
   else
    dx=1./length_scale
   endif
   if( ny > 1 ) then
    dy = ( Ly/(ny-1.d0) )/length_scale
   else
    dy = 1./length_scale
   endif
   dz = ( Lz/(nz-1.d0) )/length_scale
   sigma_x = 1.d0 * dx
   sigma_y = 1.d0 * dy
   sigma_z = 1.d0 * dz
   xl = 0.d0
   xr = Lx/length_scale
   yl = 0.d0
   yr = Ly/length_scale
   zb = 0.d0
   zt = Lz/length_scale
   first_entry=.FALSE.
  endif
 !-----------------------------------------------------------
 
 
 
 rhs    => explicit_rhs(:,:,:,id,MM0)
 diff   => implicit_rhs(:,:,:,id,N)

!! for 2d problems in yz plane w/ no rotation 
!!  ==> u is identically zero at all times
 if( nx==1 .and. Rot(1) == 0.d0 ) then  
  return
 endif
 
 !-----------------------------------------------------------------
 ! grad u 
 !        store results in tmpY(:,:,:,1-3)
 !-----------------------------------------------------------------
 order=1
 method(:)=deriv_type(id,:,order)
 call gradient(u,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),method)

 
 !-----------------------------------------------------------------
 ! rhs <==  -(u dot grad u)
 !-----------------------------------------------------------------
 if( do_nonlinear ) then
  call YBdotproduct(rhs,sign,u,v,w,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3))
 else
  rhs(:,:,:)=0.d0
 endif
  

 if( Rot(1) /= 0.d0 ) then  ! i.e.   if f0 is nonzero...
 
  !----------------------------------------------------------------
  ! for solid walls at x=0,Lx  fv produces acceleration through walls
  ! clip these accelerations (smooth approx delta fn) so that
  ! homogeneous normal pressure gradients can be used...
  !----------------------------------------------------------------
  if( deriv_type(1,1,1)=='sin' .and. deriv_type(2,1,1)=='cos' ) then  
   do i=1,nx
    xx = 1.d0 - exp( -((x(i)-xl)/sigma_x)**2 ) - exp( -((x(i)-xr)/sigma_x)**2 )
    rhs(:,i,:) = rhs(:,i,:) + Rot(1)*xx*v(:,i,:)
   enddo
  else
   rhs(:,:,:) = rhs(:,:,:) + Rot(1)*v(:,:,:)
  endif
   
 endif
 
 
 !-----------------------------------------------------------------
 ! horizontal component of coriolis acceleration...
 ! rhs <==  rhs - 'f~' * w  ; Rot(3) is "dimensionless f~"
 !   clip f~->0 az z-> rigid lids, 
 !   this sets bcs on w* and allows homogeneous bcs for pressure
 !-----------------------------------------------------------------
 if( trim(f_plane_key)=='non-traditional' ) then
 
  if( deriv_type(3,3,1)=='sin' .and. deriv_type(1,3,1)=='cos' ) then  ! w=sin(z),u=cos(z)
   do k=1,array_size(KDIM,YBLOCK,myid)   ! z is dim3 in YBLOCK format
    kg=global_z_indices(START,YBLOCK,myid) + k - 1
    xx = 1.d0 - exp( -((z(kg)-zb)/sigma_z)**2 ) - exp( -((z(kg)-zt)/sigma_z)**2 )
    rhs(:,:,k) = rhs(:,:,k) - Rot(3)*w(:,:,k)*xx
   enddo
  endif
    
 endif
 
 
 !-----------------------------------------------------------------
 ! beta efffect ...
 ! rhs <== rhs + 'beta' * (y-y0) * v   ! y is north here
 !                Rot(2) is "dimensionless beta"
 !-----------------------------------------------------------------
 if( Rot(2) /= 0 ) then
  
   y0 = y_pivot/length_scale
   
   if( deriv_type(1,1,1)=='sin' .and. deriv_type(2,1,1)=='cos' ) then  ! u=sin,v=cos in x...
   !----------------------------------------------------------------
   ! for solid walls at x=0,Lx  fv produces acceleration through walls
   ! clip these accelerations (smooth approx delta fn) so that
   ! homogeneous normal pressure gradients can be used...
   !----------------------------------------------------------------
    do j=1,ny
     do i=1,nx
      xx = 1.d0 - exp( -((x(i)-xl)/sigma_x)**2 ) - exp( -((x(i)-xr)/sigma_x)**2 )
      rhs(j,i,:) = rhs(j,i,:) + Rot(2)*xx*(y(j)-y0)*v(j,i,:)
     enddo
    enddo
   else
    do j=1,ny
     rhs(j,:,:) = rhs(j,:,:) + Rot(2)*(y(j)-y0)*v(j,:,:)
    enddo
   endif
    
 endif 
 
 !-----------------------------------------------------------------
 ! apply the generalized diffusion operator to u
 !  (use previously computed z derivative stored in tmpY(1,1,1,3))
 ! store the result in   diff  => implicit_rhs(:,:,:,id,N)
 ! NB
 !  diffusion uses tmpX(:,:,:,1-2),tmpY(:,:,:,1-2),tmpZ(:,:,:,1-2)
 !  internally for work space
 !-----------------------------------------------------------------
 call diffusion(u,tmpY(1,1,1,3),diff,id)
 

end subroutine explicit_rhs_1





subroutine explicit_rhs_2
 use decomposition_params
 use mpi_params,             only: myid 
 use etc,                    only: MM0,N
 use pde_params,             only: Rot 
 use methods_params,         only: deriv_type,do_nonlinear
 use independent_variables,  only: nx,ny,nz,Lx,Ly,Lz,x,y,z
 use dimensional_scales,     only: length_scale,y_pivot
 use dependent_variables,    only: u,v,w      
 use intermediate_variables, only: tmpY,          &
                                   explicit_rhs,  &
                                   implicit_rhs
                                  
 implicit none  
 real(kind=8)                            :: y0
 integer                                 :: i,j,k
 integer                                 :: order
 integer,parameter                       :: id=2
 integer,parameter                       :: sign=-1
 character(len=80),dimension(3)          :: method
 
 real(kind=8),dimension(:,:,:),pointer   :: rhs      ! -u dot grad v - fu 
                                                     ! (- beta*(y-y0)*u)
                                                     
 real(kind=8),dimension(:,:,:),pointer   :: diff     ! generalized {D[v]}
 
 !-----------------------------------------------------------
 !-----  params for forcing f,beta -> 0 at sidewalls etc
 !    not always needed...
 !-----------------------------------------------------------
  real(kind=8),save            :: dx,sigma_x,xl,xr
  real(kind=8),save            :: dy,sigma_y,yl,yr
  real(kind=8),save            :: dz,sigma_z,zb,zt,xx
  logical,save                 :: first_entry=.TRUE.
  
  if( first_entry ) then
   if( nx > 1 ) then
    dx = ( Lx/(nx-1.d0) )/length_scale
   else
    dx=1./length_scale
   endif
   if( ny > 1 ) then
    dy = ( Ly/(ny-1.d0) )/length_scale
   else
    dy = 1./length_scale
   endif
   dz = ( Lz/(nz-1.d0) )/length_scale
   sigma_x = 1.d0 * dx
   sigma_y = 1.d0 * dy
   sigma_z = 1.d0 * dz
   xl = 0.d0
   xr = Lx/length_scale
   yl = 0.d0
   yr = Ly/length_scale
   zb = 0.d0
   zt = Lz/length_scale
   first_entry=.FALSE.
  endif
 !-----------------------------------------------------------
 
 rhs    => explicit_rhs(:,:,:,id,MM0)
 diff   => implicit_rhs(:,:,:,id,N)
  
 
 !-----------------------------------------------------------------
 ! grad v 
 !        store results in tmpY(:,:,:,1-3)
 !-----------------------------------------------------------------
 order=1
 method(:)=deriv_type(id,:,order)
 call gradient(v,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),method)
 
 !-----------------------------------------------------------------
 ! rhs <==  -(u dot grad v)
 !-----------------------------------------------------------------
 if( do_nonlinear ) then
  call YBdotproduct(rhs,sign,u,v,w,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3))
 else
  rhs(:,:,:)=0.d0
 endif
 
 !-----------------------------------------------------------------
 ! rhs <==  rhs - 'f0' * u  ; Rot(1) is "dimensionless f0"
 !-----------------------------------------------------------------
 if( Rot(1) /= 0.d0 ) then
 
  if( deriv_type(2,2,1)=='sin' .and. deriv_type(1,2,1)=='cos' ) then  
  !----------------------------------------------------------------
  ! for solid walls at y=0,Ly  fu produces acceleration through walls
  ! clip these accelerations (smooth approx delta fn) so that
  ! homogeneous normal pressure gradients can be used...
  !----------------------------------------------------------------
   do j=1,ny
    xx = 1.d0 - exp( -((y(j)-yl)/sigma_y)**2 ) - exp( -((y(j)-yr)/sigma_y)**2 )
    rhs(j,:,:) = rhs(j,:,:) - Rot(1)*xx*u(j,:,:)
   enddo
  else
   do j=1,ny
    rhs(j,:,:) = rhs(j,:,:) - Rot(1)*u(j,:,:)
   enddo
  endif
    
 endif
 
 !-----------------------------------------------------------------
 ! beta efffect ...
 ! rhs <== rhs - 'beta' * (y-y0) * u   ! y is north here
 !                Rot(2) is "dimensionless beta"
 !-----------------------------------------------------------------
 if( Rot(2) /= 0 ) then
 
   y0 = y_pivot/length_scale
   if( deriv_type(2,2,1)=='sin' .and. deriv_type(1,2,1)=='cos' ) then  
    !----------------------------------------------------------------
    ! for solid walls at y=0,Ly  fu produces acceleration through walls
    ! clip these accelerations (smooth approx delta fn) so that
    ! homogeneous normal pressure gradients can be used...
    !----------------------------------------------------------------
    do j=1,ny
     xx = 1.d0 - exp( -((y(j)-yl)/sigma_y)**2 ) - exp( -((y(j)-yr)/sigma_y)**2 )
     rhs(j,:,:) = rhs(j,:,:) - Rot(2)*xx*(y(j)-y0)*u(j,:,:)
    enddo
   else
    do j=1,ny
     rhs(j,:,:) = rhs(j,:,:) - Rot(2)*(y(j)-y0)*u(j,:,:)
    enddo
   endif

 endif 

 !-----------------------------------------------------------------
 ! apply the generalized diffusion operator to v
 !  (use previously computed z derivative stored in tmpY(1,1,1,3))
 ! store the result in   diff  => implicit_rhs(:,:,:,id,N)
 ! NB
 !  diffusion uses tmpX(:,:,:,1-2),tmpY(:,:,:,1-2),tmpZ(:,:,:,1-2)
 !  internally for work space
 !-----------------------------------------------------------------
 call diffusion(v,tmpY(1,1,1,3),diff,id)

end subroutine explicit_rhs_2




subroutine explicit_rhs_3
 use decomposition_params
 use mpi_params,             only: myid 
 use etc,                    only: MM0,N
 use pde_params,             only: Ri,Rot,f_plane_key 
 use methods_params,         only: deriv_type,do_nonlinear
 use dependent_variables,    only: u,v,w   
 use independent_variables,  only: nx,ny,nz,Lx,Ly,Lz,x,y,z
 use intermediate_variables, only: pd,tmpY,       &
                                   explicit_rhs,  &
                                   implicit_rhs
 use dimensional_scales,     only: length_scale
 implicit none  
 integer                                 :: order,k,kg
 integer,parameter                       :: id=3
 integer,parameter                       :: sign=-1
 character(len=80),dimension(3)          :: method
 
 real(kind=8),dimension(:,:,:),pointer   :: rhs      ! -u dot grad w - Ri*pd 
                                                     ! (+ f~ * u )
 real(kind=8),dimension(:,:,:),pointer   :: diff     ! generalized {D[w]}
 
 !--------------------------------------------------------------------
 !-----  params for dealing with buoyancy, nt coriolis accelerations
 !-----  near rigid lids z=0,Lz when rho' is no flux
 !--------------------------------------------------------------------
  real(kind=8),save            :: dx,sigma_x,xl,xr
  real(kind=8),save            :: dy,sigma_y,yl,yr
  real(kind=8),save            :: dz,sigma_z,zb,zt,xx
  logical,save                 :: first_entry=.TRUE.
  
  if( first_entry ) then
   if( nx > 1 ) then
    dx = ( Lx/(nx-1.d0) )/length_scale
   else
    dx=1./length_scale
   endif
   if( ny > 1 ) then
    dy = ( Ly/(ny-1.d0) )/length_scale
   else
    dy = 1./length_scale
   endif
   dz = ( Lz/(nz-1.d0) )/length_scale
   sigma_x = 1.d0 * dx
   sigma_y = 1.d0 * dy
   sigma_z = 1.d0 * dz
   xl = 0.d0
   xr = Lx/length_scale
   yl = 0.d0
   yr = Ly/length_scale
   zb = 0.d0
   zt = Lz/length_scale
   first_entry=.FALSE.
  endif
 !-----------------------------------------------------------
 
 rhs    => explicit_rhs(:,:,:,id,MM0)
 diff   => implicit_rhs(:,:,:,id,N)
   
 
 !-----------------------------------------------------------------
 ! grad w 
 !        store results in tmpY(:,:,:,1-3)
 !-----------------------------------------------------------------
 order=1
 method(:)=deriv_type(id,:,order)
 call gradient(w,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),method)
 
 !-----------------------------------------------------------------
 ! rhs <==  -(u dot grad w)
 !-----------------------------------------------------------------
 if( do_nonlinear ) then
  call YBdotproduct(rhs,sign,u,v,w,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3))
 else
  rhs(:,:,:) = 0.d0
 endif
 
 !-----------------------------------------------------------------
 ! buoyancy term ...
 ! rhs <==  rhs - Ri * pd
 !   take g->0 as z-> rigid lids (smooth, approximate delta function) 
 !   this sets bcs on w* and allows homogeneous bcs for pressure
 !-----------------------------------------------------------------
 if( Ri /= 0.d0 ) then
 
  if( deriv_type(3,3,1)=='sin' .and. deriv_type(4,3,1)=='cos' ) then  ! w=sin(z), rho'=s1==cos(z)
   do k=1,array_size(KDIM,YBLOCK,myid)   ! z is dim3 in YBLOCK format
    kg=global_z_indices(START,YBLOCK,myid) + k - 1
    xx = 1.d0 - exp( -((z(kg)-zb)/sigma_z)**2 ) - exp( -((z(kg)-zt)/sigma_z)**2 )
    rhs(:,:,k) = rhs(:,:,k) - Ri*pd(:,:,k)*xx
   enddo
  else   
   rhs(:,:,:) = rhs(:,:,:) - Ri*pd(:,:,:)
  endif
  
 endif
 
 !-----------------------------------------------------------------
 ! rhs <==  rhs + 'f~' * u  ; Rot(3) is "dimensionless f~"
 !   take f~->0 az z-> rigid lids, 
 !   this sets bcs on w* and allows homogeneous bcs for pressure
 !-----------------------------------------------------------------
 if( trim(f_plane_key)=='non-traditional' ) then
 
  if( deriv_type(3,3,1)=='sin' .and. deriv_type(1,3,1)=='cos' ) then  ! w=sin(z),u=cos(z)
   do k=1,array_size(KDIM,YBLOCK,myid)   ! z is dim3 in YBLOCK format
    kg=global_z_indices(START,YBLOCK,myid) + k - 1
    xx = 1.d0 - exp( -((z(kg)-zb)/sigma_z)**2 ) - exp( -((z(kg)-zt)/sigma_z)**2 )
    rhs(:,:,k) = rhs(:,:,k) + Rot(3)*u(:,:,k)*xx
   enddo
  else
   rhs(:,:,:) = rhs(:,:,:) + Rot(3)*u(:,:,:)
  endif
  
 endif
 
 !-----------------------------------------------------------------
 ! apply the generalized diffusion operator to w
 !  (use previously computed z derivative stored in tmpY(1,1,1,3))
 ! store the result in   diff  => implicit_rhs(:,:,:,id,N)
 ! NB
 !  diffusion uses tmpX(:,:,:,1-2),tmpY(:,:,:,1-2),tmpZ(:,:,:,1-2)
 !  internally for work space
 !-----------------------------------------------------------------
 call diffusion(w,tmpY(1,1,1,3),diff,id)
 
 
end subroutine explicit_rhs_3




subroutine explicit_rhs_4
 use decomposition_params
 use mpi_params,             only: myid,comm,ierr
 use etc,                    only: MM0,N
 use methods_params,         only: deriv_type,do_nonlinear,do_immersed_boundary
 use dependent_variables,    only: u,v,w,s1,s1_bar 
 use intermediate_variables, only: tmpY,div_u,    &
                                   explicit_rhs,  &
                                   implicit_rhs
 use immersed_boundary,      only: nhat
                                   
 implicit none  
 integer                                 :: order
 integer                                 :: i,j,k,k_global,ig,jg
 integer,parameter                       :: id=4
 integer,parameter                       :: sign=-1
 integer,parameter                       :: uid=1,vid=2,wid=3
 integer,parameter                       :: xdir=1,ydir=2,zdir=3
 character(len=80),dimension(3)          :: method
 
 real(kind=8),dimension(:,:,:),pointer   :: rhs      ! -u dot grad s1 - w*ddz(s1_bar) 
 real(kind=8),dimension(:,:,:),pointer   :: diff     ! generalized {D[s1]}
 

 rhs    => explicit_rhs(:,:,:,id,MM0)
 diff   => implicit_rhs(:,:,:,id,N)

 
 !-----------------------------------------------------------------
 ! div(u)  (in a perfect world, it would be zero...)
 !          store results in div_u(:,:,:)
 !-----------------------------------------------------------------
 !order=1
 !method(1)=deriv_type(uid,xdir,order)  ! for 1st deriv of u wrt x
 !method(2)=deriv_type(vid,ydir,order)  ! for 1st deriv of v wrt y
 !method(3)=deriv_type(wid,zdir,order)  ! for 1st deriv of w wrt z
 !call divergence(u,v,w,div_u,method)

 
 !-----------------------------------------------------------------
 ! grad s1 
 !        store results in tmpY(:,:,:,1-3)
 !-----------------------------------------------------------------
 order=1
 method(:)=deriv_type(id,:,order)
 call gradient(s1,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),method)
 

 !-----------------------------------------------------------------
 ! rhs <==  -(u dot grad s1)  (i.e. s1') 
 !-----------------------------------------------------------------
 if( do_nonlinear ) then
  call YBdotproduct(rhs,sign,u,v,w,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3))
 endif
 
 
 !-----------------------------------------------------------------
 ! advection of time independent background gradient ...
 ! rhs <==  rhs - w*ddz(s1_bar)
 !-----------------------------------------------------------------
 if( do_nonlinear ) then  ! add to existing rhs
  do k=1,array_size(KDIM,YBLOCK,myid)   ! z is dim3 in YBLOCK format
   k_global = global_z_indices(START,YBLOCK,myid) + k - 1
   rhs(:,:,k) = rhs(:,:,k) - w(:,:,k)*s1_bar(k_global,2)   
  enddo
 else   ! rhs hasn't been set yet
  do k=1,array_size(KDIM,YBLOCK,myid)   ! z is dim3 in YBLOCK format
   k_global = global_z_indices(START,YBLOCK,myid) + k - 1
   rhs(:,:,k) = - w(:,:,k)*s1_bar(k_global,2)  
  enddo
 endif
 
 
 !-----------------------------------------------------------------
 ! apply the generalized diffusion operator to s1
 !  (use previously computed z derivative stored in tmpY(1,1,1,3))
 ! store the result in   diff  => implicit_rhs(:,:,:,id,N)
 ! NB
 !  diffusion uses tmpX(:,:,:,1-2),tmpY(:,:,:,1-2),tmpZ(:,:,:,1-2)
 !  internally for work space
 !-----------------------------------------------------------------
 call diffusion(s1,tmpY(1,1,1,3),diff,id)
 
end subroutine explicit_rhs_4



subroutine explicit_rhs_5
 use decomposition_params
 use mpi_params,             only: myid 
 use etc,                    only: MM0,N
 use methods_params,         only: deriv_type,do_nonlinear,  &
                                   do_second_scalar,do_immersed_boundary
 use dependent_variables,    only: u,v,w,s2,s2_bar      
 use intermediate_variables, only: tmpY,div_u,    &
                                   explicit_rhs,  &
                                   implicit_rhs
 use immersed_boundary,      only: nhat
 
 implicit none  
 integer                                 :: order
 integer                                 :: i,j,k,k_global,ig,jg
 integer,parameter                       :: id=5
 integer,parameter                       :: sign=-1
 character(len=80),dimension(3)          :: method
 
 real(kind=8),dimension(:,:,:),pointer   :: rhs      ! -u dot grad s2 - w*ddz(s2_bar) 
 real(kind=8),dimension(:,:,:),pointer   :: diff     ! generalized {D[s2]}
 
 if( .NOT. do_second_scalar ) stop 'inside explicit_rhs w/ do_second_scalar set FALSE'
 !-----------------------------------------------------------------
 ! div(u) has been computed in explicit_rhs_4 and stored in div_u 
 !-----------------------------------------------------------------

 rhs    => explicit_rhs(:,:,:,id,MM0)
 diff   => implicit_rhs(:,:,:,id,N)

 
 !-----------------------------------------------------------------
 ! grad s2 
 !        store results in tmpY(:,:,:,1-3)
 !-----------------------------------------------------------------
 order=1
 method(:)=deriv_type(id,:,order)
 call gradient(s2,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),method)

 
 !-----------------------------------------------------------------
 ! rhs <==  -(u dot grad s2)
 !-----------------------------------------------------------------
 if( do_nonlinear ) then
  call YBdotproduct(rhs,sign,u,v,w,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3))
 endif
 
 
 !-----------------------------------------------------------------
 ! rhs <==  rhs - w*ddz(s2_bar)
 !-----------------------------------------------------------------
 if( do_nonlinear ) then
  do k=1,array_size(KDIM,YBLOCK,myid)   ! z is dim3 in YBLOCK format
   k_global = global_z_indices(START,YBLOCK,myid) + k - 1
   rhs(:,:,k) = rhs(:,:,k) - w(:,:,k)*s2_bar(k_global,2)   
  enddo
 else
  do k=1,array_size(KDIM,YBLOCK,myid)   ! z is dim3 in YBLOCK format
   k_global = global_z_indices(START,YBLOCK,myid) + k - 1
   rhs(:,:,k) = - w(:,:,k)*s2_bar(k_global,2)   
  enddo
 endif
 
 !-----------------------------------------------------------------
 ! apply the generalized diffusion operator to s2
 !  (use previously computed z derivative stored in tmpY(1,1,1,3))
 ! store the result in   diff  => implicit_rhs(:,:,:,id,N)
 ! NB
 !  diffusion uses tmpX(:,:,:,1-2),tmpY(:,:,:,1-2),tmpZ(:,:,:,1-2)
 !  internally for work space
 !-----------------------------------------------------------------
 call diffusion(s2,tmpY(1,1,1,3),diff,id)
  

end subroutine explicit_rhs_5



subroutine diffusion(f,fz,diff_f,fid)
!------------------------------------------------------------
!  Compute generalized D[f] (D = diffusion operator).
!  Input and output data arrays arranged in YBLOCK format.
!------------------------------------------------------------
 use mpi_params,             only: myid,comm,ierr
 use decomposition_params
 use diffusion_params
 use methods_params,         only: deriv_type
 use intermediate_variables, only: tmpX,tmpY,tmpZ
 
 implicit none
 integer                       ::  order,dir,fid
 character(len=80)             ::  method
 integer,save                  ::  i,j,k,n,isgn(3)
 logical,save                  ::  first_entry=.TRUE.
 real(kind=8),parameter        ::  dzero=0.d0
 real(kind=8)                  ::  f( array_size(IDIM,YBLOCK,myid),   &
                                      array_size(JDIM,YBLOCK,myid),   &
                                      array_size(KDIM,YBLOCK,myid)  )
                                       
 real(kind=8)                  ::  fz( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid)  )
                                       
 real(kind=8)                  :: diff_f( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid)  )
 if( first_entry ) then
  n = array_size(IDIM,YBLOCK,myid)   &
     *array_size(JDIM,YBLOCK,myid)   &
     *array_size(KDIM,YBLOCK,myid)
  isgn(1) = (-1)**(p(1)-1)
  isgn(2) = (-1)**(p(2)-1)
  isgn(3) = (-1)**(p(3)-1)
  first_entry=.FALSE.
 endif
 
 !---------------------------------------------------
 !  transpose f to ZBLOCK format
 !---------------------------------------------------
  call yblock_2_zblock(f,tmpZ)  
   
 !---------------------------------------------------
 !  take z deriv of result
 !---------------------------------------------------
  order=1          ! NB methods for orders>1 are determined internally                   
                   ! triggering off the method for order=1
  dir=3            ! ddz
  method = deriv_type(fid,dir,order)
  
 ! need to pass along the actual order once the method is determined
  order=2*p(3)    ! actual order depends on p(3), start w/ f
  
  call ddz(tmpZ(1,1,1,1),tmpZ(1,1,1,2),method,order)
  
 !--------------------------------------------------------
 !  transpose d/dz( alpha*dfdz ) back to to YBLOCK format
 !       store partial result in tmpY(:,:,:,1)
 !--------------------------------------------------------
  call zblock_2_yblock(tmpZ(1,1,1,2),tmpY(1,1,1,1))
  
 
 if( array_size(IDIM,XBLOCK,myid) > 1 ) then
  !---------------------------------------------------
  !  transpose f to XBLOCK format
  !---------------------------------------------------
   call yblock_2_xblock(f,tmpX)
  
  !---------------------------------------------------
  !  take appropriate x deriv of result
  !---------------------------------------------------
   order=1      ! for x & y higher deriv methods determined
                ! internally based on 1st order method
   dir=1        ! ddx
   method = deriv_type(fid,dir,order)
   order=2*p(1) ! now specify which deriv we want
   call ddx(tmpX(1,1,1,1),tmpX(1,1,1,2),method,order)
  
  !--------------------------------------------------------
  !  transpose resulting deriv back to to YBLOCK format
  !--------------------------------------------------------
   call xblock_2_yblock(tmpX(1,1,1,2),tmpY(1,1,1,2))
  else
   call dinit(n,dzero,tmpY(1,1,1,2))
  endif
 
 if( array_size(IDIM,YBLOCK,myid) > 1 ) then
  !---------------------------------------------------
  !  take appropriate y deriv of f
  !---------------------------------------------------
   order=1      ! for x & y higher deriv methods determined
                ! internally based on 1st order method
   dir=2        ! ddy
   method = deriv_type(fid,dir,order)
   order=2*p(2) ! now specify which deriv we want
   call ddy(f,diff_f,method,order)
  else
   call dinit(n,dzero,diff_f)
  endif
  
!---------------------------------------------------
!  weight and add the partial results
!---------------------------------------------------

   diff_f(:,:,:) = isgn(1)*diff_coeffs(1,fid)*tmpY(:,:,:,2)  &    ! x term
          + isgn(2)*diff_coeffs(2,fid)*diff_f(:,:,:)         &    ! y term
          + isgn(3)*diff_coeffs(3,fid)*tmpY(:,:,:,1)              ! z term
          

 return
end subroutine diffusion




