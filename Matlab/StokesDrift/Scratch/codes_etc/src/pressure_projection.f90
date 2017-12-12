subroutine pressure_projection 
 use methods_params,           only: deriv_type,do_immersed_boundary
 use intermediate_variables,   only: div_u
 use dependent_variables,      only: u,v,w
 use etc,                      only: istep
 use mpi_params,               only: myid
 use dimensional_scales,       only: velocity_scale,length_scale
 
 implicit none
 real(kind=8),save                     :: maxdiv,last_maxdiv
 integer                               :: i
 integer,parameter                     :: id=6  ! phi
 integer,parameter                     :: uid=1,vid=2,wid=3
 integer,parameter                     :: xdir=1,ydir=2,zdir=3
 integer,parameter                     :: order=1
 character(len=80),dimension(3),save   :: method
 logical,save                          :: first_entry=.TRUE.

 if( first_entry ) then
  !-------------------------------------------------------------
  ! determine differentiation methods
  !-------------------------------------------------------------
  method(1)=deriv_type(uid,xdir,order)
  method(2)=deriv_type(vid,ydir,order)
  method(3)=deriv_type(wid,zdir,order)  
  first_entry=.FALSE.
 endif
 
   
  !-------------------------------------------------------------
  ! compute divergence (u*)
  !-------------------------------------------------------------
  call divergence(u,v,w,div_u,method)

      
  !-------------------------------------------------------------
  ! solve the poisson eqn using implicit solver 
  !-------------------------------------------------------------
  call implicit_solve(id)
 
  !-------------------------------------------------------------
  ! compute the gradient of 'pressure' & project
  !-------------------------------------------------------------
  call project
 
! call zero_flow_below_IB
  
 
end subroutine pressure_projection


subroutine project
!-------------------------------------------------------------
!-   u <- u* -dt*GradPhi(1)  
!-   v <- v* -dt*GradPhi(2)
!-   w <- w* -dt*GradPhi(3)
!-------------------------------------------------------------
 use decomposition_params
 use dependent_variables,             only: u,v,w,s1,s2
 use intermediate_variables,          only: phi,tmpY
 use independent_variables,           only: dt,nz 
 use methods_params,                  only: deriv_type,do_immersed_boundary
 use mpi_params,                      only: myid
 
 implicit none
 
 integer                                      :: k,id
 integer,parameter                            :: order=1
 real(kind=8),parameter                       :: dzero=0.d0
 
 
 character(len=80),dimension(3),save          :: method
 real(kind=8),save                            :: alpha
 real(kind=8),dimension(:,:,:),pointer,save   :: phi_x
 real(kind=8),dimension(:,:,:),pointer,save   :: phi_y
 real(kind=8),dimension(:,:,:),pointer,save   :: phi_z
 integer,save                                 :: npts
 logical,save                                 :: first_entry=.TRUE.
 
 if( first_entry ) then
  !------------------------------------------------------------
  ! pointers, method, -dt, npts fixed in time
  !------------------------------------------------------------
  phi_x => tmpY(:,:,:,1)
  phi_y => tmpY(:,:,:,2)
  phi_z => tmpY(:,:,:,3)
  id=6    ! get methods for pressure
  method(:)=deriv_type(id,:,order)
  alpha = -dt
  npts = array_size(IDIM,YBLOCK,myid)*array_size(JDIM,YBLOCK,myid)
  first_entry=.FALSE.
 endif
 
 !if( do_immersed_boundary ) then
  !call homogeneous_neumann_blend(phi)
 !endif
 
 !------------------------------------------------------------
 ! compute grad phi
 !------------------------------------------------------------ 
 call gradient(phi,phi_x,phi_y,phi_z,method) 
 
 !------------------------------------------------------------
 ! update the velocity components
 !------------------------------------------------------------ 
 call axpy(alpha,phi_x,u)  ! u <== u - dt phi_x
 call axpy(alpha,phi_y,v)  ! v <== v - dt phi_y
 call axpy(alpha,phi_z,w)  ! w <== w - dt phi_z
  
 return  ! what follows is not necessary when bcs satisfied implicitly
 !----------------------------------------------------------- 
 ! except when z coordinate is periodic, 
 ! top and bottom are rigid lids --> w=0 there
 !----------------------------------------------------------- 
 id=3   !  w
 if( trim( deriv_type(id,3,order) ) == 'sin' ) then
  if( global_z_indices(START,YBLOCK,myid) == 1 ) then
   k=1
   call dinit(npts,dzero,w(1,1,k))  ! w(:,:,k) = 0.d0    
  endif
 endif
 if( trim( deriv_type(id,3,order) ) == 'sin' ) then
  if( global_z_indices(END,YBLOCK,myid) == nz ) then
   k=array_size(KDIM,YBLOCK,myid)
   call dinit(npts,dzero,w(1,1,k))  ! w(:,:,k) = 0.d0    
  endif
 endif
 
  
end subroutine project


