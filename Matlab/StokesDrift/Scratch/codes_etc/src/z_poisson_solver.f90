!------------------------------------------------------------
! driver routine for implicit solves when both x and y
! are fourier sin or cos expanded over regularly spaced
! x-y grid, i.e. no stretching in either x or y. Both x/y
! dimensions are treated implicitly leaving 1d solves in
! the z direction
!------------------------------------------------------------
subroutine implicit_solve(id) 
!--------------------------------------------------------------- 
!
! After xy transform and AM discretization:
! 
!    xx * phi_zz + q phi = R
!
!    xx = alpha3/beta  for id=1,2,...5
!         1.d0             id = 6
!
!         alpha3 = diff_coeffs(3,id)
!         beta   = -2/dt (AM3) or -12/(5dt) AM4
!
!    q  = 1.d0 + q0/beta    for id=1,2,...5
!         q1                    id = 6
!
!         q0 = -1/alpha1 * (k)^2p1 -1/alpha2 * (k)^2p2
!
!         q1 = -k^2 - l^2  (initialized in setup_diffusion.f90)
!
!    R  =  u~ - (1/beta)*D(u~)   e.g. for u, AM3 
!
!          -1/dt div(u*)         pressure
!
!          u~  the field after explicit step, prior to diffusion
!          D(u~)  1/Re grad^2 u at previous time step e.g.
!
! On entry the implicit_rhs(:,:,:,id,N) terms are filled with
!
!         'generalized' D[u] etc at t_n etc
!---------------------------------------------------------------
 use dependent_variables
 use intermediate_variables
 use decomposition_params
 use diffusion_params,        only: diff_coeffs
 use independent_variables,   only: dt,npts=>nz,ny,nx
 use pde_params,              only: Rot
 use methods_params,          only: AM_order,deriv_type
 use etc,                     only: N,NM1,istep,runlabel
 use mpi_params,              only: myid,numprocs,comm,ierr
 
 implicit none 
 real(kind=8)                    :: beta
 real(kind=8)                    :: q
 real(kind=8)                    :: rhs_max
 real(kind=8)                    :: xx
 real(kind=8),allocatable,save   :: wrk(:)
 real(kind=8),parameter          :: dzero=0.d0
 integer                         :: id
 integer                         :: order
 integer                         :: iend
 integer                         :: dim
 integer                         :: i,j,k
 integer,parameter               :: FOR=1,INV=-1,inc=1
 character(len=80)               :: exp_type(2)
 character(len=80)               :: z_deriv_type
 logical,save                    :: first_entry=.TRUE.
 
 real(kind=8),dimension(:,:,:),pointer  :: f
 real(kind=8),dimension(:,:,:),pointer  :: rhs
 real(kind=8),dimension(:,:,:),pointer  :: rhs_hat
 real(kind=8),dimension(:,:,:),pointer  :: rhs_ZBLOCK
 real(kind=8),dimension(:,:,:),pointer  :: soln_ZBLOCK
 real(kind=8),dimension(:,:,:),pointer  :: Diff_f_n
 real(kind=8),dimension(:,:,:),pointer  :: Diff_f_nm1

 if( first_entry ) then
  allocate( wrk(npts) )
  first_entry=.FALSE.
 endif
 
 !---------------------------------------------------
 !  point at variable of interest
 !---------------------------------------------------
 if( id==1 ) then
  f => u
 elseif( id==2 ) then
  f => v
 elseif( id==3 ) then
  f => w
 elseif( id==4 ) then
  f => s1
 elseif( id==5 ) then
  f => s2
 elseif( id==6 ) then
  f => phi
 endif 
 
 
 if(id<6) then
  rhs  => tmpY(:,:,:,2)                      ! rhs in YBLOCK format
 elseif(id==6) then
  rhs => div_u                               ! rhs in YBLOCK format
 endif
 
 rhs_hat => tmpY(:,:,:,3)                    ! xy transform of rhs
 rhs_ZBLOCK  => tmpZ(:,:,:,1)                ! ZBLOCK version of rhs_hat
 soln_ZBLOCK => tmpZ(:,:,:,2)                ! ZBLOCK format soln

 
 !---------------------------------------------------
 !  set expansion type in x & y
 !  and 1st z derivative type
 !  decide whether direct spectral 
 !  solution is possible
 !---------------------------------------------------
  order=1
  exp_type(:)=deriv_type(id,1:2,order)   ! 1,2 -> x,y directions
  z_deriv_type=deriv_type(id,3,order)    ! 3 -> z direction
  
    
 !-----------------------------------------------
 ! for 2d problems in yz plane w/ no rotation 
 !  ==> u is identically zero at all times
 !-----------------------------------------------
  if( id==1 .and. nx==1 .and. Rot(1) == 0.d0 ) then
   f=0.d0
   return
  endif



  dim=3
  if(   trim(deriv_type(id,dim,1)) == 'fourier'  .or.   &
        trim(deriv_type(id,dim,1)) == 'cos'      .or.   &
        trim(deriv_type(id,dim,1)) == 'sin' ) then       
  else
   if(myid==0) write(0,*) 'SHOULD BE SPECTRAL TYPE SCHEME IN Z'
   stop
  endif 

 
 !---------------------------------------------------
 !  set beta and rhs in YBLOCK format
 !  depending on id and implicit integration method
 !---------------------------------------------------
 if( AM_order==3 ) then
   
  if( id < 6 ) then
   beta = -2.d0/dt 
   xx = 1.d0/beta
   Diff_f_n => implicit_rhs(:,:,:,id,N)         ! D[f] at N stored here
   rhs(:,:,:) = f(:,:,:) - xx*Diff_f_n(:,:,:)
  endif
  
 elseif( AM_order==4 ) then
 
  if( id < 6 ) then
   beta = -12.d0/(5.d0*dt)
   xx = 1.d0/beta
   Diff_f_n => implicit_rhs(:,:,:,id,N)          ! D[f] at N stored here
   Diff_f_nm1 => implicit_rhs(:,:,:,id,NM1)      ! D[f] at N-1 stored here
   rhs(:,:,:) = f(:,:,:) - xx*(8.d0/5.d0)*Diff_f_n(:,:,:) + xx*(1.d0/5.d0)*Diff_f_nm1(:,:,:)

  endif
  
  
 else
  stop ' problem w/ AM_order in implicit_solve '
 endif
 
 if( id==6 ) then
  xx = 1.d0/dt
  call scale(xx,rhs)  ! rhs = (1./dt)*div_u
 endif
!---------------------------------------------------
!  ====> rhs now set in YBLOCK format
!---------------------------------------------------

 !--------------------------------------------------
 ! transform the entire 3d rhs array
 ! in both x & y directions
 ! (transform_xy uses tmpY(:,:,:,1) as wrk space)
 !--------------------------------------------------
 call transform_xy(rhs,rhs_hat,FOR,exp_type)

 !--------------------------------------------------
 ! transpose the transformed rhs data to
 ! ZBLOCK format, store in rhs_ZBLOCK
 !--------------------------------------------------
 call yblock_2_zblock(rhs_hat,rhs_ZBLOCK)


!----------------------------------------------------
! loop over kx,ky vals and do solves
!----------------------------------------------------
 do i=1,array_size(JDIM,ZBLOCK,myid)  ! x/kx indices
  do j=1,array_size(KDIM,ZBLOCK,myid) ! y/ky indices

   !------------------------
   !  set Helmholtz term
   !------------------------
   if(id < 6 ) then
    q = q0(i,j,id)/beta + 1.d0   ! transformed horizontal part of (generalized D + beta)/beta
   elseif( id==6 ) then
    !----------------------------------------
    !  for pressure, 'beta=0' and 'alpha3'=1
    !----------------------------------------
    q = q1(i,j)             ! transformed horizontal part of laplacian, beta=0
   endif
   
   !----------------------------------------
   ! do the solve...
   !----------------------------------------
   call spectral_z_solver(npts,q,soln_ZBLOCK(:,i,j),rhs_ZBLOCK(:,i,j),id,wrk(1)) 
       
   
   !--------------------------------------------------------
   ! For rigid lids at z=0,Lz ... deriv_type(var,dir,order)
   ! at all z levels, set the xy mean of w* and phi to zero
   ! ===> <w> = <w*> -dt*< phi_z > will be identically zero
   !  k=l=0 easily identified via q1=0
   !
   !    bad idea for IB sidewalls, or no-slip sidewalls
   !
   !--------------------------------------------------------
   if( trim(runlabel) == 'hc_temp_bcs' ) then
    if( trim( deriv_type(3,3,1) ) .eq. 'sin' ) then
     if( q1(i,j) == 0.d0 ) then
      if( id==3 .or. id==6 ) then
       soln_ZBLOCK(:,i,j)=0.d0
      endif
     endif
    endif
   endif
   
   
   
   !--------------------------------------------------------
   ! set all solns to zero at x and y nyquist wavenumbers
   ! (for which q1 = 999.d0, generally q1=-(k^2+l^2)<=0 )
   !   actually, I'm no longer setting nyquist to 999
   !   in setup diffusion.f90
   !--------------------------------------------------------
   if( q1(i,j) > 0 ) then
    soln_ZBLOCK(:,i,j)=0.d0 
   endif

  enddo
 enddo


!---------------------------------
! transpose soln to YBLOCK format
! store result in tmpY(:,:,:,1)
! overwrite f with the result
!---------------------------------
 call zblock_2_yblock(soln_ZBLOCK,tmpY(:,:,:,1))
 
!---------------------------------
! inv xy transform the soln
! final result  ===> f in YBLOCK format
!---------------------------------
 call transform_xy(tmpY(:,:,:,1),f,INV,exp_type)
 
 return
end subroutine implicit_solve

subroutine spectral_z_solver(n,q,y,f,id,wrk)
!--------------------------------------------!
!                                            !
!     ode,bcs:  L[y(z)] = f    z in [0,L]    !
!                                            !
!     xx * d2/dz2 y + q*y =  f               !
!                                            !
!     xx = alpha3/beta  id = 1,2,...5        !
!        = 1.d0         id = 6               ! 
!                                            !
!     z regularly spaced                     !
!                                            !
!     bcs  consistent w/ sin/cos or fourier  !
!                                            !
!--------------------------------------------!
 use mpi_params,             only: myid
 use methods_params,         only: deriv_type,AM_ORDER
 use independent_variables,  only: nz,dt
 use diffusion_params,       only: diff_coeffs,p
 use differentiation_params


 implicit none
 integer                        :: id,i,j,n,ii
 integer,parameter              :: idim=3,for=1,inv=2,order=1
 
 real(kind=8)                   :: q,xx,beta,yy
 real(kind=8)                   :: y(n)
 real(kind=8)                   :: f(n)
 real(kind=8)                   :: wrk(n)
 real(kind=8)                   :: xnorm
 character(len=80)              :: method
 integer(kind=8)                :: plans(2)

 
 if( AM_order==3 ) then
  beta = -2.d0/dt
 elseif( AM_order==4 ) then
  beta = -12.d0/(5.d0*dt)
 endif
 


 !---------------------------------------
 ! get the method from the deriv type
 ! and use to select plans, normalization
 !---------------------------------------
 method = deriv_type(id,idim,order)
 if( trim(method)=='fourier' ) then
  plans(1) = fourier_plan(idim,for)
  plans(2) = fourier_plan(idim,inv)
  xnorm = 1.d0/dfloat(n)
  ii=1                   ! transform full data array
 elseif( trim(method)=='cos' ) then
  plans(1) = cos_plan(idim)
  plans(2) = cos_plan(idim)
  xnorm = 1.d0/(2.d0*(n-1.d0))
  ii=1                   ! transform full data array
 elseif( trim(method)=='sin' ) then
  plans(1) = sin_plan(idim)
  plans(2) = sin_plan(idim)
  xnorm = 1.d0/(2.d0*(n-1.d0))
  ii=2                   ! ignore zero endpts for sin transforms
  wrk(1)=0.d0
  wrk(n)=0.d0            ! make sure there are no bad vals left over
 endif

 
 !---------------------------------------
 !  take forward transform of the rhs f
 !---------------------------------------
 call dfftw_execute_r2r(plans(1),f(ii),wrk(ii))      

 if( id < 6 ) then
  yy = diff_coeffs(3,id)/beta
  do i=1,n
   xx = (q - yy*kz(i)**(2*p(3)))    !! could be high order diffusion in z
   wrk(i) = xnorm*wrk(i)/xx
  enddo
 elseif( id == 6 ) then
  do i=2,n
   xx = (q - kz(i)**2)
   wrk(i) = xnorm*wrk(i)/xx
  enddo
 endif
 
 !-----------------------------------------
 ! The only special cases are for pressure
 !-----------------------------------------
 if( id==6 ) then
  if( q==0.d0 ) then         ! kz=0 and q=0
   wrk(1)=0.d0
  else                       ! kz=0 but q nonzero
   wrk(1) = xnorm*wrk(1)/q
  endif
 endif
 
 !---------------------------------------
 !  take inverse transform of result
 !---------------------------------------
 call dfftw_execute_r2r(plans(2),wrk(ii),y(ii))
 
 if( trim(method)=='sin' ) then
  y(1)=0.d0
  y(n)=0.d0
 endif

! q is ALWAYS nonzero for u,v,w,s1,s2 for ALL k,l
! q=0 for phi when k=l=0

return
end subroutine spectral_z_solver


