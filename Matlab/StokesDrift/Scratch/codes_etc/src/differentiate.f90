subroutine gradient(f,fx,fy,fz,method)
!------------------------------------------------------------
!  Compute spatial gradient of f using derivative methods
!  specified in method(1:3). Input and output data arrays
!  arranged in YBLOCK format.
!------------------------------------------------------------
 use decomposition_params
 use mpi_params,             only: myid
 use intermediate_variables, only: tmpX,tmpZ
 implicit none
 integer,parameter             ::   order=1
 real(kind=8),parameter        ::   dzero=0.d0
 integer,save                  ::   n_xb,n_yb,n_zb
 logical,save                  ::   first_entry=.TRUE.
 character(len=80)             ::   method(3)
 real(kind=8)                  ::   f( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid)  )
                                       
 real(kind=8)                  ::  fx( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid)  )
                                       
 real(kind=8)                  ::  fy( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid)  )
                                       
 real(kind=8)                  ::  fz( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid)  )    
 if( first_entry ) then
  n_xb = array_size(IDIM,XBLOCK,myid)  &
        *array_size(JDIM,XBLOCK,myid)  &
        *array_size(KDIM,XBLOCK,myid)
        
  n_yb = array_size(IDIM,YBLOCK,myid)  &
        *array_size(JDIM,YBLOCK,myid)  &
        *array_size(KDIM,YBLOCK,myid)
        
  n_zb = array_size(IDIM,ZBLOCK,myid)  &
        *array_size(JDIM,ZBLOCK,myid)  &
        *array_size(KDIM,ZBLOCK,myid)
  first_entry=.FALSE.
 endif
 
 if( array_size(IDIM,YBLOCK,myid) > 1 ) then
  call ddy(f,fy,method(2),order)
 else
  call dinit(n_yb,dzero,fy)
 endif
 
 if( array_size(IDIM,ZBLOCK,myid) > 1 ) then
  call yblock_2_zblock(f,tmpZ)
  call ddz(tmpZ,tmpZ(1,1,1,2),method(3),order)
  call zblock_2_yblock(tmpZ(1,1,1,2),fz)
 else
  call dinit(n_zb,dzero,fz)
 endif
 
 if( array_size(IDIM,XBLOCK,myid) > 1 ) then
  call yblock_2_xblock(f,tmpX)
  call ddx(tmpX,tmpX(1,1,1,2),method(1),order) 
  call xblock_2_yblock(tmpX(1,1,1,2),fx)
 else
  call dinit(n_xb,dzero,fx)
 endif
 
 return
end subroutine gradient

subroutine divergence(u,v,w,div,method)
!------------------------------------------------------------
!  Compute divergence of [u,v,w] using derivative methods
!  specified in method(1:3). Input and output data arrays
!  arranged in YBLOCK format.
!------------------------------------------------------------
 use decomposition_params
 use mpi_params,             only: myid
 use intermediate_variables, only: tmpX,tmpY,tmpZ
 
 implicit none
 integer                       ::   i,j,k
 integer,parameter             ::   order=1
 real(kind=8),parameter        ::   dzero=0.d0
 integer,save                  ::   n
 logical,save                  ::   first_entry=.TRUE.
 character(len=80)             ::   method(3)
 real(kind=8)                  ::   u( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid)  )
                                       
 real(kind=8)                  ::   v( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid)  )
                                       
 real(kind=8)                  ::   w( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid)  )
                                       
 real(kind=8)                  :: div( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid)  )    
 
 if( first_entry ) then
  n = array_size(IDIM,YBLOCK,myid)  &
     *array_size(JDIM,YBLOCK,myid)  &
     *array_size(KDIM,YBLOCK,myid)
  first_entry=.FALSE.
 endif
 
 !----------------------------------------
 !    div <--- dv/dy
 !----------------------------------------
 if( array_size(IDIM,YBLOCK,myid) > 1 ) then
  call ddy(v,div,method(2),order)               ! ==> dv/dy in YBLOCK format
 else
  call dinit(n,dzero,div)
 endif
 
 !----------------------------------------
 !    div <--- dv/dy + dw/dz
 !----------------------------------------
 if( array_size(IDIM,ZBLOCK,myid) > 1 ) then
  call yblock_2_zblock(w,tmpZ)
  call ddz(tmpZ,tmpZ(1,1,1,2),method(3),order)       ! ==> dw/dz in ZBLOCK format
  call zblock_2_yblock(tmpZ(1,1,1,2),tmpY(1,1,1,2))  ! ==> dw/dz in YBLOCK format

  div(:,:,:) = div(:,:,:) + tmpY(:,:,:,2)



 endif
 
 !----------------------------------------
 !    div <--- dv/dy + dw/dz + dudx
 !----------------------------------------
 if( array_size(IDIM,XBLOCK,myid) > 1 ) then
  call yblock_2_xblock(u,tmpX)
  call ddx(tmpX,tmpX(1,1,1,2),method(1),order)  ! ==> du/dx in XBLOCK format
  call xblock_2_yblock(tmpX(1,1,1,2),tmpY(1,1,1,1))  ! ==> du/dx in YBLOCK format

  div(:,:,:) = div(:,:,:) + tmpY(:,:,:,1)

 endif
 
 
 return
end subroutine divergence



subroutine ddx(f,df,method,order)
!-------------------------------------------------------------------------
! assume f and df are decomposed in XBLOCK decomposition
!-------------------------------------------------------------------------
 use independent_variables,  only: nx
 use mpi_params,             only: myid
 use decomposition_params,   only: array_size,IDIM,JDIM,KDIM,XBLOCK 
 implicit none 
 integer                        ::  i,j,k
 real(kind=8)                   ::  f( array_size(IDIM,XBLOCK,myid),    &
                                       array_size(JDIM,XBLOCK,myid),    &
                                       array_size(KDIM,XBLOCK,myid)  )
 real(kind=8)                   :: df( array_size(IDIM,XBLOCK,myid),    &
                                       array_size(JDIM,XBLOCK,myid),    &
                                       array_size(KDIM,XBLOCK,myid)  )
 real(kind=8),allocatable,save  :: tmp(:)
 real(kind=8),parameter         :: dzero=0.d0
 character(len=80)              :: method
 integer                        :: order
 integer,save                   :: npts
 logical,save                   :: first_entry=.TRUE.
!--------------------------------------------------------------------------

if(first_entry) then
 if( nx .ne. array_size(IDIM,XBLOCK,myid) ) stop 'XBLOCK decomposition error, ddx'
 if(nx>1) allocate( tmp(nx) )
 npts = array_size(IDIM,XBLOCK,myid)   &
       *array_size(JDIM,XBLOCK,myid)   &
       *array_size(KDIM,XBLOCK,myid)
 first_entry=.FALSE.
endif

if(nx==1) then
 call dinit(npts,dzero,df)
 return
endif

!--------------------------------------------------------------------------
! loop through 2nd & 3rd array indices, perform stride 1
! global differentiation operation, with data local to myid
!--------------------------------------------------------------------------
  
 do k=1,array_size(KDIM,XBLOCK,myid)  ! distribute threads over z loop
  do j=1,array_size(JDIM,XBLOCK,myid)  ! loop over y=index 2 in XBLOCK decomp  
   call differentiate(f(1,j,k),df(1,j,k),IDIM,method,tmp(1),order)
  enddo
 enddo



end subroutine ddx



subroutine ddy(f,df,method,order)
!-------------------------------------------------------------------------
! assume f and df are decomposed in YBLOCK decomposition
!-------------------------------------------------------------------------
 use independent_variables,  only: ny
 use mpi_params,             only: myid
 use decomposition_params
 implicit none 
 integer                        :: i,j,k
 real(kind=8)                   ::  f( array_size(IDIM,YBLOCK,myid),    &
                                       array_size(JDIM,YBLOCK,myid),    &
                                       array_size(KDIM,YBLOCK,myid)  )
 real(kind=8)                   :: df( array_size(IDIM,YBLOCK,myid),    &
                                       array_size(JDIM,YBLOCK,myid),    &
                                       array_size(KDIM,YBLOCK,myid)  )
 real(kind=8),allocatable,save  :: tmp(:)
 real(kind=8),parameter         :: dzero=0.d0
 character(len=80)              :: method
 integer                        :: order 
 integer,save                   :: npts
 logical,save                   :: first_entry=.TRUE.
!--------------------------------------------------------------------------

if(first_entry) then
 if( ny .ne. array_size(IDIM,YBLOCK,myid) ) stop 'YBLOCK decomposition error, ddy'
 if(ny>1) allocate( tmp(ny) )
 npts = array_size(IDIM,YBLOCK,myid)   &
       *array_size(JDIM,YBLOCK,myid)   &
       *array_size(KDIM,YBLOCK,myid)
 first_entry=.FALSE.
endif

if(ny==1) then
 call dinit(npts,dzero,df)
 return
endif

!--------------------------------------------------------------------------
! loop through 2nd & 3rd array indices, perform stride 1
! global differentiation operation, with data local to myid
!--------------------------------------------------------------------------

 do k=1,array_size(KDIM,YBLOCK,myid)  ! distribute threads over z loop
  do j=1,array_size(JDIM,YBLOCK,myid)  ! loop over x=index 2 in YBLOCK decomp
   call differentiate(f(1,j,k),df(1,j,k),JDIM,method,tmp(1),order)   
 enddo
 enddo


end subroutine ddy



subroutine ddz(f,df,method,order)
!-------------------------------------------------------------------------
! assume f and df are decomposed in ZBLOCK decomposition
!-------------------------------------------------------------------------
 use independent_variables,  only: nz
 use mpi_params,             only: myid
 use decomposition_params,   only: array_size,IDIM,JDIM,KDIM,ZBLOCK 
 implicit none 
 integer                        :: i,j,k
 real(kind=8)                   ::  f( array_size(IDIM,ZBLOCK,myid),    &
                                       array_size(JDIM,ZBLOCK,myid),    &
                                       array_size(KDIM,ZBLOCK,myid)  )
 real(kind=8)                   :: df( array_size(IDIM,ZBLOCK,myid),    &
                                       array_size(JDIM,ZBLOCK,myid),    &
                                       array_size(KDIM,ZBLOCK,myid)  )
 real(kind=8),allocatable,save  :: tmp(:)
 real(kind=8),parameter         :: dzero=0.d0
 character(len=80)              :: method
 integer                        :: order
 integer,save                   :: npts
 logical,save                   :: first_entry=.TRUE.
!--------------------------------------------------------------------------

if(first_entry) then
 if( nz .ne. array_size(IDIM,ZBLOCK,myid) ) stop 'ZBLOCK decomposition error, ddz'
 if(nz>1) allocate( tmp(nz) )
 npts = array_size(IDIM,ZBLOCK,myid)   &
       *array_size(JDIM,ZBLOCK,myid)   &
       *array_size(KDIM,ZBLOCK,myid)
 first_entry=.FALSE.
endif

if(nz==1) then
 call dinit(npts,dzero,df)
 return
endif


!--------------------------------------------------------------------------
! loop through 2nd & 3rd array indices, perform stride 1
! global differentiation operation, with data local to myid
!--------------------------------------------------------------------------

 do k=1,array_size(KDIM,ZBLOCK,myid)  ! distribute threads over y loop
  do j=1,array_size(JDIM,ZBLOCK,myid)  ! loop over x=index 2 in ZBLOCK decomp
   call differentiate(f(1,j,k),df(1,j,k),KDIM,method,tmp(1),order)
  enddo
 enddo


end subroutine ddz





subroutine differentiate(f,df,dir,method,tmp,order)
 use independent_variables
 use differentiation_params
 use dimensional_scales,              only: length_scale
 implicit none 
 integer                                 :: dir       ! 1,2,3 for diff wrt to x,y,z coordinate
 integer                                 :: order
 real(kind=8)                            :: f(*)   
 real(kind=8)                            :: df(*)
 real(kind=8)                            :: tmp(*)
 character(len=80)                       :: method   ! method of 1st deriv, even when order>1
 
 integer                                 :: n 
 real(kind=8)                            :: L
 integer(kind=8)                         :: plans(2)
 character(len=80)                       :: exp_type
 logical                                 :: apply_chain_rule
 
 real(kind=8), dimension(:), pointer     :: k
 real(kind=8), dimension(:), pointer     :: kfilter
 real(kind=8), dimension(:), pointer     :: dsdx
 real(kind=8), dimension(:,:), pointer   :: A
 integer, dimension(:), pointer          :: ipiv
 
 
 !---------------------------------------------------
 ! can always point, even if some arrays are unused
 ! because the unused ones are allocated with len=1
 !---------------------------------------------------
 if(dir==1) then
  L=Lx/length_scale
  k => kx
  kfilter => kxfilter
  A => c6_Ax
  ipiv => c6_ipivx
  dsdx => s_x
  n=nx
 elseif(dir==2) then
  L=Ly/length_scale
  k => ky
  kfilter => kyfilter
  A => c6_Ay
  ipiv => c6_ipivy
  dsdx => s_y
  n=ny
 elseif(dir==3) then
  L=Lz/length_scale
  k => kz
  kfilter => kzfilter
  A => c6_Az
  ipiv => c6_ipivz
  dsdx => s_z
  n=nz
 endif

 !----------------------------------------------------------
 ! for input data of length 1, set deriv to zero and return
 !----------------------------------------------------------
 if( n==1 ) then
  df(1)=0.d0
  return
 endif
 
 !----------------------------------------------------------
 ! decide whether chain rule is needed
 !----------------------------------------------------------
 apply_chain_rule = stretch(dir)
 
 !----------------------------------------------------------
 !   take d/ds of f, store result in df
 !   note L is dimensionless spatial domain size
 !----------------------------------------------------------
 if( method=='fourier' ) then
 
  if( order>1 .and. apply_chain_rule) &
   stop 'differentiate error: calling fourier deriv with order>1 when chain rule required'
  exp_type='fourier'
  plans(1)=fourier_plan(dir,1)
  plans(2)=fourier_plan(dir,2)
  call fourier_deriv(f,df,n,order,exp_type,k,kfilter,tmp,plans)
                         
 elseif( method=='cos' ) then
 
  if( order>1 .and. apply_chain_rule) &
   stop 'differentiate error: calling cos deriv with order>1 when chain rule required'
  exp_type='cos'
  plans(1)=cos_plan(dir)
  plans(2)=sin_plan(dir)
  call fourier_deriv(f,df,n,order,exp_type,k,kfilter,tmp,plans)
  
 elseif( method=='sin' ) then
 
  if( order>1 .and. apply_chain_rule) &
   stop 'differentiate error: calling sin deriv with order>1 when chain rule required'
  exp_type='sin'
  plans(1)=sin_plan(dir)
  plans(2)=cos_plan(dir)
  call fourier_deriv(f,df,n,order,exp_type,k,kfilter,tmp,plans)
  
 elseif( method=='compact' ) then
 
  if( order>1 ) stop 'differentiate error: order>1 not implemented for compact scheme'
  call compact_deriv(f,df,ds(dir),n,A,ipiv,order)
  
 elseif( method=='cheby' ) then
 
  !if( order>1 ) stop 'differentiate error: order>1 not implemented for cheby scheme'
  plans(1)=cheby_plan(dir)
  call cheby_deriv(f,df,L,n,order,tmp,plans(1))
  
 endif
 
 
 if( apply_chain_rule ) then   
  !---------------------------------------------------
  !   df/dx = df/ds * ds/dx etc.
  !---------------------------------------------------
  df(1:n) = df(1:n) * dsdx(1:n)
 endif
  
end subroutine differentiate


