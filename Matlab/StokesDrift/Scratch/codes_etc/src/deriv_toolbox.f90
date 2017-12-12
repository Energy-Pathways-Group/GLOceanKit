subroutine fourier_deriv(f,df,n,order,exp_type,kx,kfilter,tmp,plans)
!!------------------------------------------------------------------
!! differentiate f(:) "order" times and store the result in df(:)
!! using fourier/cos/sin expansion methods
!!------------------------------------------------------------------
!!
!! n      number of grid points
!! order  which order deriv to take
!!        (not the order of convergence)
!! kx     input:   wavenumber vector
!!
!! kfilt  input:  wavenumber filter
!!
!! tmp    input:   workspace 
!!                 real(kind=8) dimension(n)
!!
!! exp_type   'fourier','cos' or 'sin' 
!!
 implicit none
 integer                     :: n,order,i,istart
 real(kind=8)                :: f(n)
 real(kind=8)                :: df(n)
 integer(kind=8)             :: plan_f
 integer(kind=8)             :: plan_i
 integer(kind=8)             :: plans(2)
 character(len=80)           :: exp_type 
 real(kind=8)                :: normfactor,xx
 real(kind=8)                :: kx(n),kfilter(n),tmp(n)
 
 integer,parameter           :: maxorder=12
 integer,save                :: sign_fourier(12)
 integer,save                :: sign_cos(12)
 integer,save                :: sign_sin(12)
 logical,save                :: do_phase_shift(12)
 logical,save                :: first_entry=.TRUE.
 
 if(n<4) stop 'Call to fourier_deriv with vector of length < 4'
 if( first_entry ) then
  !! sign and phase shift params, fourier method  
  sign_fourier(1:4) = (/-1,-1,1,1/)
  sign_fourier(5:8) = (/-1,-1,1,1/)
  sign_fourier(9:12) = (/-1,-1,1,1/)
  do_phase_shift(1:2:maxorder-1)=.TRUE.
  do_phase_shift(2:2:maxorder ) =.FALSE.
  !! sign parameters cos method, input data is assumed even
  sign_cos(1:4) = (/-1,-1,1,1/)     ! cos:  --> -sin --> -cos --> sin --> cos
  sign_cos(5:8) = (/-1,-1,1,1/)
  sign_cos(9:12) = (/-1,-1,1,1/) 
  !! sign parameters sin method, input data is assumed odd
  sign_sin(1:4) = (/1,-1,-1,1/)     ! sin:  --> cos --> -sin --> -cos --> sin
  sign_sin(5:8) = (/1,-1,-1,1/)
  sign_sin(9:12) = (/1,-1,-1,1/)
  first_entry=.FALSE.
 endif
 
 
 if( trim(exp_type) == 'fourier' ) then
 
  normfactor = (1.d0/dfloat(n)) * sign_fourier(order)
  plan_f = plans(1)
  plan_i = plans(2)
  
  !-----------------------------------------------
  ! do the forward transform 
  !-----------------------------------------------
  call dfftw_execute_r2r(plan_f,f,tmp)  !! fourier
 
  !-------------------------------------------------------------
  ! if necessary, the "i" part of -ikx
  ! location shift necessary for real--> 1/2 complex transform
  ! implementation, other ways also possible
  !-------------------------------------------------------------
  if( do_phase_shift(order) ) then  
   do i=2,n/2
    xx = tmp(i)
    tmp(i) = tmp(n-i+2)
    tmp(n-i+2) = xx
   enddo
  endif
  
  !-------------------------------------------------------------
  ! now the "-k_x" and the normalization:
  ! tmp(:) = - kx(:)*normfactor*tmp(:)
  ! ---> generalized for arbitrary orders, see above
  !-------------------------------------------------------------
   do i=1,n
    tmp(i) = tmp(i)*normfactor*(kx(i)*kfilter(i))**order
   enddo
   tmp(n/2+1) = 0.d0    !! zero the nyquist frequency
 
  !-------------------------------------------------------------
  ! do the inverse transform
  !-------------------------------------------------------------
  call dfftw_execute_r2r(plan_i,tmp,df) 
  
 elseif( trim(exp_type) == 'cos' ) then
    
  normfactor = (1.d0/(2.d0*(dfloat(n)-1.d0))) * sign_cos(order)
  plan_f=plans(1)
  if( mod(order,2)==0 ) then   ! inverse is also n point cosine transform
   plan_i=plan_f
   istart=1
  else                         ! inverse is shortened sin transform
   plan_i=plans(2)
   istart=2             
  endif
  
  call dfftw_execute_r2r(plan_f,f,tmp)        !! forward/cos transform
  
  !-----------------------------------------------
  ! wavenumber multiplication 
  !-----------------------------------------------
  do i=1,n
   tmp(i) = tmp(i)*normfactor*(kx(i)*kfilter(i))**order
  enddo
  tmp(n)=0.d0      ! set nyquist to zero
  
  !-------------------------------------------------------------
  ! do the inverse transform
  !-------------------------------------------------------------
  call dfftw_execute_r2r(plan_i,tmp(istart),df(istart))
  
  !-------------------------------------------------------------
  ! if final result is sin expandable, i.e. odd deriv of even fn
  !  ==> zero the endvals
  !-------------------------------------------------------------
  if( mod(order,2) .ne. 0) then
   df(1)=0.d0
   df(n)=0.d0
  endif
    
 elseif( trim(exp_type) == 'sin' ) then
    
  normfactor = 1.d0/(2.d0*(dfloat(n)-1.d0)) * sign_sin(order)
  plan_f=plans(1)
  if( mod(order,2)==0 ) then  ! inverse is also shortened sin transform
   plan_i=plan_f
   istart=2
  else                        ! inverse is n point cosine transform
   plan_i=plans(2)
   istart=1
  endif
  
  !-----------------------------------------------
  !  forward/sin transform
  !-----------------------------------------------
  call dfftw_execute_r2r(plan_f,f(2),tmp(2))
  
  !-----------------------------------------------
  ! wavenumber multiplication 
  !-----------------------------------------------
  do i=1,n
   tmp(i) = tmp(i)*normfactor*(kx(i)*kfilter(i))**order
  enddo
  tmp(n)=0.d0    !! nyquist
  
  !-------------------------------------------------------------
  ! do the inverse transform
  !-------------------------------------------------------------
  call dfftw_execute_r2r(plan_i,tmp(istart),df(istart))
  
  !-------------------------------------------------------------
  ! if final result is sin expandable, i.e. even deriv of odd fn
  !  ==> zero the endvals
  !-------------------------------------------------------------
  if( mod(order,2) == 0) then
   df(1)=0.d0
   df(n)=0.d0
  endif
  
 endif
 
 return 
end subroutine fourier_deriv

subroutine fourier_init(exp_type,n,Lx,kx,kfilter,plan_f,plan_i)
  use mpi_params,            only: myid
  implicit none
  character(len=80)        :: exp_type,this_type
  integer                  :: n,i,k
  integer(kind=8)          :: plan_f,plan_i,plan
     integer(kind=8)          :: plan_cos,plan_sin
  real(kind=8)             :: Lx
  real(kind=8),allocatable :: in(:)
  real(kind=8),allocatable :: out(:)
  real(kind=8)             :: kx(n)
  real(kind=8)             :: kfilter(n)
  real(kind=8)             :: rem,xx
  
  real(kind=8)             :: pi,dk,kstar
  integer                  :: kmax,fft_kind
  character(len=80)        :: reality_in
  character(len=80)        :: reality_out
  character(len=80)        :: label
  integer                  :: rank_trans
  integer                  :: howmany_rank
  integer                  :: n_trans,   n_loop
  integer                  :: is_trans, is_loop
  integer                  :: os_trans, os_loop
  include 'fftw3.f'
  
  if( n < 4 ) then
   if(myid==0) then
    write(0,*) ' ................      Warning: fourier initialization for n < 4: ',  &
    trim( exp_type)
   endif
   kx=0.d0
   kfilter=0.d0
   plan_f=-999
   plan_i=-999
   return
   
  endif
     
     
  if( trim(exp_type)=='fourier' .and. mod(n,2) /= 0 ) then
   write(0,*) 'n = ',n
   stop 'need even number of gridpoints w/ 1 end excluded for fourier implementation'
  elseif( trim(exp_type)=='cos' .and. mod(n,2) /= 1 ) then
   stop 'need odd number of gridpoints w/ both ends included for cos implementation'
  elseif( trim(exp_type)=='sin' .and. mod(n,2) /= 1 ) then
   stop 'need odd number of gridpoints w/ both ends included for sin implementation'
  endif
    
  rank_trans=1      !! 1d transforms
  howmany_rank=1    !! Any loops over additional dimensions 
                    !! specified as though data laid out
                    !! in 1d array in memory (not an issue here)
  
  pi=4.d0*datan(1.d0)
  if( trim(exp_type) .eq. 'sin' .or. trim(exp_type) .eq. 'cos') then
  
    dk = pi/Lx 
    kmax = n
    kstar = .85d0*kmax  !(8./9.)*kmax 
    rem = (kstar-kmax)
    xx = max(4.d0,rem/3.d0)
    do i=1,n              
     kfilter(i) = 1.d0! - 0.5*(1.+tanh((i-kstar)/(xx))) 
     kx(i)=(i-1.)*dk
    enddo
    
    !!=======================================================
    !!make a plan for a sin transform ignoring zero endvalues
    !!=======================================================
    
    !!describe the transforms  in dim 1
    n_trans=n-2              ! length of transform excluding zeros at ends
    is_trans=1               ! stride btwn elements in a single transform
    os_trans=is_trans        ! set output stride = input stride
        
    !!describe the loops over additional dims, i.e. 2,3
    n_loop=1                ! collapse more general scheme to 1d 
    is_loop=1               ! stride betwn transforms
    os_loop=is_loop         ! set output stride = input stride
    
    reality_in=  'real'
    reality_out= 'real'
    this_type='sin'
    call get_kind(reality_in,reality_out,this_type,fft_kind,label)
    
    allocate( in(n),out(n) )
    call dfftw_plan_guru_r2r(plan,              &
                             rank_trans,        &
                             n_trans,           &
                             is_trans,          &
                             os_trans,          &
                             howmany_rank,      &
                             n_loop,            &
                             is_loop,           &
                             os_loop,           &
                             in(2),             &
                             out(2),            &
                             fft_kind,          &
                             FFTW_EXHAUSTIVE)
    if( trim(exp_type) == 'sin' ) plan_f=plan
    if( trim(exp_type) == 'cos' ) plan_i=plan    
       plan_sin = plan_f   ! sin w/ no endpoints, testing only
    
    if(myid==0) then
     !write(0,*) ' ................       testing ',plan,trim(exp_type),n
     call dfftw_execute_r2r(plan,in(2),out(2))
    endif
 
    deallocate(in,out)
    
    
    !!=======================================================
    !!Now make a plan for a cos transform keeping endvalues
    !!=======================================================
    
    !!describe the transforms  in dim 1
    n_trans=n                ! length of transform including zeros at ends
    is_trans=1               ! stride btwn elements in a single transform
    os_trans=is_trans        ! set output stride = input stride
        
    !!describe the loops over additional dims, i.e. 2,3
    n_loop=1                ! collapse more general scheme to 1d 
    is_loop=1               ! stride betwn transforms
    os_loop=is_loop         ! set output stride = input stride
    
    reality_in=  'real'
    reality_out= 'real'
    this_type='cos'
    call get_kind(reality_in,reality_out,this_type,fft_kind,label)
    !write(0,*) 'FOURIER INITIALIZING cos plan label ',trim(label)
    allocate( in(n),out(n) )
    call dfftw_plan_guru_r2r(plan,              &
                             rank_trans,        &
                             n_trans,           &
                             is_trans,          &
                             os_trans,          &
                             howmany_rank,      &
                             n_loop,            &
                             is_loop,           &
                             os_loop,           &
                             in,                &
                             out,               &
                             fft_kind,          &
                             FFTW_EXHAUSTIVE)
    !write(0,*) 'FOURIER INITIALIZING cos plan ',plan
    if( trim(exp_type) == 'cos' ) plan_f=plan
    if( trim(exp_type) == 'sin' ) plan_i=plan    
     plan_cos = plan_f   ! testing only
    
    ! check transform/inverse transform pair for sin
    ! do i=1,n
    !  in(i) = sin(2.*pi*(i-1.)/dfloat(n-1.))
    ! enddo
    
    write(0,*) plan_sin,plan_cos
    if(myid==0) then
     !write(0,*) ' ................       testing ',plan,trim(exp_type),n
     !call dfftw_execute_r2r(plan_sin,in(2),out(2))
     !call dfftw_execute_r2r(plan_cos,out(1),in(1))
    endif
    
     !do i=1,n
     ! write(0,*) i,sin(2.*pi*(i-1.)/dfloat(n-1.)),in(i)/(2.*n)
     !enddo
    
    deallocate(in,out)
     
  elseif( trim(exp_type) == 'fourier' ) then
  
    dk = 2.d0*pi/Lx   
    kmax = n/2+1
    kstar = .85d0*kmax  !(8./9.)*kmax 
    rem = (kstar-kmax)
    xx = max(4.d0,rem/3.d0)
    do k=1,n/2+1       
     kfilter(k) = 1.d0 !- 0.5*(1.+tanh((k-kstar)/(xx)))
     kx(k)=(k-1.)*dk
    enddo
     
    do k=n,n/2+2,-1
     kx(k) = -kx( n-k+2 )     
     kfilter(k) = kfilter(n-k+2)
    enddo
   
    !!describe the transforms  in dim 1
    n_trans=n                ! length of transform
    is_trans=1               ! stride btwn elements in a single transform
    os_trans=is_trans        ! set output stride = input stride
        
    !!describe the loops over additional dims, i.e. 2,3
    n_loop=1                ! collapse more general scheme to 1d 
    is_loop=1               ! stride betwn transforms
    os_loop=is_loop         ! set output stride = input stride
    
    !! construct the forward plan
    reality_in='real'
    reality_out='halfcomplex'
    call get_kind(reality_in,reality_out,exp_type,fft_kind,label)
    allocate( in(n),out(n) )
    call dfftw_plan_guru_r2r(plan_f,            &
                             rank_trans,        &
                             n_trans,           &
                             is_trans,          &
                             os_trans,          &
                             howmany_rank,      &
                             n_loop,            &
                             is_loop,           &
                             os_loop,           &
                             in,                &
                             out,               &
                             fft_kind,          &
                             FFTW_EXHAUSTIVE)
    if(myid==0) then
     !write(0,*) ' ................       testing ',plan_f,trim(exp_type),n
     call dfftw_execute_r2r(plan_f,in,out)
    endif  
    
    !! construct the inverse plan
    reality_in='halfcomplex'
    reality_out='real'
    call get_kind(reality_in,reality_out,exp_type,fft_kind,label)     
    call dfftw_plan_guru_r2r(plan_i,            &
                             rank_trans,        &
                             n_trans,           &
                             is_trans,          &
                             os_trans,          &
                             howmany_rank,      &
                             n_loop,            &
                             is_loop,           &
                             os_loop,           &
                             in,                &
                             out,               &
                             fft_kind,          &
                             FFTW_EXHAUSTIVE)
    if(myid==0) then
     !write(0,*) ' ................       testing ',plan_i,trim(exp_type),n
     call dfftw_execute_r2r(plan_i,in,out)
    endif
    deallocate(in,out)
  endif

 return
end subroutine fourier_init

subroutine get_kind(in,out,type,fft_kind,label)
 implicit none                
 include 'fftw3.f'
 integer,intent(out)               :: fft_kind
 character (len=80),intent(out)    :: label
 character (len=11)                :: in
 character (len=11)                :: out
 character (len=11)                :: type
              
 if( trim(type)=='fourier' ) then
  if( in=='real' .and. out=='complex') then
   fft_kind=FFTW_R2HC
   label='no flag needed'
  elseif( in=='complex' .and. out=='real') then
   fft_kind=FFTW_HC2R
   label='no flag needed'
  elseif( in=='real' .and. out=='halfcomplex') then
   fft_kind=FFTW_R2HC
   label='FFTW_R2HC'
  elseif( in=='halfcomplex' .and. out=='real') then
   fft_kind=FFTW_HC2R
   label='FFTW_HC2R'
  elseif( in=='complex' .and. out=='complex') then
   fft_kind=-999
   label='no flag needed'
  else
   write(0,*) 'problem getting kind of fourier transform',in,out
   stop
  endif        
 endif
       
 if( trim(type)=='sin' ) then
  if( in=='real' .and. out=='real') then
   fft_kind=FFTW_RODFT00
   label='FFTW_RODFT00'
  elseif( in=='complex' .and. out=='complex') then
   fft_kind=-999
   label='no label needed'
  else
   write(0,*) 'problem getting kind of fourier transform',in,out
   stop
  endif  
 endif
 
 if( trim(type)=='cos' .or. trim(type)=='cheby' ) then
  if( in=='real' .and. out=='real') then
   fft_kind=FFTW_REDFT00
   label='FFTW_REDFT00'
  elseif( in=='complex' .and. out=='complex') then
   fft_kind=-999
   label='no label needed'
  else
   write(0,*) 'problem getting kind of fourier transform',in,out
   stop
  endif  
 endif         
end subroutine get_kind



subroutine compact_deriv(f,df,h,n,A,ipiv)
 use mpi_params,          only: myid
 implicit none
 integer                     :: n
 real(kind=8)                :: h
 real(kind=8)                :: f(n)
 real(kind=8)                :: df(n) 
 real(kind=8)                :: A(7,n)
 integer                     :: ipiv(n)
 

!======================================================
!  use parameters of the 6th order scheme described in
!  the Lui & Lele 2001 conference paper 
!  Direct numerical simulation of spatially developing
!  compressible, turbulent mixing layers
!  Calvin Lui and Sanjiva K. Lele 
!  American Institute of Aeronautics and Astronautics
!
! the constants in the paper are given in single precision
! real(kind=8),parameter      :: alpha_1 = 0.5381301 
! real(kind=8),parameter      :: beta_1  = 0.0666332
!
!   regard the scheme as a 2 parameter scheme w/ 
!   alpha and beta given precisely as above, 
!   constrain a,b & c as is required for 6th order scheme
!   ==> more precise values of (a,b,c) than given, consistent
!   and very close to original scheme
!
!   a1 = (1/6)*(9 + alpha - 20*beta)
!   b1 = (1/15)*(-9 +32*alpha + 62*beta)
!   c1 = (1/10)*(1 - 3*alpha + 12*beta)
!======================================================
 real(kind=8),parameter      :: a1 = 1.367577683333333d0 
 real(kind=8),parameter      :: b1 = 0.823428106666667d0
 real(kind=8),parameter      :: c1 = 0.018520810000000d0
 integer,parameter           :: nrhs=1
 integer,parameter           :: kl=2
 integer,parameter           :: ku=2
 integer,parameter           :: seven=7
 
 integer                     :: i,j
 integer                     :: ierr
 real(kind=8),save           :: c(6,7)
 real(kind=8),save           :: hinv
 
 if(n<6) stop 'Call to compact_deriv with vector of length < 6'
  hinv=1.d0/h
  c(1:6,1) = (/ -(197.d0/60.d0), -(5.d0/12.d0), (5.d0), -(5.d0/3.d0), (5.d0/12.d0), -(1.d0/20.d0) /)*hinv
  c(1:5,2) = (/ -(43.d0/96.d0), -(5.d0/6.d0), (9.d0/8.d0), (1.d0/6.d0), -(1.d0/96.d0) /)*hinv
  ! this is just a tridiag, 6th order 'interior' scheme
  c(1:2,3) = (/ (1.d0/36.d0),  (7.d0/9.d0) /)*hinv
  
  c(1:3,4) = (/ (c1/6.d0), (b1/4.d0), (a1/2.d0) /)*hinv
  
  ! this is just a tridiag, 6th order 'interior' scheme
  c(1:2,5) = (/ (1.d0/36.d0),  (7.d0/9.d0) /)*hinv
  ! note, coeffs are just the negatives of those at the left end
  c(1:5,6) = -(/ -(43.d0/96.d0), -(5.d0/6.d0), (9.d0/8.d0), (1.d0/6.d0), -(1.d0/96.d0) /)*hinv
  ! note, coeffs are just the negatives of those at the left end
  c(1:6,7) = -(/ -(197.d0/60.d0), -(5.d0/12.d0), (5.d0), -(5.d0/3.d0), (5.d0/12.d0), -(1.d0/20.d0) /)*hinv
 
 !!------------------------------------------
 !! store the rhs vector in df which is then
 !! overwritten by the derivative
 !!------------------------------------------
 
 
 i=1
 df(i)=0.d0
 do j=1,6
  df(i) = df(i) + c(j,1)*f(j)
 enddo
 
 
 i=2
 df(i)=0.d0
 do j=1,5
  df(i) = df(i) + c(j,2)*f(j)
 enddo
  
 i=3
 df(i)  = ( c(1,3)*(f(i+2)-f(i-2)) + c(2,3)*(f(i+1)-f(i-1)) )  

  
 do i=4,n-3
  df(i)  =     c(1,4)*(f(i+3)-f(i-3))             &
             + c(2,4)*(f(i+2)-f(i-2))             & 
             + c(3,4)*(f(i+1)-f(i-1))    
 enddo
 

 i=n-2
 df(i)  = ( c(1,5)*(f(i+2)-f(i-2)) + c(2,5)*(f(i+1)-f(i-1)) )
 
 i=n-1
 df(i)=0.d0
 do j=1,5
  df(i) = df(i) + c(j,6)*f(n-j+1)
 enddo
  
 i=n
 df(i)=0.d0
 do j=1,6
  df(i) = df(i) + c(j,7)*f(n-j+1)
 enddo
  
 
 !! use the factorization, the rhs and do the solve...
 call dgbtrs('n',n,kl,ku,nrhs,A,seven,ipiv,df,n,ierr)
  
 if(ierr /= 0 ) stop 'dgbtrs problem in compact_deriv'
  
end subroutine compact_deriv


subroutine compact_init_1(n,A,ipiv)
!======================================================
!  use parameters of the 6th order scheme described in
!  the Lui & Lele 2001 conference paper 
!  Direct numerical simulation of spatially developing
!  compressible, turbulent mixing layers
!  Calvin Lui and Sanjiva K. Lele 
!  American Institute of Aeronautics and Astronautics
!
! the constants in the paper are given in single precision
! real(kind=8),parameter      :: alpha_1 = 0.5381301 
! real(kind=8),parameter      :: beta_1  = 0.0666332
!
!   regard the scheme as a 2 parameter scheme w/ 
!   alpha and beta given precisely as above, 
!   constrain a,b & c as is required for 6th order scheme
!   ==> more precise values of (a,b,c) than given, consistent
!   and very close to original scheme
!
!   a1 = (1/6)*(9 + alpha - 20*beta)
!   b1 = (1/15)*(-9 +32*alpha + 62*beta)
!   c1 = (1/10)*(1 - 3*alpha + 12*beta)
!======================================================
  use mpi_params,         only: myid
  implicit none
  real(kind=8),parameter     :: alpha_1 = 0.5381301d0 
  real(kind=8),parameter     :: beta_1  = 0.0666332d0
  integer                    :: n
  integer                    :: i,row
  integer                    :: ierr
  integer, parameter         :: kl=2
  integer, parameter         :: ku=2
  integer, parameter         :: seven=7
  integer                    :: ipiv(n)
  real(kind=8),intent(inout) :: A(seven,n)  
  
  if( n < 6 ) then
   if(myid==0) then
    write(0,*) ' ................      Warning: compact initialization for n < 6: '
    A(:,:)=0.d0
    ipiv(:)=0.d0
    return
   endif
  endif

 A(:,:)=0.
 ! 1st 2 rows of banded matrix array A are reserved for fill-in
 ! 1st/last k elements for the kth super/sub diagonal not accessed

 ! superdiagonal 2
 row=3
 A(row,3:5) = (/ 0.d0, 0.d0, 0.d0/)
 A(row,6:n-1)=beta_1
 A(row,n)=0.d0
   
 ! superdiagonal 1
 row=4
 !A(row,2:4) = (/ 2.d0, 1.d0/2.d0, 1.d0/3.d0 /)   ! 3rd order at nodes 1,n
 !A(row,n-1:n)=(/ 1.d0/3.d0, 1.d0/6.d0 /)         ! 5th order at nodes 2,n-2
 
 A(row,2:4) = (/ 5.d0, 3.d0/4.d0, 1.d0/3.d0 /)    ! 6th order at all nodes
 A(row,5:n-2)=alpha_1 
 A(row,n-1:n)=(/ 1.d0/3.d0, 1.d0/8.d0 /)          ! 6th order at all nodes
   
 ! diagonal 
 row=5
 A(row,1:n)=1.d0
  

! subdiagonal 1
  do i=1,n
   A(6,i)=A(4,n-i+1)
  enddo
! subdiagonal 2
  do i=1,n
   A(7,i)=A(3,n-i+1)
  enddo 
 
 ! factor the banded matrix A using LAPACK routine DGBTRF
 ! A is overwritten with its LU factorization
 call dgbtrf(n,n,kl,ku,A,seven,ipiv,ierr)

 !check error condition
 if(ierr > 0 ) then
  write(0,*) '... exact singularity detected: dgbtrf, compact_init_1 '
  stop
 endif
 if(ierr < 0 ) then
  write(0,*) ierr ,'th input argument to dgbtrf, compact_init_1'
  stop
 endif

 return
end subroutine compact_init_1


subroutine cheby_deriv(f,df,L,n,order,tmp,cplan)
!------------------------------------------------------
! differentiate f(:) and store the result in df(:)
! using chebyshev expansion method
! assume f given for the n gauss lobatto pts
! x(:) in [-1,1]
!------------------------------------------------------
!
! n      number of grid points
!        assume 
!        x in [0,L] on Gauss Lobatto pts
!
! tmp    input:   workspace 
!                 real(kind=8) dimension(n)
!
!------------------------------------------------------
 use mpi_params,          only: myid
 implicit none
 integer                     :: n,i
 integer                     :: order
 integer,parameter           :: inc=1
 real(kind=8)                :: L
 real(kind=8)                :: f(n)
 real(kind=8)                :: df(n)
 real(kind=8)                :: tmp(n)
 real(kind=8)                :: xx
 integer(kind=8)             :: cplan
 
 if(n<4) stop 'Call to cheby_deriv with vector of length < 4'
 !------------------------------------------------------
 ! Stretching & Normalization
 !------------------------------------------------------
 ! the -2 comes about from defining x
 ! in [0,1] rather than from [1,-1] i.e.
 ! x=(1/2)*(1-cos(theta)) rather than x=cos(theta)
 
      xx = (-2.d0/L)**dfloat(order)
 
 ! the L comes from redefining the domain from [0,L]
 ! instead of from [0,1] for consistency of approach
 ! with the other differentiation methods
 !------------------------------------------------------
 
 !------------------------------------------------------
 ! cosine transform the input data
 !------------------------------------------------------
 call dfftw_execute_r2r(cplan,f,df) 
 
  call cheby_recursion_1(n,df,tmp)
 
   do i=2,order
    call dcopy(n,tmp,inc,df,inc)
    call cheby_recursion_1(n,df,tmp)
   enddo
 
 !------------------------------------------------------
 ! inverse transform the result
 !------------------------------------------------------
 call dfftw_execute_r2r(cplan,tmp,df)
 
 
 
 !------------------------------------------------------
 ! do the discrete fft normalization as well 
 !------------------------------------------------------
     xx = xx/(dfloat(n)-1.d0)
 
  call dscal(n,xx,df,inc)
     
 return 
end subroutine cheby_deriv





subroutine cheby_init(n,cplan)
  use mpi_params,       only: myid
  implicit none 
  integer                  :: n
  integer(kind=8)          :: cplan
  real(kind=8),allocatable :: in(:)
  real(kind=8),allocatable :: out(:)
  
  integer                  :: fft_kind
  character(len=80)        :: exp_type
  character(len=80)        :: reality_in
  character(len=80)        :: reality_out
  integer                  :: rank_trans
  integer                  :: howmany_rank
  integer                  :: n_trans,   n_loop
  integer                  :: is_trans, is_loop
  integer                  :: os_trans, os_loop
  include 'fftw3.f'
  
  if( n < 4 ) then
   if(myid==0) then
    write(0,*) ' ................      Warning: cheby initialization for n < 4: '
    cplan=-999
    return
   endif
  endif
  
  if( mod(n,2) /= 1 ) then
   stop 'need odd number of gridpoints w/ both ends included for cheby implementation'
  endif
    
  rank_trans=1      !! 1d transforms
  howmany_rank=1    !! Any loops over additional dimensions 
                    !! specified as though data laid out
                    !! in 1d array in memory (not an issue here)
  
    
  !!=======================================================
  !! Make a plan for a cos transform keeping endvalues
  !!=======================================================
    
  !!describe the transforms  in dim 1
  n_trans=n                ! length of transform including zeros at ends
  is_trans=1               ! stride btwn elements in a single transform
  os_trans=is_trans        ! set output stride = input stride
        
  !!describe the loops over additional dims, i.e. 2,3
  n_loop=1                ! collapse more general scheme to 1d 
  is_loop=1               ! stride betwn transforms
  os_loop=is_loop         ! set output stride = input stride
    
  reality_in  = 'real'
  reality_out = 'real'
  exp_type    = 'cos'
  fft_kind    = FFTW_REDFT00
  allocate( in(n),out(n) )
  call dfftw_plan_guru_r2r(cplan,             &
                           rank_trans,        &
                           n_trans,           &
                           is_trans,          &
                           os_trans,          &
                           howmany_rank,      &
                           n_loop,            &
                           is_loop,           &
                           os_loop,           &
                           in,                &
                           out,               &
                           fft_kind,          &
                           FFTW_EXHAUSTIVE)
  deallocate(in,out)
   
 return
end subroutine cheby_init


subroutine cheby_recursion_1(n,a,b)
  implicit none  
  integer       :: k,n
  real(kind=8)  :: a(n)    ! cheby coeffs for function f
  real(kind=8)  :: b(n)    ! cheby coeffs for function f' 
    
  
  b(n) = a(n)*.5d0
  b(n-1) = (dfloat(n)-1.d0)*a(n)

  do k=n-3,1,-1   
   b(k+1) = b(k+3) + (dfloat(k)+1.d0)*a(k+2)
  enddo
  
  b(1) = a(2) + b(3)
   
 end subroutine cheby_recursion_1
 
 
 
 
 subroutine fd2_deriv(f,df,h,n)
!!===================================================
!! differentiate f(:), store result in df(:)
!! n      number of grid points
!! h      uniform grid spacing

!!===================================================
!! xlf90 -qfree -qsuffix=f=f90 -c fd2.f90
!! gfortran -c fd2.f90
!!===================================================
 implicit none
 integer                     :: n
 real(kind=8)                :: f(n),df(n),h
  
 integer                     :: i
 real(kind=8)                :: xx
 
 xx=1.d0/(2.d0*h)
 do i=2,n-1
  df(i) = ( f(i+1)-f(i-1) )*xx
 enddo
 df(1) = ( -3.d0*f(1)+4.d0*f(2)-f(3) )*xx
 df(n) = ( f(n-2)-4.d0*f(n-1)+3.d0*f(n) )*xx
 
 
  
end subroutine fd2_deriv
!!===================================================
!!===================================================

