subroutine SetupDerivs
 use differentiation_params
 use etc,                           only: logfile
 use mpi_params,                    only: myid,comm,ierr
 use independent_variables,         only: nx,ny,nz,Lx,Ly,Lz,x,y,z,stretch
 use methods_params,                only: deriv_type
 use dimensional_scales,            only: length_scale
 implicit none 
 
 real(kind=8),dimension(:),pointer     :: k
 real(kind=8),dimension(:),pointer     :: kf
 
 integer                               :: n
 integer                               :: order
 real(kind=8)                          :: L,pi 
 integer                               :: idim
 integer                               :: id,i
 real(kind=8),allocatable              :: test_fn(:)
 real(kind=8),allocatable              :: test_deriv(:) 
 real(kind=8),allocatable              :: wrk(:)
 integer(kind=8)                       :: dummy_plan
 character(len=80)                     :: exp_type,method
 logical                               :: DO_DERIV_TESTS=.FALSE.

 
 
 if(myid==0) then
  write(0,*) ' ................'
  write(0,*) ' ................     hello world from SetupDerivs  (a bit slow, fftw3 exhaustive) '
  open(1,file=logfile,position='append')
  write(1,*) '  '
  write(1,*) '  '
  write(1,*) ' =========================================================== '
  write(1,*) ' =========================================================== '
  write(1,*) '                  SetupDerivs Report:' 
  write(1,*) '============================================================ '
  write(1,*) '============================================================ '
  write(1,*) '  '
 endif
 
 
 ! Allocate all the differentiation arrays up front
 allocate( kx(nx),ky(ny),kz(nz) )
 allocate(kxfilter(nx), kyfilter(ny), kzfilter(nz) )
  
 do idim=1,3
 
  ! point at the appropriate arrays, depending on the dimension
  if( idim==1 ) then
   n=nx
   L=Lx/length_scale  ! this used to be 1.0 and the Lx/length_scale was done with chain rule
   k => kx
   kf => kxfilter
  elseif( idim==2 ) then
   n=ny
   L=Ly/length_scale   ! this used to be 1.0 and the Lx/length_scale was done with chain rule
   k => ky
   kf => kyfilter
  elseif( idim==3 ) then
   n=nz
   L=1.d0  ! Lz/Lz
   k => kz
   kf => kzfilter
  endif
 
 ! loop over the variables and look at differentiation method,
 ! initializing the method if necessary
 order=1
 do id=1,6
  if(myid==0) write(0,*) ' ................      setting up derivs for id,idim ',id,idim
  if( trim(deriv_type(id,idim,order)) == 'fourier'   &
      .and. .not. fourier_done(idim) ) then
      
   exp_type = 'fourier'
   call fourier_init(exp_type,n,L,k,kf,              &
                     fourier_plan(idim,1),           &
                     fourier_plan(idim,2)   )
   fourier_done(idim) = .TRUE.
                     
  elseif( trim(deriv_type(id,idim,order)) == 'sin'    &
      .and. .not. sin_done(idim) ) then 
   
   exp_type = 'sin'
   call fourier_init(exp_type,n,L,k,kf,               &
                     sin_plan(idim),dummy_plan )
   sin_done(idim) = .TRUE.
   
   exp_type = 'cos'
   call fourier_init(exp_type,n,L,k,kf,               &
                     cos_plan(idim),dummy_plan )
   cos_done(idim) = .TRUE.

                         
  elseif( trim(deriv_type(id,idim,order)) == 'cos'    &
      .and. .not. cos_done(idim) ) then
        
   exp_type = 'cos'
   call fourier_init(exp_type,n,L,k,kf,               &
                     cos_plan(idim),dummy_plan )
   cos_done(idim) = .TRUE.
   
   exp_type = 'sin'
   call fourier_init(exp_type,n,L,k,kf,               &
                     sin_plan(idim),dummy_plan )
   sin_done(idim) = .TRUE.
                            
  endif
 enddo  ! end loop over variables
 
enddo  ! end loop over spatial dimensions 
 
 
 if( .not. fourier_done(1)  .and.  &
     .not. cos_done(1)      .and.  &
     .not. sin_done(1) ) then
      deallocate( kx,kxfilter )
      allocate( kx(1),kxfilter(1) )
 endif
     
 if( .not. fourier_done(2)  .and.  &
     .not. cos_done(2)      .and.  &
     .not. sin_done(2) ) then
      deallocate( ky,kyfilter )
      allocate( ky(1),kyfilter(1) )
 endif
     
 if( .not. fourier_done(3)  .and.  &
     .not. cos_done(3)      .and.  &
     .not. sin_done(3) ) then
      deallocate( kz,kzfilter )
      allocate( kz(1),kzfilter(1) )
 endif
 
 if(myid==0) then
  do i=0,1
   write(i,*) ' ................      fourier setup in x ',fourier_done(1),fourier_plan(1,:)
   write(i,*) ' ................      fourier setup in y ',fourier_done(2),fourier_plan(2,:)
   write(i,*) ' ................      fourier setup in z ',fourier_done(3),fourier_plan(3,:)
   write(i,*) ' ................      sin setup in x     ',sin_done(1),sin_plan(1)
   write(i,*) ' ................      sin setup in y     ',sin_done(2),sin_plan(2)
   write(i,*) ' ................      sin setup in z     ',sin_done(3),sin_plan(3)
   write(i,*) ' ................      cos setup in x     ',cos_done(1),cos_plan(1)
   write(i,*) ' ................      cos setup in y     ',cos_done(2),cos_plan(2)
   write(i,*) ' ................      cos setup in z     ',cos_done(3),cos_plan(3)
   write(i,*) ' ................'
  enddo
  write(1,*) ' -----> SetupDerivs routine exiting normally  <---------- '
  close(1)
 endif
 

if( .NOT. DO_DERIV_TESTS ) return
 
!----------------------------------------------------------
! test cos differentiation in z direction
!----------------------------------------------------------
 if( myid==0 .and. cos_done(3) ) then
  pi=4.d0*datan(1.d0) 
  allocate( test_fn(nz),test_deriv(nz),wrk(nz) )

  do i=1,nz
   test_fn(i)=cos(pi*z(i))
  enddo

  method='cos'
  idim=3
  order=1
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,nz
   test_fn(i)=-pi*sin(pi*z(i))
  enddo

  open(1,file='output/cos_test_z')
  do i=1,nz
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') z(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
  
  if( .not. stretch(3) ) then
  !reset test function and test 2nd derivative
  do i=1,nz
   test_fn(i)=cos(pi*z(i))
  enddo

  method='cos'
  idim=3
  order=2
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,nz
   test_fn(i)=-pi*pi*cos(pi*z(i))
  enddo

  open(1,file='output/cos_test_zz')
  do i=1,nz
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') z(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
 endif
 deallocate( test_fn,test_deriv,wrk )
 endif
 
!----------------------------------------------------------
! test sin differentiation in z direction
!----------------------------------------------------------
 if( myid==0 .and. sin_done(3) ) then
  pi=4.d0*datan(1.d0) 
  allocate( test_fn(nz),test_deriv(nz),wrk(nz) )

  do i=1,nz
   test_fn(i)=sin(pi*z(i))
  enddo

  method='sin'
  idim=3
  order=1
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,nz
   test_fn(i)=pi*cos(pi*z(i))
  enddo

  open(1,file='output/sin_test_z')
  do i=1,nz
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') z(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
  
  if( .not. stretch(3) ) then
  !reset test function and test 2nd derivative
  do i=1,nz
   test_fn(i)=sin(pi*z(i))
  enddo

  method='sin'
  idim=3
  order=2
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,nz
   test_fn(i)=-pi*pi*sin(pi*z(i))
  enddo

  open(1,file='output/sin_test_zz')
  do i=1,nz
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') z(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1) 
 endif
 deallocate( test_fn,test_deriv,wrk )
 endif
 
!----------------------------------------------------------
! test fourier differentiation in z direction
!----------------------------------------------------------
 if( myid==0 .and. fourier_done(3) .and. nz>4 ) then
  pi=4.d0*datan(1.d0) 
  allocate( test_fn(nz),test_deriv(nz),wrk(nz) )

  do i=1,nz
   test_fn(i) = cos(pi*z(i)) + sin(pi*z(i))
  enddo

  method='fourier'
  idim=3
  order=1
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,nz
   test_fn(i)=-pi*sin(pi*z(i))
  enddo

  open(1,file='output/fourier_test_z')
  do i=1,nz
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') z(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
  
  if( .not. stretch(3) ) then
  !reset test function and test 2nd derivative
  do i=1,nz
   test_fn(i) = cos(pi*z(i)) + sin(pi*z(i))
  enddo

  method='fourier'
  idim=3
  order=2
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,nz
   test_fn(i)=-pi*pi*( cos(pi*z(i)) + sin(pi*z(i)) )
  enddo

  open(1,file='output/fourier_test_z')
  do i=1,nz
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') z(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
 endif
 deallocate( test_fn,test_deriv,wrk )
 endif
 
 
 
 
 
 
!----------------------------------------------------------
! test cos differentiation in x direction
!----------------------------------------------------------
 if( myid==0 .and. cos_done(1) ) then
  pi=4.d0*datan(1.d0)
  allocate( test_fn(nx),test_deriv(nx),wrk(nx) )

  do i=1,nx
   test_fn(i)=cos(pi*x(i))
  enddo

  method='cos'
  idim=1
  order=1
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  ! exact 1st derivative
  do i=1,nx
   test_fn(i)=-(pi)*sin(pi*x(i))
  enddo

  open(1,file='output/cos_test_x')
  do i=1,nx
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') x(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
  
  if( .not. stretch(1) ) then
  ! reset test function and take 2nd derivative
  do i=1,nx
   test_fn(i)=cos(pi*x(i))
  enddo
  method='cos'   !====> method always refers to 1st deriv method
  idim=1
  order=2
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  ! exact 2nd derivative
  do i=1,nx
   test_fn(i)=-(pi**2)*cos(pi*x(i))
  enddo

  open(1,file='output/cos_test_xx')
  do i=1,nx
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') x(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
 endif
 deallocate( test_fn,test_deriv,wrk )
 endif
 
!----------------------------------------------------------
! test sin differentiation in x direction
!----------------------------------------------------------
 if( myid==0 .and. sin_done(1) ) then
  pi=4.d0*datan(1.d0) 
  allocate( test_fn(nx),test_deriv(nx),wrk(nx) )

  do i=1,nx
   test_fn(i)=sin(pi*x(i))
  enddo

  method='sin'
  idim=1
  order=1
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  ! exact 1st derivative
  do i=1,nx
   test_fn(i)= (pi)*cos(pi*x(i))
  enddo

  open(1,file='output/sin_test_x')
  do i=1,nx
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') x(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
  
  if( .not. stretch(1) ) then
  ! reset test function and take 2nd derivative
  do i=1,nx
   test_fn(i)=sin(pi*x(i))
  enddo
  method='sin'   !====> method always refers to 1st deriv method
  idim=1
  order=2
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  ! exact 2nd derivative
  do i=1,nx
   test_fn(i)=-(pi**2)*sin(pi*x(i))
  enddo

  open(1,file='output/sin_test_xx')
  do i=1,nx
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') x(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
 endif
 deallocate( test_fn,test_deriv,wrk )
 endif
 
!----------------------------------------------------------
! test fourier differentiation in x direction
!----------------------------------------------------------
 if( myid==0 .and. fourier_done(1) .and. nx>4 ) then
  pi=4.d0*datan(1.d0) 
  allocate( test_fn(nx),test_deriv(nx),wrk(nx) )

  do i=1,nx
   test_fn(i)=cos(pi*x(i)) + sin(pi*x(i))
  enddo

  method='fourier'
  idim=1
  order=1
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  ! exact 1st derivative
  do i=1,nx
   test_fn(i) = -pi*sin(pi*x(i))  + pi*cos(pi*x(i))
  enddo

  open(1,file='output/fourier_test_x')
  do i=1,nx
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') x(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
  
  if( .not. stretch(1) ) then
  ! reset test function and take 2nd derivative
  do i=1,nx
   test_fn(i)=cos(pi*x(i)) + sin(pi*x(i))
  enddo
  method='fourier'
  idim=1
  order=2
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  ! exact 2nd derivative
  do i=1,nx
   test_fn(i) = -pi**2*( cos(pi*x(i))  + sin(pi*x(i)) )
  enddo

  open(1,file='output/fourier_test_xx')
  do i=1,nx
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') x(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
 endif
 deallocate( test_fn,test_deriv,wrk )
 endif
 
 
 
 
 
 
 
 
 
!----------------------------------------------------------
! test cos differentiation in y direction
!----------------------------------------------------------
 if( myid==0 .and. cos_done(2) ) then
  pi=4.d0*datan(1.d0) 
  allocate( test_fn(ny),test_deriv(ny),wrk(ny) )

  do i=1,ny
   test_fn(i)=cos(pi*y(i))
  enddo

  method='cos'
  idim=2
  order=1
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,ny
   test_fn(i)=-pi*sin(pi*y(i))
  enddo

  open(1,file='output/cos_test_y')
  do i=1,ny
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') y(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
  
  if( .not. stretch(2) ) then
  ! reset test function and take 2nd derivative
  do i=1,ny
   test_fn(i)=cos(pi*y(i))
  enddo

  method='cos'
  idim=2
  order=2
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,ny
   test_fn(i)=-(pi**2)*cos(pi*y(i))
  enddo

  open(1,file='output/cos_test_yy')
  do i=1,ny
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') y(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
 endif
 deallocate( test_fn,test_deriv,wrk )
 endif
 
!----------------------------------------------------------
! test sin differentiation in y direction
!----------------------------------------------------------
 if( myid==0 .and. sin_done(2) ) then
  pi=4.d0*datan(1.d0) 
  allocate( test_fn(ny),test_deriv(ny),wrk(ny) )

  do i=1,ny
   test_fn(i)=sin(pi*y(i))
   test_deriv(i)=0.d0
   wrk(i)=0.d0
  enddo
  
   !write(0,*) 'CHECKING SIN & INV TRANSFORM ', sin_plan(2),cos_plan(2)
   !call dfftw_execute(sin_plan(2),test_fn(2),test_deriv(2))
   !call dfftw_execute(cos_plan(2),test_deriv,wrk)
   !do i=1,ny
   ! write(10,*) i,test_fn(i),wrk(i)
   !enddo   
    

  method='sin'
  idim=2
  order=1
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,ny
   test_fn(i)=pi*cos(pi*y(i))
  enddo

  open(1,file='output/sin_test_y')
  do i=1,ny
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') y(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
  
  if( .not. stretch(2) ) then
  ! reset test function and take 2nd derivative
  do i=1,ny
   test_fn(i)=sin(pi*y(i))
  enddo

  method='sin'
  idim=2
  order=2
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,ny
   test_fn(i)=-(pi**2)*sin(pi*y(i))
  enddo

  open(1,file='output/sin_test_yy')
  do i=1,ny
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') y(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
 endif
 deallocate( test_fn,test_deriv,wrk )
 endif
 
!----------------------------------------------------------
! test fourier differentiation in y direction
!----------------------------------------------------------
 if( myid==0 .and. fourier_done(2) .and. ny>4 ) then
  pi=4.d0*datan(1.d0)
  allocate( test_fn(ny),test_deriv(ny),wrk(ny) )

  do i=1,ny
   test_fn(i)=cos(pi*y(i)) + sin(pi*y(i))
  enddo

  method='fourier'
  idim=2
  order=1
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,ny
   test_fn(i) = pi*(-sin(pi*y(i)) + cos(pi*y(i)) )
  enddo

  open(1,file='output/fourier_test_y')
  do i=1,ny
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') y(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
  
  if( .not. stretch(2) ) then
  ! reset test function and take 2nd derivative
  do i=1,ny
   test_fn(i)=cos(pi*y(i)) + sin(pi*y(i))
  enddo

  method='fourier'
  idim=2
  order=2
  call differentiate(test_fn,test_deriv,idim,method,wrk,order)

  do i=1,ny
   test_fn(i) = -pi*pi*(cos(pi*y(i)) + sin(pi*y(i)) )
  enddo

  open(1,file='output/fourier_test_yy')
  do i=1,ny
   wrk(i) = test_fn(i) - test_deriv(i)
   write(1,'(4g14.6)') y(i),test_fn(i),test_deriv(i),wrk(i)
  enddo
  close(1)
 endif
 deallocate( test_fn,test_deriv,wrk )
 endif
 
  
 return
end subroutine SetupDerivs
