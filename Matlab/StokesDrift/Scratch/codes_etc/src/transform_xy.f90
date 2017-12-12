subroutine transform_xy(in,out,dir,exp_type)
!------------------------------------------------------
! Transform ENTIRE 3d YBLOCK data in x and y
! dir indicates forward or inverse transforms
! in and out can be the same array 
!------------------------------------------------------
 use mpi_params,             only: myid
 use decomposition_params
 use differentiation_params
 use etc
 use independent_variables,  only: nx,ny 
 use intermediate_variables, only: tmpX,tmpY
 
 implicit none     
 real(kind=8)                  ::  in( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid) )
 
 real(kind=8)                  :: out( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid) )
  
 real(kind=8),save             :: xfac,yfac,normfac
 integer                       :: i,j,k,dim,dim_start
 integer                       :: dir,mydir,i0
 character(len=80)             :: exp_type(2)
 logical,save                  :: first_entry=.TRUE.
 
 
  if( nx==1 ) then
   dim_start=2      ! don't need x transform plans if nx=1
  else
   dim_start=1
  endif
        
  do dim=dim_start,2  !! x & y dimensions 
   ! re-use the sin plan from deriv setup
   if( trim(exp_type(dim))=='sin' ) then
    if( sin_done(dim) ) then
     xy_plan(dim,1)=sin_plan(dim)  ! for in dim   n-2 pt sin transform
     xy_plan(dim,2)=sin_plan(dim)  ! inv in dim   ""
    else
     write(0,*) dim,sin_done
     stop 'need sin plan in transform_xy'
    endif
  
   ! re-use the cos plan from deriv setup
   elseif( trim(exp_type(dim))=='cos' ) then
    if( cos_done(dim) ) then
     xy_plan(dim,1)=cos_plan(dim)  ! for in dim   n pt cos transform
     xy_plan(dim,2)=cos_plan(dim)  ! inv in dim   ""
    else
     write(0,*) dim
     stop 'need cos plan in transform_xy'
    endif
  
   ! re-use the fourier plan from deriv setup
   elseif( trim(exp_type(dim))=='fourier' ) then
    if( fourier_done(dim) ) then
     xy_plan(dim,1)=fourier_plan(dim,1)  ! for in dim
     xy_plan(dim,2)=fourier_plan(dim,2)  ! inv in dim
    else
     write(0,*) dim
     stop 'need fourier plan in transform_xy'
    endif
   endif  
  enddo
 
  !---------------------
  ! x normalization
  !---------------------
  if( nx > 1 ) then 
   if( trim(exp_type(1)) == 'cos' ) then
    xfac = 1.d0/(2.d0*(nx-1.d0))
   elseif( trim(exp_type(1)) == 'sin' ) then
    xfac = 1.d0/(2.d0*(nx-1.d0))
   elseif( trim(exp_type(1)) == 'fourier' ) then
    xfac = 1.d0/dfloat(nx)
   endif   
  elseif(nx==1) then
   xfac=1.d0
  endif
 
  !---------------------
  ! y normalization
  !---------------------
  if( ny > 1 ) then   
   if( trim(exp_type(2)) == 'cos' ) then
    yfac = 1.d0/(2.d0*(ny-1.d0))
   elseif( trim(exp_type(2)) == 'sin' ) then
    yfac = 1.d0/(2.d0*(ny-1.d0))
   elseif( trim(exp_type(2)) == 'fourier' ) then
    yfac = 1.d0/dfloat(ny)
   endif   
  elseif( ny == 1 ) then
   yfac=1.d0
  endif
  normfac = xfac*yfac
 
 !! array index corresponding to forward/inverse transform direction
 if(dir==1) then
  mydir=1
 elseif(dir==-1) then
  mydir=2
 else
  stop ' dir is wrong in xy transform'
 endif
 
 
 !-------------------------------------------------------
 ! do the x transforms
 !-------------------------------------------------------
 if( trim(exp_type(1)) == 'sin' ) then
  i0=2    ! exclude data endpoints for sin transform
 else
  i0=1
 endif
 
 if( nx > 1 ) then   
  !----------------------------------------
  !  transpose to XBLOCK format
  !----------------------------------------
  call yblock_2_xblock(in,tmpX)
  
  !------------------------------
  !  loop and x transform
  !------------------------------
 
   do k=1,array_size(KDIM,XBLOCK,myid)  
    do j=1,array_size(JDIM,XBLOCK,myid)
     call dfftw_execute_r2r(xy_plan(1,mydir),   &
                            tmpX(i0,j,k,1),      &
                            tmpX(i0,j,k,2))
    enddo
   enddo
   
  if( trim(exp_type(1)) == 'sin' ) then
   tmpX(1,:,:,2)=0.d0    ! "mean"/endpoint location in extended array
   tmpX(nx,:,:,2)=0.d0   ! "nyquist"/endpoint location in extended array
  endif


  !----------------------------------------
  !  transpose back to YBLOCK format
  !  store result in tmpY(:,:,:,1)
  !----------------------------------------
  if( ny > 1 ) then      ! result goes in tmpY
   call xblock_2_yblock(tmpX(1,1,1,2),tmpY(1,1,1,1))
  elseif( ny == 1 ) then  ! result goes directly into 'out'
   call xblock_2_yblock(tmpX(1,1,1,2),out)
   goto 999    ! don't do y transforms
  endif
  
 endif
 

 !-------------------------------------------
 ! y transform the results in tmpY(:,:,:,1)
 ! and store the results in out
 !-------------------------------------------
 if( trim(exp_type(2)) == 'sin'  ) then
  i0=2    ! exclude data endpoints for sin transform
 else
  i0=1
 endif
   
  if( nx > 1 ) then       ! input is in tmpY 

      do k=1,array_size(KDIM,YBLOCK,myid)  
       do j=1,array_size(JDIM,YBLOCK,myid)
        call dfftw_execute_r2r(xy_plan(2,mydir),   &
                               tmpY(i0,j,k,1),      &
                                out(i0,j,k) )              
       enddo
      enddo
      


  elseif( nx==1 ) then   ! input still in 'in'
 

      do k=1,array_size(KDIM,YBLOCK,myid)  
       do j=1,array_size(JDIM,YBLOCK,myid)
        call dfftw_execute_r2r(xy_plan(2,mydir),   &
                               in(i0,j,k),          &
                               out(i0,j,k) )

       enddo 
      enddo
            
      
      
  endif
  
 if( trim(exp_type(2)) == 'sin' ) then
  out(1,:,:)=0.d0    ! "mean"/endpoint location in extended array
  out(ny,:,:)=0.d0   ! "nyquist"/endpoint location in extended array
 endif

 

 
999 continue
if( dir==-1 ) then
 !-------------------------------------------
 ! only normalize the inverse transforms
 !-------------------------------------------
 call scale(normfac,out)

endif 

end subroutine transform_xy




subroutine test_transform_xy(in,out)
!------------------------------------------------------
! Transform ENTIRE 3d YBLOCK data in x and y
! dir indicates forward or inverse transforms
! in and out can be the same array 
!------------------------------------------------------
 use mpi_params,             only: myid
 use decomposition_params
 use intermediate_variables, only: tmpY
 use dimensional_scales,     only: length_scale
 use independent_variables,  only: y,Ly,ny
 use methods_params,         only: deriv_type
 
 implicit none     
 real(kind=8)                  ::  in( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid) )
 
 real(kind=8)                  :: out( array_size(IDIM,YBLOCK,myid),   &
                                       array_size(JDIM,YBLOCK,myid),   &
                                       array_size(KDIM,YBLOCK,myid) )
  
 integer                       :: dir,i,j,k,dim,dim_start
 character(len=80)             :: exp_type(2)
 real(kind=8)                  :: pi,ky
 
 
 !   test the 2d transform of a function that is sin expanded in y and (1pt) fourier expanded in x
 pi=4.*datan(1.d0)
 ky=length_scale*1*pi/Ly   ! dless wavenumber
 
 do k=1,array_size(KDIM,YBLOCK,myid)
  do i=1,array_size(JDIM,YBLOCK,myid)
   do j=1,array_size(IDIM,YBLOCK,myid)
   
    if( trim(deriv_type(2,2,1))=='sin' ) then   ! ! v, y dir, 1st deriv
     in(j,i,k) = sin(ky*y(j))
    elseif( trim(deriv_type(2,2,1))=='cos' ) then
     in(j,i,k) = cos(ky*y(j))
    else
     in(j,i,k) = .75*cos(2*ky*y(j)) + .25*sin(2*ky*y(j))
    endif
    
   enddo
  enddo
 enddo
 
 
  dir = 1   ! 1=for  -1=inv
  
  exp_type(1)=deriv_type(2,1,1)   ! v, x dir, 1st deriv
  exp_type(2)=deriv_type(2,2,1)   ! v, y dir, 1st deriv
  call transform_xy(in,out,dir,exp_type)
  
  dir=-1
  call transform_xy(out,in,dir,exp_type)
 
 k=(array_size(KDIM,YBLOCK,myid)-1)/2   ! pick a "random" level to test
 i=1
 do k=1,array_size(KDIM,YBLOCK,myid)
 do j=1,array_size(IDIM,YBLOCK,myid)
  if( abs( in(j,i,k) - sin(ky*y(j)) ) .gt. 1.e-15 ) then
   write(0,*) 'problem in transform_xy for v ',j,i,k,myid
   stop
  endif
 enddo
 enddo
 
 if(myid==0) write(0,*) '....... test of transform_xy looks fine '
 
 return
end subroutine test_transform_xy


