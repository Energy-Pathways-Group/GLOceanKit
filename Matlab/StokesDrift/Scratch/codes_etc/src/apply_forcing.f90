subroutine apply_forcing
 use mpi_params,                      only: myid
 use decomposition_params
 use dimensional_scales,              coriolis=>f
 use independent_variables,           only: x,y,z,t_secs
 use intermediate_variables
 use methods_params
 use etc,                             only: MM0
 
 implicit none
 integer                                      :: id,i0,i1
 integer,parameter                            :: inc=1
 real(kind=8),save                            :: xx(5)
 real(kind=8),save                            :: inv_length_scale
 
 real(kind=8),dimension(:),pointer,save       :: xvals
 real(kind=8),dimension(:),pointer,save       :: yvals
 real(kind=8),dimension(:),pointer,save       :: zvals
 real(kind=8),dimension(:,:,:),pointer,save   :: F
 real(kind=8),dimension(:,:,:),pointer        :: rhs
 integer,save                                 :: npts(3),npvs
 logical,save                                 :: first_entry=.TRUE.
  
 if( do_forcing ) then
 
 if( first_entry ) then
  npts(1) = array_size(JDIM,YBLOCK,myid)  ! # of x pts in YBLOCK
  npts(2) = array_size(IDIM,YBLOCK,myid)  ! # of y pts in YBLOCK
  npts(3) = array_size(KDIM,YBLOCK,myid)  ! # of z pts in YBLOCK
  inv_length_scale = 1.d0/length_scale
  xx(1:3) = 1./(velocity_scale/time_scale)
  xx(4)  =  1./(scalar_scale(1)/time_scale)
  if( do_second_scalar ) then
   xx(5)  =  1./(scalar_scale(2)/time_scale)
   npvs=5
  else
   npvs=4
  endif
  
  !------------------------------------
  ! set some pointers that don't change
  !------------------------------------
  F => tmpY(:,:,:,1)
  i0 = global_x_indices(START,YBLOCK,myid)
  i1 = global_x_indices(END,YBLOCK,myid)
  xvals => x(i0:i1)  

  i0 = global_y_indices(START,YBLOCK,myid)
  i1 = global_y_indices(END,YBLOCK,myid)
  yvals => y(i0:i1) 

  i0 = global_z_indices(START,YBLOCK,myid)
  i1 = global_z_indices(END,YBLOCK,myid) 
  zvals => z(i0:i1) 
  first_entry=.FALSE.
 endif
 
 
 !------------------------------------------
 ! user gets dimensional coord values
 ! so scale the dimensionless values
 !------------------------------------------
  call dscal(npts(1),length_scale,xvals,inc)
  call dscal(npts(2),length_scale,yvals,inc)
  call dscal(npts(3),length_scale,zvals,inc)
 
 !----------------------------------------- 
 ! First pass, assume each field is forced 
 ! but then let user have finer control over
 ! which variables need to be forced.
 ! time independent forcing fields are not
 ! saved, calls are needed each time step
 !-----------------------------------------
 do id=1,npvs
 
  rhs => explicit_rhs(:,:,:,id,MM0)
  
  if( forcing_key(id) ) then
   
   call user_forcing                                  &
         (xvals,yvals,zvals,t_secs,F,id,              &
          npts(1),npts(2),npts(3),forcing_key(id)) 
   
   !------------------------------------
   ! nondimensionalize returned values
   ! add to previously computed values
   !------------------------------------
   call axpy(xx(id),F,rhs) 
   
  endif
 enddo  
 
 !----------------------------------------------
 ! Generally, code uses non-dimensional 
 ! values ==> re-scale
 !----------------------------------------------
  call dscal(npts(1),inv_length_scale,xvals,inc)
  call dscal(npts(2),inv_length_scale,yvals,inc)
  call dscal(npts(3),inv_length_scale,zvals,inc)
 endif  ! end if do_forcing block
   
 return
end subroutine apply_forcing
