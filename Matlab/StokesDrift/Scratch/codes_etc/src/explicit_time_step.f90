subroutine explicit_step
!!  Routine to take an explicit time step via Euler,
!!  2nd or third order Adams Bashforth method. The differential
!!  equations are as follows:
!!
!!  d/dt u   = rhs1
!!  d/dt v   = rhs2
!!  d/dt w   = rhs3
!!  d/dt s1  = rhs4
!!  d/dt s2  = rhs5
!!
!!  notation:
!!  dt   time increment
!!  MM0   integer index identifying the time t    storage locations
!!  MM1   integer index identifying the time t-dt storage locations
!!  MM2   integer index identifying the time t-2*dt storage locations
!!  MM3   integer index identifying the time t-3*dt storage locations
!!
!!  explicit_rhs(1,MM0),
!!  explicit_rhs(1,MM1),
!!  explicit_rhs(1,MM2),
!!  explicit_rhs(1,MM3)  
!!       ====>   rhs for u at t,t-dt,t-2*dt,,t-3*dt
!!  etc for explicit_rhs(2,:),
!!          explicit_rhs(3,:),
!!          explicit_rhs(4,:) 
!!          explicit_rhs(5,:)

!!
!!  outputs:
!!  calculated values of u at time t+dt
!!  calculated values of v at time t+dt
!!  calculated values of w at time t+dt
!!  calculated values of scalar(1) at time t+dt (overwrite)
!!  calculated values of scalar(2) at time t+dt (overwrite)
!!=======================================================================
!!=======================================================================
 use pde_params,              only: Re,Pr
 use methods_params,          only: AB_order,do_second_scalar
 use decomposition_params
 use independent_variables,   only: dt,nx
 use intermediate_variables,  only: explicit_rhs
 use diffusion_params
 use dependent_variables
 use pde_params,              only: Rot
 
 use etc
 use mpi_params 
 implicit none

 integer                :: i,j,k,k_global
 real(kind=8), save     :: dt2,dt12,dt24
 real(kind=8)           :: xx
 logical, save          :: first_entry = .TRUE.

 if( first_entry ) then
  dt2 = dt/2.d0
  dt12= dt/12.d0
  dt24= dt/24.d0
  first_entry = .FALSE.
 endif
  
  
 if( step_flag == 'euler') then
 
  if(myid==0) write(0,*) '... 1st order Euler step '

  if( nx>1 .or. Rot(1) .ne. 0.d0 ) then
   u(:,:,:) = u(:,:,:)  +  dt*explicit_rhs(:,:,:,1,MM0)
  endif
  v(:,:,:) = v(:,:,:)  +  dt*explicit_rhs(:,:,:,2,MM0)
  w(:,:,:) = w(:,:,:)  +  dt*explicit_rhs(:,:,:,3,MM0)
  s1(:,:,:) = s1(:,:,:) + dt*explicit_rhs(:,:,:,4,MM0)
  if( do_second_scalar ) then
   s2(:,:,:) = s2(:,:,:) + dt*explicit_rhs(:,:,:,5,MM0)
  endif
 

  if(AB_order>1) step_flag='AB2'
  
 elseif( step_flag == 'AB2' ) then
  
  if(myid==0) write(0,*) '... 2nd order AB step '

  if( nx>1 .or. Rot(1) .ne. 0.d0 ) then
   u(:,:,:)  =  u(:,:,:)  +  dt2*( 3.d0*explicit_rhs(:,:,:,1,MM0)    &
                          -        explicit_rhs(:,:,:,1,MM1) )
  endif
                     
  v(:,:,:)  =  v(:,:,:)  +  dt2*( 3.d0*explicit_rhs(:,:,:,2,MM0)    &
                         -        explicit_rhs(:,:,:,2,MM1) )
                     
  w(:,:,:)  =  w(:,:,:)  +  dt2*( 3.d0*explicit_rhs(:,:,:,3,MM0)    &
                         -        explicit_rhs(:,:,:,3,MM1) )
                     
  s1(:,:,:) = s1(:,:,:)  +  dt2*( 3.d0*explicit_rhs(:,:,:,4,MM0)    &
                         -        explicit_rhs(:,:,:,4,MM1) )
  
  if( do_second_scalar ) then
   s2(:,:,:) = s2(:,:,:)  +  dt2*( 3.d0*explicit_rhs(:,:,:,5,MM0)    &
                          -        explicit_rhs(:,:,:,5,MM1) )
  endif
  if(AB_order>2) step_flag='AB3'

 elseif( step_flag == 'AB3' ) then
 
  if(myid==0 .and. istep==2) write(0,*) '... 3rd order AB step '
  
   if( nx>1 .or. Rot(1) .ne. 0.d0 ) then
    u(:,:,:)  =  u(:,:,:)  +  dt12*( 23.d0*explicit_rhs(:,:,:,1,MM0)    &
                           -16.d0*explicit_rhs(:,:,:,1,MM1)             &
                           + 5.d0*explicit_rhs(:,:,:,1,MM2) )
   endif
  
   v(:,:,:)  =  v(:,:,:)  +  dt12*( 23.d0*explicit_rhs(:,:,:,2,MM0)    &
                          -16.d0*explicit_rhs(:,:,:,2,MM1)             &
                          + 5.d0*explicit_rhs(:,:,:,2,MM2) )
 
   w(:,:,:)  =  w(:,:,:)  +  dt12*( 23.d0*explicit_rhs(:,:,:,3,MM0)    &
                          -16.d0*explicit_rhs(:,:,:,3,MM1)             &
                          + 5.d0*explicit_rhs(:,:,:,3,MM2) )
                    
  s1(:,:,:)  =  s1(:,:,:)  +  dt12*( 23.d0*explicit_rhs(:,:,:,4,MM0)    &
                           -16.d0*explicit_rhs(:,:,:,4,MM1)             &
                           + 5.d0*explicit_rhs(:,:,:,4,MM2) )                   
  if( do_second_scalar ) then
   s2(:,:,:)  =  s2(:,:,:)  + dt12*( 23.d0*explicit_rhs(:,:,:,5,MM0)    &
                            -16.d0*explicit_rhs(:,:,:,5,MM1)            &
                            + 5.d0*explicit_rhs(:,:,:,5,MM2) )
  endif
  if(AB_order>3) step_flag='AB4'

 elseif( step_flag == 'AB4' ) then

  if(myid==0) write(0,*) '... 4th order AB step '

   if( nx>1 .or. Rot(1) .ne. 0.d0 ) then
    u(:,:,:)  =  u(:,:,:)  +  dt24*( 55.d0*explicit_rhs(:,:,:,1,MM0)   &
                           -59.d0*explicit_rhs(:,:,:,1,MM1)            &
                           +37.d0*explicit_rhs(:,:,:,1,MM2)            &
                           - 9.d0*explicit_rhs(:,:,:,1,MM3)     )
   endif
  
   v(:,:,:)  =  v(:,:,:)  +  dt24*( 55.d0*explicit_rhs(:,:,:,2,MM0)   &
                          -59.d0*explicit_rhs(:,:,:,2,MM1)            &
                          +37.d0*explicit_rhs(:,:,:,2,MM2)            &
                          - 9.d0*explicit_rhs(:,:,:,2,MM3)     )
  
   w(:,:,:)  =  w(:,:,:)  +  dt24*( 55.d0*explicit_rhs(:,:,:,3,MM0)   &
                          -59.d0*explicit_rhs(:,:,:,3,MM1)            &
                          +37.d0*explicit_rhs(:,:,:,3,MM2)            &
                          - 9.d0*explicit_rhs(:,:,:,3,MM3)     ) 
 
   s1(:,:,:)  =  s1(:,:,:)  +  dt24*( 55.d0*explicit_rhs(:,:,:,4,MM0) &
                            -59.d0*explicit_rhs(:,:,:,4,MM1)          &
                            +37.d0*explicit_rhs(:,:,:,4,MM2)          &
                            - 9.d0*explicit_rhs(:,:,:,4,MM3)     )
  
  if( do_second_scalar ) then
   s2(:,:,:)  =  s2(:,:,:)  +  dt24*( 55.d0*explicit_rhs(:,:,:,5,MM0) &
                            -59.d0*explicit_rhs(:,:,:,5,MM1)          &
                            +37.d0*explicit_rhs(:,:,:,5,MM2)          &
                            - 9.d0*explicit_rhs(:,:,:,5,MM3)     )
  endif
  endif
 
 !!=================================================================
 !!=================================================================
 !! Handle the (1/PrRe) d2/dz2 ( s_bar ) terms analytically
 !! since they are time independent and shouldn't be advanced
 !! using an explicit method. Another way would be to add the
 !! terms to the Diffusive terms before doing the implicit solves,
 !! but again, since we can do them analytically, why not?
 !!=================================================================
 !!=================================================================
 !xx=1.d0/(Re*Pr(1))
 xx = diff_coeffs(3,4)    ! z direction, variable id = 4

  if( .not. high_order_z ) then
   do k=1,array_size(KDIM,YBLOCK,myid)
    k_global = global_z_indices(START,YBLOCK,myid) + k - 1
    s1(:,:,k) = s1(:,:,k)  + dt * s1_bar(k_global,3)*xx
   enddo
  endif
 
 if( .NOT. do_second_scalar ) return

  !xx=1.d0/(Re*Pr(2))
  xx = diff_coeffs(3,5)    ! z direction, variable id = 5

   if( .not. high_order_z ) then
    do k=1,array_size(KDIM,YBLOCK,myid)
     k_global = global_z_indices(START,YBLOCK,myid) + k - 1
     s2(:,:,k) = s2(:,:,k)  + dt * s2_bar(k_global,3)*xx
    enddo
   endif
   
 return
end subroutine explicit_step


