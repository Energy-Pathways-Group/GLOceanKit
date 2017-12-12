subroutine PreliminaryTasks
 use etc
 use dimensional_scales
 use pde_params
 use methods_params
 use mpi_params
 use dependent_variables
 use independent_variables
 use particles
 implicit none  
 integer          :: id,idim,kend
 integer          :: order1=1
 integer          :: order2=2
 real(kind=8)     :: pi

 
 pi=4.d0*datan(1.d0)
 if(myid==0) then
  open(1,file=logfile,position='append')
  write(1,*) ' =========================================================== '
  write(1,*) ' =========================================================== '
  write(1,*) '                  PreliminaryTasks Report:'
  write(1,*) ' =========================================================== '
  write(1,*) ' =========================================================== '
  write(1,*) '  ' 
  write(0,*) ' ................'
  write(0,*) ' ................     hello world from PreliminaryTasks'
 endif 
 
 !------------------------------------------------------------
 ! invoke some specialized IB logic for special runlabels
 !------------------------------------------------------------
 if( do_immersed_boundary ) then
  if( trim(runlabel)=='hc_temp_bcs' ) then
   if(myid==0) write(0,*) ' ................      special IB logic ',trim(runlabel)
   hc_IB_test = .TRUE.          ! hc_IB IB test case (upper flat boundary)
   no_slip_no_flux_IB=.FALSE.   ! generic (no slip, no scalar flux) IB case  
  elseif( trim(runlabel)=='stokes2' ) then
   stokes2 = .TRUE.             ! stokes2 IB test case
   no_slip_no_flux_IB=.FALSE.   ! generic (no slip, no scalar flux) IB case
  endif
 endif
 
 if( vertical_coriolis ) then
  f_plane_key='non-traditional'
  if(myid==0) then
   write(0,*) ' ................'
   write(0,*) '                   =====> USING NONTRADITIONAL APPROXIMATION FOR CORIOLIS'
   write(0,*) ' ................'
  endif
 else
  f_plane_key='traditional'
 endif
 
 if( .NOT. do_nonlinear ) then
  if(myid==0) then
   write(0,*) ' ................'
   write(0,*) '                   =====> SUPPRESSING NONLINEAR TERMS U DOT GRAD PHI FOR ALL VARIABLES PHI'
   write(0,*) ' ................'
  endif
 endif

 if( .NOT. do_second_scalar ) then
  if(myid==0) then
   write(0,*) ' ................'
   write(0,*) '                   =====> SUPPRESSING CALCULATION OF SECOND SCALAR FIELD'
   write(0,*) ' ................'
  endif
 endif
 
 !!===========================================================
 !! look for inconsistencies in bcs & methods
 !!===========================================================
 do id=1,6
  do idim=1,3
   do kend=1,2
   
    if( trim(deriv_type(id,idim,order1)) == 'fourier' ) then
     if( trim(bc_type(id,idim,kend)) .ne. 'periodic' ) &
     stop 'fourier/periodic mismatch'
    endif
    
    if( trim(deriv_type(id,idim,order1)) == 'sin' ) then
     if( trim(bc_type(id,idim,kend)) .ne. 'dirichlet' ) then
      write(0,*) id,idim,kend
      write(0,*) trim(deriv_type(id,idim,order1))
      write(0,*) trim(bc_type(id,idim,kend))
      stop 'sin/dirichlet mismatch'
     endif
    endif
    
    if( trim(deriv_type(id,idim,order1)) == 'cos' ) then
     if( trim(bc_type(id,idim,kend)) .ne. 'neumann' ) then
      write(0,*) id,idim,kend
      write(0,*) trim(deriv_type(id,idim,order1))
      write(0,*) trim(bc_type(id,idim,kend))
      stop 'cos/neumann mismatch'
     endif
    endif
        
   enddo
  enddo
 enddo
 

 !!===========================================================
 !! derivative types for first derivs have been defined 
 !! via user input, now set second derivatives
 !!===========================================================
 do id=1,6
  do idim=1,3
   if( trim(deriv_type(id,idim,order1)) == 'fourier' ) then
    deriv_type(id,idim,order2)='fourier'
   elseif( trim(deriv_type(id,idim,order1)) == 'sin' ) then
    deriv_type(id,idim,order2)='cos'
   elseif( trim(deriv_type(id,idim,order1)) == 'cos' ) then
    deriv_type(id,idim,order2)='sin'
   elseif( trim(deriv_type(id,idim,order1)) == 'cheby' ) then
    deriv_type(id,idim,order2)='cheby'
   elseif( trim(deriv_type(id,idim,order1)) == 'compact' ) then
    deriv_type(id,idim,order2)='compact'
   endif
  enddo
 enddo
 
        
!!===========================================================
!! explicit methods start @ 1st order
!!===========================================================
 step_flag='euler'
 particle_step_flag='euler'
 
!!===========================================================
!! I'll set the explicit & implicit
!! timestepping orders for now
!!===========================================================
 AB_order=3         !! k=3 step method error=O(h^4)
 AM_order=AB_order  !! k=2 step method error=O(h^4)
 
!!===========================================================
!! initial values for cyclic counters
!!===========================================================
 MM0=1;MM1=2;MM2=3;MM3=4  !! for AB timestepping
 PM0=1;PM1=2;PM2=3;PM3=4  !! for AB stepping particles/grid
 N=1;NM1=2                !! for AM timestepping
 M_oldest=AB_order
 P_oldest=PM3
 N_oldest=AM_order-2
 istart = nint(t0/dt)
 iend = istart + nint( (tf-t0)/dt )
 t_secs = t0 
 istep = istart
 integrate = .TRUE.
 
!!===========================================================
!!Dimensional scales
!!===========================================================
 velocity_scale = u0
 length_scale = Lz
 time_scale = length_scale/u0
 density_scale = length_scale*DGRAD
 pressure_scale = rho_0*velocity_scale**2
 bfreq = sqrt((g/rho_0)*DGRAD) 
 omega_earth = 2.d0*pi/(23.d0*3600.d0 + 56.d0*60.d0 + 4.d0)
 latitude=asin( f(1)/(2.d0*omega_earth) )
 
!!===========================================================
!!Dimensionless parameters.
!!
!! N.B.
!! if a scalar is pd, scale it consistently with
!! how I do the perturbation density in the w eqn,
!! ====> in this case, ignore user supplied scale
!!===========================================================
 if( scalar_kind(1)=='r' ) then
  scalar_scale(1)=length_scale*dgrad
 endif
 if( scalar_kind(2)=='r' ) then
  scalar_scale(2)=length_scale*dgrad
 endif
 if( scalar_kind(1)=='r' .and. scalar_kind(2)=='r' ) then
  stop 'cant have both scalars set to density'
 endif
 
 !! f, beta,  stored in f(1), f(2)
 !!  dless versions (+ nontraditional f~) in Rot(1), Rot(2), Rot(3)
 Rot(1) = f(1)*length_scale/U0
 Rot(2) = f(2)*(length_scale*length_scale/u0) 
 
 Rot(3) = 2.d0*omega_earth*cos(latitude)*length_scale/U0
 
 N2 = (g/rho_0)*dgrad
 if( scalar_kind(1)=='p' .and. scalar_kind(2)=='p' ) then
  Ri=0.d0
 elseif( scalar_kind(1)=='p' .and. .not. do_second_scalar ) then
  Ri=0.d0
 else
  Ri = N2*length_scale**2/(U0**2)
 endif
 if( nu > 0.d0 ) then
  Re = U0*length_scale/nu
 else
  Re = -999   ! special value used internally
 endif
 Pr(1) = nu/kappa(1)
 Pr(2) = nu/kappa(2)
 
 !!===========================================================
 !!   ===> from here on... dt is DIMENSIONLESS
 !!===========================================================
 dt = dt/time_scale   !! time step now dimensionless

 
 !!===========================================================
 !!    make sure f~=0 of not 'non-traditional'		 
 !!===========================================================
 if(f_plane_key .ne. 'non-traditional') then
  Rot(3)=0.d0
 endif
 
 
 !!===========================================================
 !!    Write out parameters characterizing simulation.
 !!===========================================================
 if( myid == 0 ) then
  do id=0,1
   write(id,*) '................       bulk velocity scale U                   [m/s]:  ',U0
   write(id,*) '................       characteristic density gradient scale [kg/m4]:  ',dgrad
   write(id,*) '................       s1 r=rho, t=temp, s=salinity, p=passive      :   ',trim(scalar_kind(1))
   if( do_second_scalar ) & 
   write(id,*) '................       s2 r=rho, t=temp, s=salinity, p=passive      :   ',trim(scalar_kind(2)) 
   write(id,*) '................       characteristic range of scalar 1 [e.g. deg C]:  ',scalar_scale(1)
   if( do_second_scalar ) &
   write(id,*) '................       characteristic range of scalar 2   [e.g. psu]:  ',scalar_scale(2)
   write(id,*) '................       characteristic density                [kg/m3]:  ',rho_0 
   write(id,*) '................       latitude in degrees                          :  ',latitude*360/(2.*pi)
   write(id,*) '................       coriolis parameter f                    [1/s]:  ',f(1) 
   write(id,*) '................       coriolis parameter f~                   [1/s]:  ',Rot(3)*velocity_scale/length_scale
   write(id,*) '................       beta plane parameter                   [1/ms]:  ',f(2)
   if( f(1) .ne. 0 )  &
   write(id,*) '................       inertial period                          [s] :  ',2.*pi/f(1)
   if( f(2) .ne. 0 )  &
   write(id,*) '................       beta plane pivot                         [m] :  ',y_pivot
   write(id,*) '................       gravitational acceleration g           [m/s2]:  ',g
   write(id,*) '................       molecular diffusivity of scalar 1      [m2/s]:  ',kappa(1)
   if( do_second_scalar ) &
   write(id,*) '................       molecular diffusivity of scalar 2      [m2/s]:  ',kappa(2)
   write(id,*) '................       viscosity                              [m2/s]:  ',nu
   write(id,*) '................       time step dt                             [s] :  ',dt*time_scale
   write(id,*) '................       initial time                             [s] :  ',t0
   write(id,*) '................       final time                               [s] :  ',tf
   write(id,*) '................       Lz/U0                                    [s] :  ',time_scale
   write(id,*) '................       number of steps                              :  ',iend-istart+1
   write(id,*) '................       inverse Rossby number                fL/U    :  ',Rot(1) 
   write(id,*) '................       dimensionless beta parameter                 :  ',Rot(2)
   write(id,*) '................       dimensionless f~ parameter                   :  ',Rot(3)
   write(id,*) '................       bulk Richardson number (NL/U)^               :  ',Ri
   write(id,*) '................       bulk Reynolds number UL/nu                   :  ',Re
   write(id,*) '................       Prandtl number  nu/kappa_1                   :  ',Pr(1)
   if( do_second_scalar ) &
   write(id,*) '................       Prandtl number  nu/kappa_2                   :  ',Pr(2)
   if( N2 .ne. 0 )  then
    write(id,*) '................       N                                      [1/s] :  ',sqrt(N2)
    write(id,*) '................       buoyancy period                          [s] :  ',2.*pi/sqrt(N2)
   endif
   if( f(1) .ne. 0 )  then
    write(id,*) '................       N/f                                      [1] :  ',sqrt(N2)/f(1)
    write(id,*) '................       sqrt(nu/f)                               [m] :  ',sqrt(nu)/f(1) 
   endif
   write(id,*) '................'
  enddo
  write(1,*) ' -----> PreliminaryTasks routine exiting normally  <---------- '
  write(1,*) '  '
  write(1,*) '  '
  close(1)
 endif
  
 return
end subroutine PreliminaryTasks


