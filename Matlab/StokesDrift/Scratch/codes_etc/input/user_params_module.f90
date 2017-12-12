!!=============================================
!!   Module for user parameters and variables
!!=============================================

  module user_params
   real(kind=8),parameter       :: delta_U = 0.05      !! velocity difference [m/s] 
   real(kind=8),parameter       :: delta_rho = 0.50    !! density difference [kg/m3] 
   real(kind=8),parameter                  :: delta = 0.2        !! thickness of density interface [m]
   real(kind=8),parameter       :: h = 0.10            !! thickness of velocity interface [m]
   real(kind=8),parameter       :: d = 117.1875            !! initial separation between lagrangian particles [m]
   real(kind=8),parameter       :: amp_pert = 1.e-4    !! ampl of density perturbations in ics [kg/m3]
   real(kind=8)                 :: z0                  !! height of max shear  [m]             
   real(kind=8)                 :: zr                  !! height of max density gradient  [m]  (z0-d)
!!=============================================
!! variables used in forcing subroutine
   real(kind=8),parameter       :: z_damping = 1.e-2  !! damping rate at bottom of domain
   real(kind=8),parameter       :: h_damping = 1.e-9  !! damping rate in the horizontal
!!=============================================
!!=============================================
!! variables used in user_ics
   real(kind=8)                 :: up_x0 = 0.5        !! location of eddy center as a fraction of Lx
   real(kind=8)                 :: up_y0 = 0.5        !! location of eddy center as a fraction of Ly
   real(kind=8)                 :: up_z0 = 1.         !! location of eddy center as a fraction of Lz
   real(kind=8)                 :: up_U0 = 3.e-5         !! amplitude of eddy at the surface
   real(kind=8)                 :: up_alpha = 1.e-9   !! [1/m^2]
   real(kind=8)                 :: up_beta = 4.e-6    !! [1/m^2]
!!  info for restarting simulations 
!!=============================================
logical, parameter              :: advect_tracer_linear=.TRUE.
logical, parameter              :: up_restart=.TRUE.
character(len=80), parameter    :: up_rs_basename='RESTART/XYZ_039400' !relative to $BASE directory
integer, parameter              :: up_rs_islice = 1 ! time slice (important if restart file is from appended file)
logical, parameter              :: up_subtract_s1_bar = .FALSE. !TRUE if s1_bar is added to s1 in restart files
logical, parameter              :: up_subtract_s2_bar = .FALSE. !TRUE if s2_bar is added to s2 in restart files

end module user_params

! "up_" stands for "user parameter"... to ensure that internal variable names don't get inadvertently duplicated
