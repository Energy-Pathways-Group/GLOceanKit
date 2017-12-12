!!=============================================
!!   Data Modules
!!=============================================

  module mpi_params
   integer,parameter     :: maxprocs=1024
   integer               :: myid
   integer               :: comm
   integer               :: ierr
   integer               :: numprocs
   integer               :: nreqs
   integer               :: source
   integer               :: dest
   integer               :: tag   
  end module mpi_params

  
  module decomposition_params
   integer                   :: np   ! np=p1*p2
   integer                   :: p1
   integer                   :: p2
  !------------------------------------------------------------------------------------------
  ! mem_order(1,XBLOCK) = the fortran array index corresponding to x in XBLOCK decomposition
  ! mem_order(2,XBLOCK) = the fortran array index corresponding to y in XBLOCK decomposition
  ! mem_order(3,XBLOCK) = the fortran array index corresponding to z in XBLOCK decomposition 
  !------------------------------------------------------------------------------------------
   integer                   :: mem_order(3,3)  
  !-----------------------------------------------------------------------------------------
  !  e.g. layout(XCOORD,ZBLOCK) defines the number of processors used in the decomposition of
  !                             x when the global data set is arranged in ZBLOCKS
  !-----------------------------------------------------------------------------------------
   integer                   :: layout(3,3)
  !-----------------------------------------------------------------------------------------
  !  e.g. proc_row(YBLOCK,pid) defines the row assignment of pid in the n1xn2 processor grid
  !                            when the global data set is arranged in YBLOCKS
  !-----------------------------------------------------------------------------------------
   integer,allocatable       :: proc_row(:,:)
  !-----------------------------------------------------------------------------------------
  !  e.g. proc_col(YBLOCK,pid) defines the col assignment of pid in the n1xn2 processor grid
  !                            when the global data set is arranged in YBLOCKS
  !-----------------------------------------------------------------------------------------
   integer,allocatable       :: proc_col(:,:)
  !-----------------------------------------------------------------------------------------
  !  e.g. array_size(JDIM,XBLOCK,pid) defines the number of grid points in the local data array
  !                                   on pid in the JDIM storage index in the XBLOCK decomposition
  !-----------------------------------------------------------------------------------------
   integer,allocatable       :: array_size(:,:,:) 
  !-----------------------------------------------------------------------------------------
  !  e.g. global_x_indices(END,XBLOCK,pid) defines the the LAST global x index for the local
  !                                   data array on pid in the XBLOCK decomposition
  !-----------------------------------------------------------------------------------------
   integer,allocatable       :: global_x_indices(:,:,:)
   integer,allocatable       :: global_y_indices(:,:,:)
   integer,allocatable       :: global_z_indices(:,:,:)
  !-----------------------------------------------------------------------------------------
  !  these are simply conventions and should never change
  !  (to keep me from having to use 1,2,3 in different contexts)
  !-----------------------------------------------------------------------------------------
   integer,parameter         :: XCOORD=1
   integer,parameter         :: YCOORD=2
   integer,parameter         :: ZCOORD=3
   integer,parameter         :: XBLOCK=1
   integer,parameter         :: YBLOCK=2
   integer,parameter         :: ZBLOCK=3
   integer,parameter         :: IDIM=1   ! convention: fortran arrays are indexed (IDIM,JDIM,KDIM)
   integer,parameter         :: JDIM=2   ! i.e. IDIM is fastest varying array index, then JDIM,KDIM
   integer,parameter         :: KDIM=3
   integer,parameter         :: START=1
   integer,parameter         :: END=2

   integer,allocatable       :: jk_indices_YBLOCK(:,:) 
   integer,allocatable       :: ijk_indices_YBLOCK(:,:)  
   integer,allocatable       :: ijk_indices_ZBLOCK(:,:)
   
  end module decomposition_params
 
 
  module independent_variables
  !---------------------------------------------------------------
  ! except for cheby, s is equally space computational coordinate
  !---------------------------------------------------------------
   integer                         :: nx     ! global number of grid points
   integer                         :: ny     ! global number of grid points
   integer                         :: nz     ! global number of grid points
   real(kind=8)                    :: Lx     ! always [m]
   real(kind=8)                    :: Ly     ! always [m]
   real(kind=8)                    :: Lz     ! always [m]
   real(kind=8)                    :: dt     ! nondimensionalized after user input, time_scale
   real(kind=8)                    :: t0     ! ""
   real(kind=8)                    :: tf     ! ""
   real(kind=8)                    :: t_secs ! dimensional value [s]
 !---------------------------------------------------------------
 ! except for cheby, s is equally space computational coordinate
 ! these functions define the grid stretching for chain rule etc
 !---------------------------------------------------------------
   real(kind=8)                    :: ds(3)       ! ignored in cheby directions
   real(kind=8),allocatable,target :: x(:)        ! x(s), nx values
   real(kind=8),allocatable,target :: x_s(:)      ! dx/ds
   real(kind=8),allocatable,target :: x_ss(:)     ! d2x/ds2
   real(kind=8),allocatable,target :: y(:)        ! y(s), ny values
   real(kind=8),allocatable,target :: y_s(:)      ! dy/ds
   real(kind=8),allocatable,target :: y_ss(:)     ! d2y/ds2
   real(kind=8),allocatable,target :: z(:)        ! z(s), nz values
   real(kind=8),allocatable,target :: z_s(:)      ! dz/ds
   real(kind=8),allocatable,target :: z_ss(:)     ! d2z/ds2
   real(kind=8),allocatable,target :: s_of_x(:)   ! 's' values in x direction
   real(kind=8),allocatable,target :: s_of_y(:)   ! 's' values in y direction
   real(kind=8),allocatable,target :: s_of_z(:)   ! 's' values in z direction
   real(kind=8),allocatable,target :: s_x(:)      ! ds/dx = d/dx (s_of_x)
   real(kind=8),allocatable,target :: s_y(:)      ! ds/dy = d/dy (s_of_y)
   real(kind=8),allocatable,target :: s_z(:)      ! ds/dz = d/dz (s_of_z)
   real(kind=8)                    :: cmap_params(5,3)  ! 5 parameter mapping in each direction
  !---------------------------------------------------------------------
  ! convenient to have some logical variables to avoid unnecessary work
  !---------------------------------------------------------------------
   logical                         :: stretch(3)        ! is chain rule required?
   logical                         :: xdim_periodic     ! is x periodic over Lx
   logical                         :: ydim_periodic     ! is y periodic over Ly
   logical                         :: zdim_periodic     ! is z periodic over Lz
  end module independent_variables




  module differentiation_params
   
   real(kind=8),allocatable,target :: kx(:)
   real(kind=8),allocatable,target :: ky(:)
   real(kind=8),allocatable,target :: kz(:)
   real(kind=8),allocatable,target :: kxfilter(:)
   real(kind=8),allocatable,target :: kyfilter(:)
   real(kind=8),allocatable,target :: kzfilter(:)
   integer(kind=8)                 :: fourier_plan(3,2)
   integer(kind=8)                 :: cos_plan(3)
   integer(kind=8)                 :: sin_plan(3)
   
   integer(kind=8)                 :: xy_plan(2,2)
   logical                         :: fourier_done(3)=.FALSE.
   logical                         :: sin_done(3)=.FALSE.
   logical                         :: cos_done(3)=.FALSE.
   
   !----------------------------------------------------------------
   !  differentiate is a routine from my differentiation_toolbox
   !  that is more general than needed in this fluids code...
   !  for compatability, the following need to be declared but
   !  the fluids code will not allocate memory or use these arrays
   !----------------------------------------------------------------
     real(kind=8),allocatable,target :: c6_Ax(:,:)
     real(kind=8),allocatable,target :: c6_Ay(:,:)
     real(kind=8),allocatable,target :: c6_Az(:,:)
     integer,allocatable,target      :: c6_ipivx(:)
     integer,allocatable,target      :: c6_ipivy(:)
     integer,allocatable,target      :: c6_ipivz(:)
     real(kind=8),allocatable,target :: A_filter_x(:,:)
     real(kind=8),allocatable,target :: A_filter_y(:,:)
     real(kind=8),allocatable,target :: A_filter_z(:,:)
     integer,allocatable,target      :: ipiv_filter_x(:)
     integer,allocatable,target      :: ipiv_filter_y(:)
     integer,allocatable,target      :: ipiv_filter_z(:)
     real(kind=8),target             :: alpha_f_x
     real(kind=8),target             :: alpha_f_y
     real(kind=8),target             :: alpha_f_z
     integer(kind=8)                 :: cheby_plan(3)
     logical                         :: cheby_done(3)=.FALSE.
     logical                         :: compact_done(3)=.FALSE.
   !----------------------------------------------------------------
   
  end module differentiation_params


  
  module dimensional_scales  
   real(kind=8)             :: f(2)  !! f and beta
   real(kind=8)             :: latitude
   real(kind=8)             :: y_pivot
   real(kind=8)             :: omega_earth
   real(kind=8)             :: g
   real(kind=8)             :: rho_0, rho0
   real(kind=8)             :: nu
   real(kind=8)             :: dgrad !! char val |d/dz rho|, > 0
   real(kind=8)             :: bfreq
   real(kind=8)             :: N2
   real(kind=8)             :: kappa(2)
   real(kind=8)             :: scalar_scale(2)
   real(kind=8)             :: pressure_scale
   real(kind=8)             :: time_scale
   real(kind=8)             :: velocity_scale, u0
   real(kind=8)             :: length_scale
   real(kind=8)             :: density_scale
  end module dimensional_scales




  module pde_params
   character(len=80)               :: bc_type(6,3,2)
   character(len=80)               :: f_plane_key
   real(kind=8)                    :: Ri
   real(kind=8)                    :: Rot(3)
   real(kind=8)                    :: Re
   real(kind=8)                    :: Pr(2)
  end module pde_params




  module methods_params
   character(len=80)        :: deriv_type(6,3,2)
   integer                  :: AB_order
   integer                  :: AM_order
   logical                  :: do_forcing
   logical                  :: do_immersed_boundary
   logical                  :: do_particles
   logical                  :: do_nonlinear      
   logical                  :: do_second_scalar
   logical                  :: vertical_coriolis
   logical                  :: forcing_key(5)=.TRUE. 
   !-------------------------------------------------------------------
   ! flags for specialized IB logic, detect if necessary via runlabel
   !-------------------------------------------------------------------
   logical                  :: stokes2 = .FALSE.                ! runlabel = 'stokes2'
   logical                  :: hc_IB_test = .FALSE.             ! runlabel = 'hc_temp_bcs'
   logical                  :: moving_cylinder= .FALSE.         ! runlabel = 'flow_around_cylinder'
   logical                  :: no_slip_no_flux_IB=.TRUE.        ! generic case 
  end module methods_params
 
 
 
 
  module diffusion_params
   integer                         :: p(3)                  ! 1/2 order of diffusion operators
   real(kind=8)                    :: T_diss(3)             ! dissipation time scale, used to define coeffs
   real(kind=8)                    :: delta(3)              ! decay factor = exp(-delta) after T_diss in x,y,z
   logical                         :: high_order_operators
   logical                         :: high_order_z=.TRUE.   ! enable if requested and methods are OK
   real(kind=8)                    :: diff_coeffs(3,5)      ! dless coeffs(DIR,VARIABLE_ID)
  end module diffusion_params
  
 
 
 

   		 
  module dependent_variables
   real(kind=8),allocatable,target      :: u(:,:,:)
   real(kind=8),allocatable,target      :: v(:,:,:)
   real(kind=8),allocatable,target      :: w(:,:,:)
   real(kind=8),allocatable,target      :: s1(:,:,:)
   real(kind=8),allocatable,target      :: s2(:,:,:)
   real(kind=8),allocatable,target      :: s1_bar(:,:)   !! ===>(nz,3)
   real(kind=8),allocatable,target      :: s2_bar(:,:)
   character(len=1)                     :: scalar_kind(2)
  end module dependent_variables


  module intermediate_variables
  !--------------------------------------------
  !  tmp arrays in each decomposition format
  !--------------------------------------------
   real(kind=8),allocatable,dimension(:,:,:,:),target    :: tmpX         ! XBLOCK format
   real(kind=8),allocatable,dimension(:,:,:,:),target    :: tmpY         ! YBLOCK format
   real(kind=8),allocatable,dimension(:,:,:,:),target    :: tmpZ         ! ZBLOCK format
  !-------------------------------------------- 
  ! 3d arrays in YBLOCK format
  !--------------------------------------------
   real(kind=8),allocatable,dimension(:,:,:),target      :: pd           ! YBLOCK format
   real(kind=8),allocatable,dimension(:,:,:),target      :: phi          ! YBLOCK format
   real(kind=8),allocatable,dimension(:,:,:),target      :: div_u        ! YBLOCK format
   real(kind=8),allocatable,dimension(:,:,:,:,:),target  :: explicit_rhs ! YBLOCK format
   real(kind=8),allocatable,dimension(:,:,:,:,:),target  :: implicit_rhs ! YBLOCK format
  !-------------------------------------------- 
  ! 2d arrays for boundary condition data
  !--------------------------------------------
   real(kind=8),allocatable,dimension(:,:,:),target      :: gamma_xy     ! 2d global, (x,y,2)
   real(kind=8),allocatable,dimension(:,:,:),target      :: gamma_yx     ! 2d global, (y,x,2)
   real(kind=8),allocatable,dimension(:,:,:),target      :: bvals        ! k-collapsed ZBLOCK, 2
  !-------------------------------------------- 
  ! 2d arrays for implicit Helmholtz terms
  !--------------------------------------------
   real(kind=8),allocatable,dimension(:,:,:),target      :: q0        ! k-collapsed ZBLOCK + fid
   real(kind=8),allocatable,dimension(:,:),target        :: q1        ! k-collapsed ZBLOCK
   
  end module intermediate_variables
  
  
  
  
   module etc
   character(len=80)     :: runlabel
   character(len=80)     :: logfile='output/logfile'
   character(len=80)     :: memoryfile='output/memlog'
   character(len=80)     :: message
   character(len=80)     :: step_flag
   integer               :: istep
   integer               :: istart
   integer               :: iend
   integer               :: MM0,MM1,MM2,MM3  !! for AB timestepping
   integer               :: PM0,PM1,PM2,PM3  !! for AB stepping particles/grid
   integer               :: N,NM1            !! for AM timestepping
   integer               :: M_oldest
   integer               :: P_oldest
   integer               :: N_oldest 
   logical               :: integrate=.TRUE.
   real(kind=8)          :: mymax(6,2)
   real(kind=8)          :: mymin(6,2)
  end module etc
  
  
  
  
  module particles
   logical                           :: user_defines_particles=.FALSE.
   character(len=80)                 :: particle_step_flag   
   real(kind=8),allocatable          :: positions(:,:)
   real(kind=8),allocatable,target   :: uvels    (:,:)
   real(kind=8),allocatable,target   :: vvels    (:,:)
   real(kind=8),allocatable,target   :: wvels    (:,:)  
   integer,allocatable               :: my_labels(:)
   integer                           :: nparticles
   integer                           :: j_particle_time
   integer                           :: istart_trajectories=0
   integer                           :: particle_write_inc=1
   integer                           :: my_first
   integer                           :: my_last
  end module particles
  
 
 module io_params
   ! data read in from input/io_params
   integer,parameter             :: maxsets=2048
   integer                       :: num_file_sets             ! user supplied
   character(len=80)             :: filename_root(maxsets)    ! user supplied
   character(len=80)             :: mode(maxsets)             ! user supplied
   integer                       :: ilocs(3,maxsets)          ! user supplied    
   integer                       :: jlocs(3,maxsets)          ! user supplied   
   integer                       :: klocs(3,maxsets)          ! user supplied
   integer                       :: variable_key(9,maxsets)   ! user supplied
   integer                       :: nsteps(maxsets)           ! user supplied      
   
   ! nonsingleton indices arranged 1st in the contiguous array 'outdata' that
   ! contains only and all the data to be actually written to the netcdf file
   integer                       :: dimid_x(maxsets)  ! 'outdata' array dimension for x in nc file
   integer                       :: dimid_y(maxsets)  ! 'outdata' array dimension for x in nc file
   integer                       :: dimid_z(maxsets)  ! 'outdata' array dimension for x in nc file
   integer                       :: dimid_t(maxsets)  ! 'outdata' array dimension for x in nc file
   
   ! these arrays are used in the actual netcdf write statements
   ! to define how to interpret the layout of the 'outdata' array
   integer                       :: nspace(maxsets)        ! number of nonsingleton space dimensions
   integer                       :: time_counter(maxsets)  ! time index
   integer                       :: count(4,maxsets)
   
   ! filename constructed given its root, num of dimensions and topdir
   character(len=80)             :: fullname(maxsets)
   
   ! indices, counts of local YBLOCK data to extract for writing
   integer                       :: my_x0(maxsets)    ! initial local index for x coordinate
   integer                       :: my_x1(maxsets)    ! final   local index for x coordinate
   integer                       :: my_xinc(maxsets)  ! increment for x coordinate sampling
   integer                       :: my_nx(maxsets)    ! number of x coord vals to be written
   
   integer                       :: my_y0(maxsets)    ! initial local indey for y coordinate
   integer                       :: my_y1(maxsets)    ! final   local indey for y coordinate
   integer                       :: my_yinc(maxsets)  ! increment for y coordinate sampling
   integer                       :: my_ny(maxsets)    ! number of y coord vals to be written
   
   integer                       :: my_z0(maxsets)    ! initial local indez for z coordinate
   integer                       :: my_z1(maxsets)    ! final   local indez for z coordinate
   integer                       :: my_zinc(maxsets)  ! increment for z coordinate sampling
   integer                       :: my_nz(maxsets)    ! number of z coord vals to be written
   
   ! key specifying whether there is data to be written 
   logical                       :: do_write(maxsets)
   logical                       :: write_s1_bar(maxsets)
   logical                       :: write_s2_bar(maxsets)
    
  end module io_params
  
  
 
 
  module interpolation
   integer,target                       :: interp_order_XBLOCK(3)
   real(kind=8),allocatable,target      :: tx_XBLOCK(:)
   real(kind=8),allocatable,target      :: ty_XBLOCK(:)
   real(kind=8),allocatable,target      :: tz_XBLOCK(:)
   real(kind=8),allocatable,target      :: wrk_XBLOCK(:)
   
   integer,target                       :: interp_order_YBLOCK(3)
   real(kind=8),allocatable,target      :: tx_YBLOCK(:)
   real(kind=8),allocatable,target      :: ty_YBLOCK(:)
   real(kind=8),allocatable,target      :: tz_YBLOCK(:)
   real(kind=8),allocatable,target      :: wrk_YBLOCK(:)
   
   integer,target                       :: interp_order_ZBLOCK(3)
   real(kind=8),allocatable,target      :: tx_ZBLOCK(:)
   real(kind=8),allocatable,target      :: ty_ZBLOCK(:)
   real(kind=8),allocatable,target      :: tz_ZBLOCK(:)
   real(kind=8),allocatable,target      :: wrk_ZBLOCK(:)
  end module interpolation
  
  
  

  module immersed_boundary 
   real(kind=8),allocatable      :: dist(:,:,:)        ! (y,x,z)
   real(kind=8),allocatable      :: nhat(:,:,:,:)      ! (y,x,z,3)  1,2,3=>x,y,z
   real(kind=8)                  :: sigma_d,sigma_n,TOL
   
   integer,allocatable           :: num_remote_image_points(:)
   integer                       :: total_num_remote_image_points
   integer                       :: my_start_index
   
   real(kind=8),allocatable      :: xyz_solid(:,:)
   real(kind=8),allocatable      :: xyz_image(:,:)
   real(kind=8),allocatable      :: image_vals(:)
   real(kind=8),allocatable      :: xtmp(:),fval(:)
   integer,allocatable           :: ijk_solid(:,:)
   integer,allocatable           :: owner_solid(:)
   integer,allocatable           :: owner_image(:)
   integer,allocatable           :: itmp(:)
   
   real(kind=8)                  :: xpt,ypt,zpt
   logical                       :: stationary_IB=.TRUE.
   real(kind=8)                  :: U_cyl,W_cyl
   
  end module immersed_boundary


  module timing
   real(kind=8)                    :: t_start_time_step
   real(kind=8)                    :: t_end_time_step
   real(kind=8)                    :: t_time_step
   real(kind=8)                    :: t_total=0.d0
   real(kind=8),external           :: mpi_wtime
  end module timing

