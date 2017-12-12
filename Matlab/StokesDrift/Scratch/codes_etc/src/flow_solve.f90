!!-------------------------------------------------------------!!
!!  Stratified flow solver.                                    !!
!!                                                             !!
!!    MPI parallel parallel version                            !!
!!    Two-dimensional data/'processor' decomposition           !!
!!                                                             !!
!!    Spatial differentiation via fourier, sin, cos expansions !!
!!                                                             !!
!!    1 or 2 scalars, configured as combinations of density,   !!
!!    temperature, salinity or passive tracers.                !!
!!                                                             !!
!!    2nd, 3rd or 4th order Adams Bashforth time stepping      !!
!!    for advection, buoyancy, rotation, and forcing.          !!
!!                                                             !!
!!    3rd or 4th order Adams Moulton implicit time stepping    !!
!!    for the diffusive terms.                                 !!
!!                                                             !!
!!    Generalized diffusion operators:                         !!
!!      'hyper' diffusion in x,y,z w/ independent coeffs       !!
!!      default settings: laplacian w/ constant coeff          !!
!!                                                             !!
!!    Immersed boundary z = h(x,y)                             !!
!!                                                             !!
!!    Rotation (f or beta plane approximations)                !!
!!             (horizontal component of Coriolis optional)     !!
!!                                                             !!
!!                                                             !!
!!    External software:                                       !!
!!     MPI               (message passing, optional)           !!
!!     BLAS/LAPACK       (sequential Ax=b solves etc)          !!
!!     FFTW3             (sequential spectral transforms)      !!
!!     PD3FFT            (MPI data transposes)                 !!
!!     NETCDF            (output data format)                  !!
!!                                                             !!
!!    Author: K. B. Winters                                    !!
!!            Scripps Institution of Oceanography              !!
!!                         and                                 !!
!!            Mechanical and Aerospace Engineering             !!
!!            University of California, San Diego              !!
!!-------------------------------------------------------------!!  
program flow_solve
use etc
use mpi_params
use timing
call initialize
 do while (integrate)
  call write_results
   t_start_time_step = mpi_wtime()
    call eval_explicit_rhs
     call apply_forcing
      call explicit_step
       call implicit_step
        call immersed_bdry  
         call pressure_projection
          call immersed_bdry
           call advance_particles
            call toggle
             enddo
end program flow_solve



subroutine eval_explicit_rhs
use methods_params, only: do_second_scalar
implicit none 
 call equation_of_state
  call explicit_rhs_1
   call explicit_rhs_2
    call explicit_rhs_3
     call explicit_rhs_4
      if( do_second_scalar ) then
       call explicit_rhs_5
      endif
end subroutine eval_explicit_rhs



subroutine implicit_step
use methods_params, only: do_second_scalar
implicit none
  call implicit_solve(1)
   call implicit_solve(2)
    call implicit_solve(3)
     call implicit_solve(4)
      if( do_second_scalar ) then
       call implicit_solve(5)
      endif
end subroutine implicit_step



