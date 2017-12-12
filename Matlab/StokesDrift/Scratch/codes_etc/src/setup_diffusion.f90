subroutine SetupDiffusion
!--------------------------------------------------------------------
!  set up some coefficients etc associated with diffusion operators
!--------------------------------------------------------------------
use diffusion_params
use decomposition_params
use dimensional_scales, only       : viscosity=>nu,kappa,length_scale,velocity_scale
use independent_variables, only    : nx,ny,nz,Lx,Ly,Lz,z
use pde_params, only               : Re,Pr
use mpi_params, only               : myid,comm,ierr
use etc, only                      : logfile
use intermediate_variables, only   : q0,q1
use differentiation_params, only   : kx,ky
use methods_params, only           : do_second_scalar,deriv_type

implicit none
integer                           :: i,j,k
integer                           :: dim,fid
integer                           :: npvs
integer                           :: i_global,j_global
real(kind=8)                      :: nu(3),pi,zval

if(myid==0) then
  write(0,*) ' ................'
  write(0,*) ' ................     hello world from SetupDiffusion'
  open(1,file=logfile,position='append')
  write(1,*) '  '
  write(1,*) '  '
  write(1,*) ' =========================================================== '
  write(1,*) ' =========================================================== '
  write(1,*) '                  SetupDiffusion Report:' 
  write(1,*) '============================================================ '
  write(1,*) '============================================================ '
  write(1,*) '  '
 endif
 
 if( do_second_scalar ) then
  npvs=5
 else
  npvs=4
 endif
 pi=4.d0*datan(1.d0)
 
 
!--------------------------------------------------------------------
!  default values for standard laplacian diffusion
!  with constant fluid parameters nu, kappa(1) and kappa(2) [m^2/s]
!--------------------------------------------------------------------  
  if( .not. high_order_operators ) then
   high_order_z=.FALSE.
   p(:) = (/1,1,1/)
   T_diss(:) = 1.d0    ! [s]
   delta(3) = viscosity*T_diss(3)/(Lz/(pi*nz))**2
   delta(2) = viscosity*T_diss(2)/(Ly/(pi*ny))**2
   delta(1) = viscosity*T_diss(1)/(Lx/(pi*nx))**2
   !------------------------------------------------------
   !  isotropic, nu, kappa1, kappa2, dimensionless form
   !------------------------------------------------------
   do dim=1,3
    if( Re .ne. -999 ) then
     diff_coeffs(dim,1:3) = 1.d0/Re
     diff_coeffs(dim,4) = 1.d0/(Re*Pr(1))
     if( do_second_scalar )  &
      diff_coeffs(dim,5) = 1.d0/(Re*Pr(2))
    else
     diff_coeffs(dim,1:3) = 0.d0
     diff_coeffs(dim,4) = 0.d0
     if( do_second_scalar )  &
      diff_coeffs(dim,5) = 0.d0
    endif
enddo
  
  
!--------------------------------------------------------------------
!  if high_order_operators is turned on, find out whether
!  its legal to apply high order logic in the z direction
!--------------------------------------------------------------------  
  elseif( high_order_operators ) then
   
   !--------------------------------------------------------------------------
   !  revert to 2nd deriv diffusion in z if methods are not appropriate
   !--------------------------------------------------------------------------
   if( .not. high_order_z ) then
    p(3)=1     ! disable high order treatment for z direction only
   endif
   
   !--------------------------------------------------------------------
   ! derived quantities:
   !--------------------------------------------------------------------
    nu(1) = delta(1)*( Lx/(pi*nx) )**(2*p(1)) / T_diss(1)
    nu(2) = delta(2)*( Ly/(pi*ny) )**(2*p(2)) / T_diss(2)
    if(high_order_z) then
     nu(3) = delta(3)*( Lz/(pi*nz) )**(2*p(3)) / T_diss(3)
    else
     nu(3) = viscosity
     T_diss(3) = 1.d0
     delta(3) = nu(3)*T_diss(3)/(Lz/(pi*nz))**2  ! just so I can see what it is
    endif
   
   !--------------------------------------------------------------------
   ! corresponding dimensionless coefficients
   !--------------------------------------------------------------------
    do dim=1,2
     diff_coeffs(dim,1:3) = nu(dim)/(velocity_scale*(length_scale)**(2*p(dim)-1) )
    enddo
    
    dim=3
    if(high_order_z) then
     diff_coeffs(dim,1:3) = nu(dim)/(velocity_scale*(length_scale)**(2*p(dim)-1) )
    else
     diff_coeffs(dim,1:3) = 1.d0/Re
    endif
    
    diff_coeffs(:,4)=diff_coeffs(:,1)*(kappa(1)/viscosity)
    if( do_second_scalar )  &
     diff_coeffs(:,5)=diff_coeffs(:,1)*(kappa(2)/viscosity)
        
  endif   ! end of high_order_diffusion true or false block
      


99 continue  
  !--------------------------------------------------------------------
  !  Helmholtz terms needed for implicit solves:
  !     q0(x,y,fid) = -alpha*(kx)^(2p) - alpha*(ky)^2p
  !                   alpha,p chosen for each direction/ variable
  !
  !     q1(x,y) = -(kx)^2 - (ky)^2   (used for pressure poisson eqn)
  !
  !     k-collapsed ZBLOCK format  (i,j) <--> (x,y)  
  !                               (ZBLOCK is (z,x,y))
  !--------------------------------------------------------------------
    do i=1,array_size(JDIM,ZBLOCK,myid)
     i_global = global_x_indices(START,ZBLOCK,myid) + i - 1
     do j=1,array_size(KDIM,ZBLOCK,myid)
      j_global = global_y_indices(START,ZBLOCK,myid) + j - 1
      
      q1(i,j) = - (kx(i_global))**2     &
                - (ky(j_global))**2
      
      !--------------------------------------------------------------------
      ! set q1 to special value for nyquist wavenumbers in x and y
      !--------------------------------------------------------------------
      if( deriv_type(6,1,1) == 'fourier'  &
         .and. i_global==nx/2+1           &
         .and. nx>1 ) then  ! pressure, x, 1st deriv
       !q1(i,j) = 999.d0
      elseif( deriv_type(6,1,1) == 'sin'  &
         .or. deriv_type(6,1,1) == 'cos'  &
         .and. i_global==nx               &
         .and. nx>1 ) then
        !q1(i,j) = 999.d0
       endif
       
      if( deriv_type(6,2,1) == 'fourier'  &
         .and. j_global==ny/2+1           &
         .and. ny>1 ) then  ! pressure, y, 1st deriv
       !q1(i,j) = 999.d0
      elseif( deriv_type(6,2,1) == 'sin'  &
         .or. deriv_type(6,2,1) == 'cos'  &
         .and. j_global==ny               &
         .and. ny > 1 ) then
        !q1(i,j) = 999.d0
       endif
       
      
      do fid=1,npvs
       q0(i,j,fid) = - diff_coeffs(1,fid) * (kx(i_global))**(2.d0*p(1))     &
                     - diff_coeffs(2,fid) * (ky(j_global))**(2.d0*p(2))
      enddo
            
     enddo
    enddo
  
  
  
  !--------------------------------------------------------------------
  ! report:
  !--------------------------------------------------------------------
   if( myid==0) then
    do k=0,1   
     write(k,*)            ' ................      high_order_z ',high_order_z
     write(k,*)            ' ................      high_order_operators ',high_order_operators
     write(k,'(a,3e18.6)') ' ................      nu, kappa(2) [m2/s] ', viscosity,kappa(:)
     write(k,'(a,3e18.6)') ' ................      Nyquist scales  [m] ', (nx*pi/Lx),(ny*pi/Ly),(nz*pi/Lz)
     write(k,'(a,3i6)')    ' ................      order of diff operators ', 2*p(:)
     write(k,'(a,3e18.6)') ' ................      x decay time scale T_diss [s],[min],[hrs] ', T_diss(1),T_diss(1)/60.,T_diss(1)/3600.
     write(k,'(a,3e18.6)') ' ................      y decay time scale T_diss [s],[min],[hrs] ', T_diss(2),T_diss(2)/60.,T_diss(2)/3600.
     write(k,'(a,3e18.6)') ' ................      z decay time scale T_diss [s],[min],[hrs] ', T_diss(3),T_diss(3)/60.,T_diss(3)/3600.
     write(k,'(a,3e18.6)') ' ................      delta factors prescribed:  ',delta(:)
     write(k,'(a,1e18.6)') ' ................      decay factor exp(-delta) after T_diss at the Nyquist spatial scales:  '
     write(k,'(a,3e18.6)') ' ................                    ', exp(-delta(1)),exp(-delta(2)),exp(-delta(3))
     write(k,*) ' ................'
     write(k,*) ' -----> SetupDiffusion routine exiting normally  <---------- '
     write(k,*) ' ................'
    enddo    
    close(1)
   endif
 
 call mpi_barrier(comm,ierr)
return 
end subroutine SetupDiffusion
