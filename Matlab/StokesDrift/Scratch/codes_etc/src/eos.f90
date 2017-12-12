subroutine equation_of_state
!-----------------------------------------------------------------------
!     Routine to calculate d'less pd from values of d'less s1 and s2
!     and their reference profiles s1_bar(z) and s2_bar(z)
!     The "perturbation" density is defined relative to the reference
!     profile rho_bar(z) which is determined by s1_bar and s2_bar.
!
!     EOS from MILLERO et.al. (1980) DEEP-SEA RES.
!     The scalar s1=temparature in deg C and 
!     the scalar s2 is salinity in psu (i.e. ~ 36)
!-----------------------------------------------------------------------
 use dimensional_scales
 use decomposition_params
 use pde_params,             only: Ri
 use mpi_params,             only: myid
 use dependent_variables,    only: s1,s2,s1_bar,s2_bar,scalar_kind
 use independent_variables,  only: nx,ny,z
 use intermediate_variables, only: pd
 use methods_params,         only: do_second_scalar
 
 
 implicit none 
 integer                        :: i,j,k,k_global
 real(kind=8),save              :: xx,yy,PMpa
 real(kind=8)                   :: total_density
 real(kind=8)                   :: temperature,salinity
 real(kind=8)                   :: scale_factor
 real(kind=8),allocatable,save  :: rho_bar(:)
 logical,save                   :: first_entry=.TRUE.
  
 if( first_entry ) then
  PMpa = 0.0   !   compute potential density, not in situ density

  if( trim(scalar_kind(1)) .eq. 't'                    &
     .and. trim(scalar_kind(2)) .eq. 's' ) then
           
            xx=1.0    ! use this to keep  Salinity
            yy=1.0    ! use this to keep Temperature
            
  elseif(trim(scalar_kind(1)) .eq. 't'              &
        .and. trim(scalar_kind(2)) .ne. 's' ) then
           
            yy=1.0    ! use this to keep Temperature
            xx=0.0    ! use this to neglect Salinity
            
  elseif(trim(scalar_kind(1)) .ne. 't'              &
        .and. trim(scalar_kind(2)) .eq. 's' ) then
           
            xx=1.0    ! use this to keep  Salinity
            yy=0.0    ! use this to neglect Temperature
                                
  endif
  xx = xx/1000.d0  ! for units consistency w/ fortran routine rho.f
  
  if( do_second_scalar .and. Ri /= 0.d0 ) then
   allocate( rho_bar(array_size(KDIM,YBLOCK,myid)) )
    do k=1,array_size(KDIM,YBLOCK,myid)  
     k_global = global_z_indices(START,YBLOCK,myid) + k - 1 
     Temperature = yy*(s1_bar(k,1)*scalar_scale(1))   !! [deg C]
     Salinity =    xx*(s2_bar(k,1)*scalar_scale(2))   !! [psu]
     call rho(Salinity,Temperature,Pmpa,rho_bar(k))   !! ambient density [kg/m3]
    enddo
  endif
  
  first_entry=.FALSE.
 endif
 
 if( Ri==0.d0 ) then
  return    ! don't need to access pd            
 endif
 

 if( trim(scalar_kind(1)) .eq. 'r' ) then      
  call copy(s1,pd)   ! pd <-- s1
  return  
 endif
      
 if( trim(scalar_kind(2)) .eq. 'r' ) then     
  call copy(s2,pd)   ! pd <-- s2
  return                        
 endif

 if( trim(scalar_kind(1)) .eq. 'r'  .AND. trim(scalar_kind(2)) .eq. 'r') then
  stop 'both scalars set to r=density!!'
 endif

if( do_second_scalar ) then
 scale_factor = 1.d0/density_scale     
 do k=1,array_size(KDIM,YBLOCK,myid) 
 
  k_global = global_z_indices(START,YBLOCK,myid) + k - 1 
     
  do j=1,array_size(JDIM,YBLOCK,myid)    ! 2nd index is x in YBLOCK
   do i=1,array_size(IDIM,YBLOCK,myid)   ! 1st index is y in YBLOCK
           
    !!now repeat the calculations for total s1 & s2(x,y,z)
    Temperature = yy*( (s1_bar(k_global,1)  +  &
                        s1(i,j,k))*scalar_scale(1) )  ! dimensional, total Temp.
 
    Salinity    = xx*( (s2_bar(k_global,1)  +  &
                        s2(i,j,k))*scalar_scale(2) )  ! dimensional, total Salinity                        

    call rho(Salinity,Temperature,PMpa,total_density) ! total density [kg/m3]
           
    pd(i,j,k) = (total_density - rho_bar(k))*scale_factor
    
   enddo
  enddo
 enddo  
endif

end subroutine equation_of_state
 
 
 
 
 
 
 
!-----------------------------------------------------------------------  
!   "rho.f"

!   Parameter   :   equation of state
!   Source      :   MILLERO et.al. (1980) DEEP-SEA RES.
!                   Vol 27A, pp 255-264
!   Units       :   Kg / m**3
!   Input       :   s (concentration units), t (deg C), p (MPa)
!   Range       :   s (0 to 0.040), t (-2 to 40), p (MPa)
!   Ck value    :   K (0, 0, 100)      = 2.297721e4
!                   K (0.035, 0, 100)  = 2.499200e4
!                   K (0, 25, 100)     = 2.54051e4
!                   K (0.035, 25, 100) = 2.710895e4
!   Coded       :   Ngoc Dang, July 1982

subroutine rho (s, t, p, density)
 real(kind=8) s, t, p
 real(kind=8) a, b, c, d, e, k, t1, t2, t3, t4, ta, tb
 real(kind=8) k0, k0w, alpha0, alpha, rho0, aw, bw, density 


!  convert input units
 s = 1000.0d0 * s
 p = 10.0d0 * p

!  compute alpha0 = alpha (s, t, 0)
 t1 = 9.99841594d-1 + t * (6.793952d-5 + t * (-9.095290d-6 + t * &
 (1.001685d-7 + t * (-1.120083d-9 + t * 6.536332d-12))))
     
 t2 = 8.25917d-4 + t * (-4.4490d-6 + t * (1.0485d-7 + t *    &   
 (-1.258d-9 + t * 3.315d-12)))
             
 t3 = -6.33761d-6 + t * (2.8441d-7 + t * (-1.6871d-8 + t *   &   
 2.83258d-10))
             
 t4 = 5.4705d-7 + t * (-1.97975d-8 + t * (1.6641d-9 - t *    &   
 3.1203d-11))
     
 rho0 = t1 + t2 * s + t3 * s**(1.5) + t4 * s**2
 alpha0 = 1.0d0 / rho0

!  compute the pressure coefficient

!   the first term
 k0w = 1.965221d4 + t * (1.484206d2 + t * (-2.327105d0 + t *   &  
 (1.360477d-2 - t * 5.155288d-5)))
              
 ta = 5.46746d1 + t * (-6.03459d-1 + t * (1.09987d-2 - t *     &    
 6.167d-5))
             
 tb = 7.944d-2 + t * (1.6483d-2 - t * 5.3009d-4)
 k0 = k0w + ta * s + tb * s**(1.5)

!   the second term
 aw = 3.239908d0 + t * (1.43713d-3 + t * (1.16092d-4 - t *  &     
 5.77905d-7))
     
 c = 2.2838d-3 + t * (-1.0981d-5 - t * 1.6078d-6)
 d = 1.9107d-4
 a = aw + c * s + d * s**(1.5)

!   the third term
 bw = 8.50935d-5 + t * (-6.12293d-6 + t * 5.2787d-8)
 e = -9.9348d-7 + t * (2.0816d-8 + t * 9.1697d-10)
 b = bw + e * s

 k = k0 + a * p + b * p**2

!  compute density
 alpha = alpha0 * (1.0d0 - p / k)
 density = 1.0d0 / alpha

!  convert to kg / m**3
 density = 1000.0d0 * density
 s = s / 1000.0d0
 p = p / 10.0d0

return
end
