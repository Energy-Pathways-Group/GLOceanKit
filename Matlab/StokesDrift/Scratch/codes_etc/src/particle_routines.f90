!=============================================================== 
!===============================================================
  subroutine advance_particles
   use particles
   use mpi_params
   use dependent_variables,     only: u,v,w
   use etc,                     only: PM0,P_oldest

   use interpolation
   
   
   
   implicit none
   include 'mpif.h' 
   real(kind=8),dimension(:),pointer     :: tmp
   integer                               :: i
   

   if(nparticles<=0) return
   
   call wrap_particles       ! if necessary due to periodicity
   call mark_my_particles    ! if "i" owned by myid, my_labels(i)<=1
   !call find_s_values        ! convert all particles x/y/z to s vals
   
   !-------------------------------------------------------
   ! do the interpolation step: i.e. 
   ! given u, 
   ! calculate the e3 array 
   !   ==> e3 stored in tmpY(:,:,:,1)
   !-------------------------------------------------------
   call Spline_Interp_YB(u)
   
   
   !-------------------------------------------------------
   ! do the evaluation step: 
   !  use e3 and the s position values of particles to
   !  compute the u values at the particle positions
   ! ===> for particles owned by some other processor,
   !      the returned value is ZERO <===
   !-------------------------------------------------------
   call Spline_Eval_YB(uvels(1,PM0))
   
   !----------------------
   ! repeat for v and w:
   !----------------------
   call Spline_Interp_YB(v)
   call Spline_Eval_YB(vvels(1,PM0))
   
   call Spline_Interp_YB(w)
   call Spline_Eval_YB(wvels(1,PM0))
   
  !-------------------------------------------------------
  ! At this point positions(i,:) are nonzero only if
  ! particle i resides on processor myid. Similarly
  ! for uvels(i,PM0),vvels(i,PM0) etc. Values for
  ! PM1, PM2 etc should be identical for all processors.
  !-------------------------------------------------------
  
  
  !-------------------------------------------------------
  ! Given the velocities, time step the particle positions
  !-------------------------------------------------------
  call particle_step
  
 !------------------------------------------------------- 
 ! Now globally_update_particle_data
 ! so that all processors have the data for each
 ! particle (if one crosses into "MY" area, I need
 ! to know its current position and previous velocities)
 !-------------------------------------------------------
   tmp => uvels(:,P_oldest)   ! this is free space after particle_step
   
   tmp(:) = uvels(:,PM0)
   call mpi_allreduce(tmp,uvels(1,PM0),nparticles,                 &
                      mpi_double_precision,mpi_sum,comm,ierr )
   
   tmp(:) = vvels(:,PM0)
   call mpi_allreduce(tmp,vvels(1,PM0),nparticles,                 &
                      mpi_double_precision,mpi_sum,comm,ierr )
   
   tmp(:) = wvels(:,PM0)
   call mpi_allreduce(tmp,wvels(1,PM0),nparticles,                 &
                      mpi_double_precision,mpi_sum,comm,ierr )
   
   do i=1,3
    tmp(:) = positions(:,i)
    call mpi_allreduce(tmp,positions(1,i),nparticles,                 &
                       mpi_double_precision,mpi_sum,comm,ierr )
   enddo
   
 !------------------------------------------------------- 
 ! All processors now have updated positions and velocities
 ! for all particles, whether they reside locally or not
 !-------------------------------------------------------
end subroutine advance_particles


subroutine wrap_particles
!--------------------------------------
! periodic wrapping where necessary
!--------------------------------------
 use independent_variables, only: x,y,z,nx,ny,nz,  &
                                  xdim_periodic,   &
                                  ydim_periodic,   &
                                  zdim_periodic,   &
                                  Lx,Ly,Lz
 use dimensional_scales,    only: length_scale
 use particles,             only: positions,nparticles
 
 implicit none
 integer                       :: i,idim
 real(kind=8)                  :: hb2,shift
 real(kind=8)                  :: x0,x1,y0,y1,z0,z1
 
 !-----------------------------------------------------------
 !  keep particle positions within, for example,
 !   x(1) - dx/2   ( slightly < 0 )
 !      and
 !   x(n) + dx/2 = Lx - dx/2  ( slightly less than Lx )
 !   so that they are always within the extrapolation range
 !   of the B-spline interpolation routines
 !-----------------------------------------------------------
 
 
 if(xdim_periodic .and. nx > 1) then
  idim=1
  shift = Lx/length_scale
  hb2 = 0.5*( x(2)-x(1) )
  x0 = -hb2
  x1 = (Lx/length_scale) - hb2
  do i=1,nparticles
   if( positions(i,idim) < x0 ) then
    positions(i,idim) = positions(i,idim) + shift
   elseif( positions(i,idim) > x1 ) then
    positions(i,idim) = positions(i,idim) - shift
   endif
  enddo
 endif
  
 if(ydim_periodic .and. ny > 1) then
  idim=2
  shift = Ly/length_scale
  hb2 = 0.5*( y(2)-y(1) )
  y0 = -hb2
  y1 = (Ly/length_scale) - hb2
  do i=1,nparticles
   if( positions(i,idim) < y0 ) then
    positions(i,idim) = positions(i,idim) + shift
   elseif( positions(i,idim) > y1 ) then
    positions(i,idim) = positions(i,idim) - shift
   endif
  enddo
 endif
 
 if(zdim_periodic .and. nz > 1) then
  idim=3
  shift = Lz/length_scale
  hb2 = 0.5*( z(2)-z(1) )
  z0 = -hb2
  z1 = (Lz/length_scale) - hb2
  do i=1,nparticles
   if( positions(i,idim) < z0 ) then
    positions(i,idim) = positions(i,idim) + shift
   elseif( positions(i,idim) > z1 ) then
    positions(i,idim) = positions(i,idim) - shift
   endif
  enddo
 endif
 
 
end subroutine wrap_particles



subroutine mark_my_particles
!--------------------------------------
! id particles in my extended YBLOCK 
!--------------------------------------
 use mpi_params,            only: myid
 use independent_variables, only: x,z,nx,nz
 use decomposition_params
 use particles,             only: positions,     &
                                  my_labels,     &
                                  nparticles
                                   
 implicit none
 integer                       :: i
 integer,save                  :: i0,i1
 real(kind=8),save             :: dx,dz
 real(kind=8),save             :: x0,x1,z0,z1
 logical,save                  :: first_entry=.TRUE.
 
 if( first_entry ) then
 
  i0 = global_x_indices(START,YBLOCK,myid)
  if( i0 == 1 .and. nx > 1 ) then
   dx = x(i0+1)-x(i0)
  elseif( i0 > 1 .and. nx > 1 ) then
   dx = x(i0)-x(i0-1)
  elseif( nx==1 ) then
   dx = 0.d0
  endif
  x0 = x(i0) - dx/2.d0
  
  i1 = global_x_indices(END,YBLOCK,myid)
  if( i1 == nx .and. nx > 1 ) then
   dx = x(i1)-x(i1-1)
  elseif( i1 < nx .and. nx > 1 ) then
   dx = x(i1+1)-x(i1)
  elseif( nx==1 ) then
   dx = 0.d0
  endif
  x1 = x(i1) + dx/2.d0
  
  i0 = global_z_indices(START,YBLOCK,myid)
  if( i0 == 1 .and. nz > 1 ) then
   dz = z(i0+1)-z(i0)
  elseif( i0 > 1 .and. nz > 1 ) then
   dz = z(i0)-z(i0-1)
  elseif( nz==1 ) then
   dz = 0.d0
  endif
  z0 = z(i0) - dz/2.d0
  
  i1 = global_z_indices(END,YBLOCK,myid)
  if( i1 == nz .and. nz > 1 ) then
   dz = z(i1)-z(i1-1)
  elseif( i1 < nz .and. nz > 1 ) then
   dz = z(i1+1)-z(i1)
  elseif( nz==1 ) then
   dz = 0.d0
  endif
  z1 = z(i1) + dz/2.d0
    
  first_entry=.FALSE.
 endif
 
 
 do i=1,nparticles
 
  if( positions(i,3) >= z0 .and.      &
      positions(i,3) <  z1 )  then
      
      my_labels(i) = 1   !  particle in my z range
      
      
      if( nx>1 .and. positions(i,1) < x0 ) then    ! to left of my x range
       my_labels(i) = 0
       positions(i,1:3) = 0.d0
      endif
      
      if( nx>1 .and. positions(i,1) >= x1 ) then   ! to right of my x range
       my_labels(i) = 0
       positions(i,1:3) = 0.d0
      endif
                 
   else
       
    my_labels(i) = 0
    positions(i,1:3) = 0.d0
          
   endif
  enddo
    
end subroutine mark_my_particles


subroutine particle_step
 use particles
 use independent_variables,   only: dt
 use etc,                     only: PM0,PM1,PM2,PM3
 !---------------------------------------------------
 ! At this point positions(i,:) are nonzero only if
 ! particle i resides on processor myid. Similarly
 ! for uvels(i,PM0),vvels(i,PM0) etc. Values for
 ! PM1, PM2 etc should be identical for all processors.
 !---------------------------------------------------
  
  implicit none
  integer                  :: i
  real(kind=8), save       :: dt2,dt12,dt24
  logical,      save       :: first_entry = .TRUE.
   
  if( first_entry ) then
   dt2=dt/2.d0
   dt12=dt/12.d0
   dt24=dt/24.d0
   first_entry = .FALSE.
  endif
  
  !---------------------------------------------------
  ! Integrate forward one time step. 
  !---------------------------------------------------
  do i=1,nparticles
  
  if( my_labels(i)==1 ) then  !! only deal with those on proc. myid
        
   if( particle_step_flag .eq. 'euler' ) then
    positions(i,1)=positions(i,1) + dt*uvels(i,PM0)    
    positions(i,2)=positions(i,2) + dt*vvels(i,PM0)
    positions(i,3)=positions(i,3) + dt*wvels(i,PM0)
    particle_step_flag='AB2'    
    
   elseif( particle_step_flag .eq. 'AB2' ) then
     positions(i,1)=positions(i,1) + dt2*( 3.*uvels(i,PM0)       &
                                            - uvels(i,PM1) )
                                            
     positions(i,2)=positions(i,2) + dt2*( 3.*vvels(i,PM0)       &
                                            - vvels(i,PM1) )
                                            
     positions(i,3)=positions(i,3) + dt2*( 3.*wvels(i,PM0)       &
                                            - wvels(i,PM1) )
     particle_step_flag='AB3'
    
   elseif( particle_step_flag .eq. 'AB3' ) then   
     positions(i,1)=positions(i,1) + dt12*(23.*uvels(i,PM0)        &
                                         - 16.*uvels(i,PM1)        &
                                         +  5.*uvels(i,PM2) )
                                         
     positions(i,2)=positions(i,2) + dt12*(23.*vvels(i,PM0)        &
                                         - 16.*vvels(i,PM1)        &
                                         +  5.*vvels(i,PM2) ) 
                                         
     positions(i,3)=positions(i,3) + dt12*(23.*wvels(i,PM0)        &
                                         - 16.*wvels(i,PM1)        &
                                         +  5.*wvels(i,PM2) )
     particle_step_flag='AB4'
   
   elseif( particle_step_flag .eq. 'AB4' ) then    
     positions(i,1)=positions(i,1) + dt24*(55.*uvels(i,PM0)        &
                                         - 59.*uvels(i,PM1)        & 
                                         + 37.*uvels(i,PM2)        & 
                                         -  9.*uvels(i,PM3) )
                                         
     positions(i,2)=positions(i,2) + dt24*(55.*vvels(i,PM0)        & 
                                         - 59.*vvels(i,PM1)        & 
                                         + 37.*vvels(i,PM2)        & 
                                         -  9.*vvels(i,PM3) )
                                         
     positions(i,3)=positions(i,3) + dt24*(55.*wvels(i,PM0)        & 
                                         - 59.*wvels(i,PM1)        & 
                                         + 37.*wvels(i,PM2)        & 
                                         -  9.*wvels(i,PM3) )
   endif
   
  endif
   
 enddo
end subroutine particle_step


