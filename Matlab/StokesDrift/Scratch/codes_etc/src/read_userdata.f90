subroutine ReadUserData
 use etc
 use dimensional_scales
 use pde_params
 use diffusion_params
 use methods_params
 use mpi_params
 use decomposition_params,   only: np,p1,p2 
 use io_params  
 use dependent_variables
 use independent_variables
 use intermediate_variables
 use particles
 implicit none  

 include 'mpif.h' 
 integer           :: processor_doing_reads
 integer           :: id
 integer           :: endpt
 integer           :: i
 integer           :: j
 integer           :: idim
 complex(kind=8)   :: dummy_beta_complex
 character(len=80) :: file1,file2,file3,file4,file5

 if(myid==0) write(0,*) ' ................'
 if(myid==0) write(0,*) ' ................     hello world from read_user_data'
 

 file1='input/problem_params'
 file2='input/io_params'
 file3='input/BoundaryConditions'
 file4='input/DifferentiationMethods'
 file5='input/high_order_diffusion_params'

  
  if(myid==0) open(1,file=file1,position='rewind') 
  
  if(myid==0) read(1,*) runlabel
  call mpi_bcast(runlabel,80,MPI_CHARACTER,0,comm,ierr)
     
  if(myid==0) read(1,*) do_nonlinear
  call mpi_bcast(do_nonlinear,1,MPI_LOGICAL,0,comm,ierr)
     
  if(myid==0) read(1,*) do_second_scalar
  call mpi_bcast(do_second_scalar,1,MPI_LOGICAL,0,comm,ierr)
     
  if(myid==0) read(1,*) vertical_coriolis
  call mpi_bcast(vertical_coriolis,1,MPI_LOGICAL,0,comm,ierr)
     
  if(myid==0) read(1,*) p1
  call mpi_bcast(p1,1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) p2
  call mpi_bcast(p2,1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) nx
  call mpi_bcast(nx,1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) ny
  call mpi_bcast(ny,1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) nz
  call mpi_bcast(nz,1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) dt
  call mpi_bcast(dt,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) t0
  call mpi_bcast(t0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) tf
  call mpi_bcast(tf,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) Lx
  call mpi_bcast(Lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) Ly
  call mpi_bcast(Ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) Lz
  call mpi_bcast(Lz,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) scalar_kind(1)
  call mpi_bcast(scalar_kind(1),1,MPI_CHARACTER,0,comm,ierr)
  
  if(myid==0) read(1,*) scalar_kind(2)
  call mpi_bcast(scalar_kind(2),1,MPI_CHARACTER,0,comm,ierr)
  
  if(myid==0) read(1,*) do_forcing
  call mpi_bcast(do_forcing,1,MPI_LOGICAL,0,comm,ierr)
  
  if(myid==0) read(1,*) do_immersed_boundary
  call mpi_bcast(do_immersed_boundary,1,MPI_LOGICAL,0,comm,ierr)
  
  if(myid==0) read(1,*) g
  call mpi_bcast(g,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) f(1)          !! f
  call mpi_bcast(f(1),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) then   !  read(1,*) f(2)          !! beta
   read(1,*) dummy_beta_complex
   f(2) = real( dummy_beta_complex )            !! beta
   if( aimag( dummy_beta_complex ) == 0 ) then
    y_pivot = Ly/2.d0                           !! default pivot location
   else
    y_pivot = aimag( dummy_beta_complex )*Ly    !! fraction*Ly, pivot for beta plane
   endif
  endif
  call mpi_bcast(f(2),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  call mpi_bcast(y_pivot,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) rho_0
  call mpi_bcast(rho_0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) nu
  call mpi_bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) kappa(1)
  call mpi_bcast(kappa(1),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) kappa(2)
  call mpi_bcast(kappa(2),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) high_order_operators
  call mpi_bcast(high_order_operators,1,MPI_LOGICAL,0,comm,ierr)
  
  if(myid==0) read(1,*) dgrad
  call mpi_bcast(dgrad,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) scalar_scale(1)
  call mpi_bcast(scalar_scale(1),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) scalar_scale(2)
  call mpi_bcast(scalar_scale(2),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) u0
  call mpi_bcast(u0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  if(myid==0) read(1,*) nparticles
  call mpi_bcast(nparticles,1,MPI_INTEGER,0,comm,ierr)
  
  if(nparticles>0) then
   if(myid==0) read(1,*) particle_write_inc
   call mpi_bcast(particle_write_inc,1,MPI_INTEGER,0,comm,ierr)
  endif
  
 if(myid==0) close(1)
 if(np .ne. p1*p2)       stop ' processor grid error, p1*p2 NE np '     
 if(numprocs .ne. np)    stop ' processor grid error, numprocs NE np '    
 if(myid==0) write(0,*) ' ................      ',trim(file1),' read'

 
 if(myid==0) open(1,file=file2,position='rewind')
 if(myid==0) read(1,*) num_file_sets
 call mpi_bcast(num_file_sets,1,MPI_INTEGER,0,comm,ierr)
 
 
 do i=1,num_file_sets
  if(myid==0) read(1,*) filename_root(i)
  call mpi_bcast(filename_root(i),80,MPI_CHARACTER,0,comm,ierr)
  
  if(myid==0) read(1,*) mode(i)
  call mpi_bcast(mode(i),80,MPI_CHARACTER,0,comm,ierr)
  
  if(myid==0) read(1,*) nsteps(i)
  call mpi_bcast(nsteps(i),1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) ilocs(1,i)
  call mpi_bcast(ilocs(1,i),1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) ilocs(2,i)
  call mpi_bcast(ilocs(2,i),1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) ilocs(3,i)
  call mpi_bcast(ilocs(3,i),1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) jlocs(1,i)
  call mpi_bcast(jlocs(1,i),1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) jlocs(2,i)
  call mpi_bcast(jlocs(2,i),1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) jlocs(3,i)
  call mpi_bcast(jlocs(3,i),1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) klocs(1,i)
  call mpi_bcast(klocs(1,i),1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) klocs(2,i)
  call mpi_bcast(klocs(2,i),1,MPI_INTEGER,0,comm,ierr)
  
  if(myid==0) read(1,*) klocs(3,i)
  call mpi_bcast(klocs(3,i),1,MPI_INTEGER,0,comm,ierr)
  
  do j=1,5 !! u,v,w,s1,s2
   if(myid==0) read(1,*) variable_key(j,i)
   call mpi_bcast(variable_key(j,i),1,MPI_LOGICAL,0,comm,ierr)
  enddo
  
  if( variable_key(4,i) == 1 ) write_s1_bar(i)=.FALSE.
  if( variable_key(4,i) == 2 ) write_s1_bar(i)=.TRUE.
  
  if( variable_key(5,i) == 1 ) write_s2_bar(i)=.FALSE.
  if( variable_key(5,i) == 2 ) write_s2_bar(i)=.TRUE.
  
  if(myid==0) write(0,*) ' ................        output fileset number ',i,' read'
  
  !! set these to 1 during testing
  variable_key(6,i)=0   ! divustar
  
  !! set permanently to 1
  variable_key(7,i)=1   ! phi
  
  !! not pd though
  variable_key(8,i)=0   ! pd
  
  !! quick sanity check ....
  if(ilocs(2,i) > nx ) stop 'io_params data > nx'
  if(jlocs(2,i) > ny ) stop 'io_params data > ny'
  if(klocs(2,i) > nz ) stop 'io_params data > nz'
 enddo 
 if(myid==0) write(0,*) ' ................      ',trim(file2),' read'
 if(myid==0) close(1)



 if(myid==0) open(1,file=file3,position='rewind')
 do id=1,6  !! including pressure
  do idim=1,3
   do endpt=1,2
    if(myid==0) read(1,*) bc_type(id,idim,endpt)
    call mpi_bcast(bc_type(id,idim,endpt),80,MPI_CHARACTER,0,comm,ierr)
   enddo
  enddo
 enddo
 if(myid==0) write(0,*) ' ................      ',trim(file3),' read'
if(myid==0)  close(1)


 
 if(myid==0) open(1,file=file4,position='rewind')
  if(myid==0) read(1,*) deriv_type(1,1,1)   ! u,x,1st deriv
     call mpi_bcast(deriv_type(1,1,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(1,2,1)   ! u,y,1st deriv
     call mpi_bcast(deriv_type(1,2,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(1,3,1)   ! u,z,1st deriv
     call mpi_bcast(deriv_type(1,3,1),80,MPI_CHARACTER,0,comm,ierr)

  if(myid==0) read(1,*) deriv_type(2,1,1)   ! v,x,1st deriv
     call mpi_bcast(deriv_type(2,1,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(2,2,1)   ! v,y,1st deriv
     call mpi_bcast(deriv_type(2,2,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(2,3,1)   ! v,z,1st deriv
     call mpi_bcast(deriv_type(2,3,1),80,MPI_CHARACTER,0,comm,ierr)

  if(myid==0) read(1,*) deriv_type(3,1,1)   ! w,x,1st deriv
     call mpi_bcast(deriv_type(3,1,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(3,2,1)   ! w,y,1st deriv
     call mpi_bcast(deriv_type(3,2,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(3,3,1)   ! w,z,1st deriv
     call mpi_bcast(deriv_type(3,3,1),80,MPI_CHARACTER,0,comm,ierr)

  if(myid==0) read(1,*) deriv_type(4,1,1)   ! s1,x,1st deriv
     call mpi_bcast(deriv_type(4,1,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(4,2,1)   ! s1,y,1st deriv
     call mpi_bcast(deriv_type(4,2,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(4,3,1)   ! s1,z,1st deriv
     call mpi_bcast(deriv_type(4,3,1),80,MPI_CHARACTER,0,comm,ierr)

  if(myid==0) read(1,*) deriv_type(5,1,1)   ! s2,x,1st deriv
     call mpi_bcast(deriv_type(5,1,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(5,2,1)   ! s2,y,1st deriv
     call mpi_bcast(deriv_type(5,2,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(5,3,1)   ! s2,z,1st deriv
     call mpi_bcast(deriv_type(5,3,1),80,MPI_CHARACTER,0,comm,ierr)

  if(myid==0) read(1,*) deriv_type(6,1,1)   ! phi,x,1st deriv
     call mpi_bcast(deriv_type(6,1,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(6,2,1)   ! phi,y,1st deriv
     call mpi_bcast(deriv_type(6,2,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) read(1,*) deriv_type(6,3,1)   ! phi,z,1st deriv
     call mpi_bcast(deriv_type(6,3,1),80,MPI_CHARACTER,0,comm,ierr)
  if(myid==0) write(0,*) ' ................      ',trim(file4),' read'
 if(myid==0) close(1)

 if( high_order_operators ) then
  if(myid==0) open(1,file=file5,position='rewind')
   if(myid==0) read(1,*) p(1)
   call mpi_bcast(p(1),1,MPI_INTEGER,0,comm,ierr)
   
   if(myid==0) read(1,*) p(2)
   call mpi_bcast(p(2),1,MPI_INTEGER,0,comm,ierr)
   
   if(myid==0) read(1,*) p(3)
   call mpi_bcast(p(3),1,MPI_INTEGER,0,comm,ierr)
   
   if(myid==0) read(1,*) T_diss(1)
   call mpi_bcast(T_diss(1),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
   
   if(myid==0) read(1,*) T_diss(2)
   call mpi_bcast(T_diss(2),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
   
   if(myid==0) read(1,*) T_diss(3)
   call mpi_bcast(T_diss(3),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
   
   if(myid==0) read(1,*) delta(1)
   call mpi_bcast(delta(1),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
   
   if(myid==0) read(1,*) delta(2)
   call mpi_bcast(delta(2),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
   
   if(myid==0) read(1,*) delta(3)
   call mpi_bcast(delta(3),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(myid==0) close(1)
 endif

 
 
 idim=1
 if( bc_type(1,idim,1) == 'periodic' ) then
  bc_type(:,idim,:)='periodic'
  xdim_periodic=.TRUE.
 else
  xdim_periodic=.FALSE.
 endif
 
 
 idim=2
 if( bc_type(1,idim,1) == 'periodic' ) then
  bc_type(:,idim,:)='periodic'
  ydim_periodic=.TRUE.
 else
  ydim_periodic=.FALSE.
 endif
 
 
 idim=3
 if( bc_type(1,idim,1) == 'periodic' ) then
  bc_type(:,idim,:)='periodic'
  zdim_periodic=.TRUE.
 else
  zdim_periodic=.FALSE.
 endif
 

 
if(myid==0) write(0,*) ' ................     '


 call mpi_barrier(comm,ierr)
 return
end subroutine ReadUserData
