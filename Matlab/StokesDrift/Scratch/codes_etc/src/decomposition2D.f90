subroutine decomposition2d
  use mpi_params,               only: myid,comm,ierr
  use decomposition_params
  use independent_variables,    only: nx,ny,nz
  
  implicit none
  integer                          :: pid,i,j,k,npts,ii
  
  if(myid==0) then
   write(0,*) ' ................'
   write(0,*) ' ................     hello world from decomposition2d'
   write(0,*) ' ................      2D Processor Grid:'
   write(0,*) ' ................',np,'=',p1,'x',p2
   write(0,*) ' ................'
  endif

  !-------------------
  ! nx is split by p1
  !-------------------
  if( nx < p1 ) then
   if(myid==0) write(0,*) ' decomposition problem: nx < p1 ',nx,p1
   stop
  endif
  
  !-------------------
  ! nz is split by p2
  !-------------------
  if( nz < p2 ) then
   if(myid==0) write(0,*) ' decomposition problem: nz < p2 ',nz,p2
   stop
  endif
  
  !-------------------------------
  ! ny is split by both p1 and p2
  !-------------------------------
  if( ny < p1 ) then
   if(myid==0) write(0,*) ' decomposition problem: ny < p1 ',ny,p1
   stop
  endif
  if( ny < p2 ) then
   if(myid==0) write(0,*) ' decomposition problem: ny < p2 ',ny,p2
   stop
  endif


  allocate( proc_row(3,0:np-1) )             ! (?BLOCK, pid)
  allocate( proc_col(3,0:np-1) )             ! (?BLOCK, pid)
  allocate( array_size(3,3,0:np-1) )         ! (array index, ?BLOCK, pid )
  allocate( global_x_indices(2,3,0:np-1) )   ! (START/END , ?BLOCK, pid )
  allocate( global_y_indices(2,3,0:np-1) )   ! (START/END , ?BLOCK, pid )
  allocate( global_z_indices(2,3,0:np-1) )   ! (START/END , ?BLOCK, pid )
   
 
  !--------------------------------------------------------------------------
  ! define decomposition of coordinates across processors
  !--------------------------------------------------------------------------
  layout(:,XBLOCK) = (/1,p1,p2/)   ! y coord split by p1, zcoord split by p2
  layout(:,YBLOCK) = (/p1,1,p2/)   ! x coord split by p1, zcoord split by p2
  layout(:,ZBLOCK) = (/p1,p2,1/)   ! x coord split by p1, ycoord split by p2
  
  !--------------------------------------------------------------------------
  ! map spatial coordinates to storage indices for each decomposition
  !--------------------------------------------------------------------------
  ! mem_order(1) = the fortran array index corresponding to x
  ! mem_order(2) = the fortran array index corresponding to y
  ! mem_order(3) = the fortran array index corresponding to z 
  !-----------------------------------------------------------
  mem_order(:,XBLOCK) = (/1,2,3/)    ! x is first index,  y is second, z third  ==> (x,y,z)
  mem_order(:,YBLOCK) = (/2,1,3/)    ! x is second index, y is first,  z third  ==> (y,x,z)
  mem_order(:,ZBLOCK) = (/2,3,1/)    ! x is second index, y is third,  z first  ==> (z,x,y)
  
  
  if(myid==0) open(1,file='output/decomposition_details',position='rewind')
  
  !--------------------------------------------------------------------------
  !! XBLOCK definitions
  !--------------------------------------------------------------------------
  if(myid==0) write(1,*) '-------------------------------------'
  if(myid==0) write(1,*) 'XBLOCK: '
  if(myid==0) write(1,*) '-------------------------------------'
  if(myid==0) write(1,*) '         pid        row          col '
  if(myid==0) write(1,*) '-------------------------------------'
  do pid=0,np-1
   proc_row(XBLOCK,pid) = mod(pid,p1)
   proc_col(XBLOCK,pid) = floor(pid/float(p1))
   if(myid==0) write(1,*) pid,proc_row(XBLOCK,pid),proc_col(XBLOCK,pid)
  enddo
  if(myid==0) write(1,*) ' '
  
  if(myid==0) write(1,*) '----------------------------------------------------------------'
  if(myid==0) write(1,*) ' local array sizes/indices '
  if(myid==0) write(1,*) '         n elements      start index       end index      coord'
  if(myid==0) write(1,*) '----------------------------------------------------------------'
  do pid=0,np-1  
   
   array_size(IDIM,XBLOCK,pid) = nx
   global_x_indices(START,XBLOCK,pid)=1
   global_x_indices(END,XBLOCK,pid) =  global_x_indices(START,XBLOCK,pid)   &
                                    +  array_size(IDIM,XBLOCK,pid)-1
   
   if( proc_row(XBLOCK,pid) < p1-1 ) then
    array_size(JDIM,XBLOCK,pid) = floor(ny/float(p1))
   else
    array_size(JDIM,XBLOCK,pid) = floor(ny/float(p1)) + mod(ny,p1)
   endif
   global_y_indices(START,XBLOCK,pid)= proc_row(XBLOCK,pid)*floor(ny/float(p1)) + 1
   global_y_indices(END,XBLOCK,pid) =  global_y_indices(START,XBLOCK,pid)   &
                                    +  array_size(JDIM,XBLOCK,pid)-1
   
   if( proc_col(XBLOCK,pid) < p2-1 ) then
    array_size(KDIM,XBLOCK,pid) = floor(nz/float(p2))
   else
    array_size(KDIM,XBLOCK,pid) = floor(nz/float(p2)) + mod(nz,p2)
   endif
   global_z_indices(START,XBLOCK,pid)= proc_col(XBLOCK,pid)*floor(nz/float(p2)) + 1
   global_z_indices(END,XBLOCK,pid) =  global_z_indices(START,XBLOCK,pid)   &
                                    +  array_size(KDIM,XBLOCK,pid)-1
   
   if(myid==0) write(1,*) pid
   if(myid==0) write(1,*) '  1st index  ',array_size(IDIM,XBLOCK,pid),global_x_indices(:,XBLOCK,pid),'           x'
   if(myid==0) write(1,*) '  2nd index  ',array_size(JDIM,XBLOCK,pid),global_y_indices(:,XBLOCK,pid),'           y'
   if(myid==0) write(1,*) '  3rd index  ',array_size(KDIM,XBLOCK,pid),global_z_indices(:,XBLOCK,pid),'           z'
  enddo
  if(myid==0) write(1,*) ' '
  
  
  !--------------------------------------------------------------------------
  !! YBLOCK definitions
  !--------------------------------------------------------------------------
  if(myid==0) write(1,*) '-------------------------------------'
  if(myid==0) write(1,*) 'YBLOCK: '
  if(myid==0) write(1,*) '-------------------------------------'
  if(myid==0) write(1,*) '         pid        row          col '
  if(myid==0) write(1,*) '-------------------------------------'
  do pid=0,np-1
   proc_row(YBLOCK,pid) = mod(pid,p1)
   proc_col(YBLOCK,pid) = floor(pid/float(p1))
   if(myid==0) write(1,*) pid,proc_row(YBLOCK,pid),proc_col(YBLOCK,pid)
  enddo
  if(myid==0) write(1,*) ' '
  
  if(myid==0) write(1,*) '----------------------------------------------------------------'
  if(myid==0) write(1,*) ' local array sizes/indices '
  if(myid==0) write(1,*) '         n elements      start index       end index      coord'
  if(myid==0) write(1,*) '----------------------------------------------------------------'
  do pid=0,np-1  
   
   array_size(IDIM,YBLOCK,pid) = ny
   global_y_indices(START,YBLOCK,pid)=1
   global_y_indices(END,YBLOCK,pid) =  global_y_indices(START,YBLOCK,pid)   &
                                    +  array_size(IDIM,YBLOCK,pid)-1
   
   if( proc_row(YBLOCK,pid) < p1-1 ) then
    array_size(JDIM,YBLOCK,pid) = floor(nx/float(p1))
   else
    array_size(JDIM,YBLOCK,pid) = floor(nx/float(p1)) + mod(nx,p1)
   endif
   global_x_indices(START,YBLOCK,pid)= proc_row(YBLOCK,pid)*floor(nx/float(p1)) + 1
   global_x_indices(END,YBLOCK,pid) =  global_x_indices(START,YBLOCK,pid)   &
                                    +  array_size(JDIM,YBLOCK,pid)-1
   
   if( proc_col(YBLOCK,pid) < p2-1 ) then
    array_size(KDIM,YBLOCK,pid) = floor(nz/float(p2))
   else
    array_size(KDIM,YBLOCK,pid) = floor(nz/float(p2)) + mod(nz,p2)
   endif
   global_z_indices(START,YBLOCK,pid)= proc_col(YBLOCK,pid)*floor(nz/float(p2)) + 1
   global_z_indices(END,YBLOCK,pid) =  global_z_indices(START,YBLOCK,pid)   &
                                    +  array_size(KDIM,YBLOCK,pid)-1
   
   !------------------------------------------------------------
   !  create index arrays for myid's YBLOCK
   !------------------------------------------------------------

   if(pid==myid) then
    npts=array_size(JDIM,YBLOCK,pid)*array_size(KDIM,YBLOCK,pid)
    allocate( jk_indices_YBLOCK(npts,2) )
    do i=1,npts
     call get_2d_indices(i,array_size(JDIM,YBLOCK,pid),j,k)
     jk_indices_YBLOCK(i,1)=j
     jk_indices_YBLOCK(i,2)=k
    enddo
   endif

   if(pid==myid) then
    npts=array_size(IDIM,YBLOCK,pid)*array_size(JDIM,YBLOCK,pid)*array_size(KDIM,YBLOCK,pid)
    allocate( ijk_indices_YBLOCK(npts,3) )
    do ii=1,npts
     call get_3d_indices(ii,array_size(IDIM,YBLOCK,pid),array_size(JDIM,YBLOCK,pid),i,j,k)
     ijk_indices_YBLOCK(ii,1)=i
     ijk_indices_YBLOCK(ii,2)=j
     ijk_indices_YBLOCK(ii,3)=k
    enddo
   endif


   if(myid==0) write(1,*) pid
   if(myid==0) write(1,*) '  1st index  ',array_size(IDIM,YBLOCK,pid),global_y_indices(:,YBLOCK,pid),'           y'
   if(myid==0) write(1,*) '  2nd index  ',array_size(JDIM,YBLOCK,pid),global_x_indices(:,YBLOCK,pid),'           x'
   if(myid==0) write(1,*) '  3rd index  ',array_size(KDIM,YBLOCK,pid),global_z_indices(:,YBLOCK,pid),'           z'
  enddo
  if(myid==0) write(1,*) ' '

  
  
  !--------------------------------------------------------------------------
  !! ZBLOCK definitions
  !--------------------------------------------------------------------------
  if(myid==0) write(1,*) '-------------------------------------'
  if(myid==0) write(1,*) 'ZBLOCK: '
  if(myid==0) write(1,*) '-------------------------------------'
  if(myid==0) write(1,*) '         pid        row          col '
  if(myid==0) write(1,*) '-------------------------------------'
  do pid=0,np-1
   proc_row(ZBLOCK,pid) = mod(pid,p1)
   proc_col(ZBLOCK,pid) = floor(pid/float(p1))
   if(myid==0) write(1,*) pid,proc_row(ZBLOCK,pid),proc_col(ZBLOCK,pid)
  enddo
  if(myid==0) write(1,*) ' '
  
  if(myid==0) write(1,*) '----------------------------------------------------------------'
  if(myid==0) write(1,*) ' local array sizes/indices '
  if(myid==0) write(1,*) '         n elements      start index       end index      coord'
  if(myid==0) write(1,*) '----------------------------------------------------------------'
  do pid=0,np-1  
   
   array_size(IDIM,ZBLOCK,pid) = nz
   global_z_indices(START,ZBLOCK,pid)=1
   global_z_indices(END,ZBLOCK,pid) =  global_z_indices(START,ZBLOCK,pid)   &
                                    +  array_size(IDIM,ZBLOCK,pid)-1
   
   if( proc_row(ZBLOCK,pid) < p1-1 ) then
    array_size(JDIM,ZBLOCK,pid) = floor(nx/float(p1))
   else
    array_size(JDIM,ZBLOCK,pid) = floor(nx/float(p1)) + mod(nx,p1)
   endif
   global_x_indices(START,ZBLOCK,pid)= proc_row(ZBLOCK,pid)*floor(nx/float(p1)) + 1
   global_x_indices(END,ZBLOCK,pid) =  global_x_indices(START,ZBLOCK,pid)   &
                                    +  array_size(JDIM,ZBLOCK,pid)-1
   
   if( proc_col(ZBLOCK,pid) < p2-1 ) then
    array_size(KDIM,ZBLOCK,pid) = floor(ny/float(p2))
   else
    array_size(KDIM,ZBLOCK,pid) = floor(ny/float(p2)) + mod(ny,p2)
   endif
   global_y_indices(START,ZBLOCK,pid)= proc_col(ZBLOCK,pid)*floor(ny/float(p2)) + 1
   global_y_indices(END,ZBLOCK,pid) =  global_y_indices(START,ZBLOCK,pid)   &
                                    +  array_size(KDIM,ZBLOCK,pid)-1
   

   !------------------------------------------------------------
   !  create index arrays for myid's ZBLOCK
   !------------------------------------------------------------

   if( pid==myid ) then
    npts=array_size(IDIM,ZBLOCK,pid)*array_size(JDIM,ZBLOCK,pid)*array_size(KDIM,ZBLOCK,pid)
    allocate( ijk_indices_ZBLOCK(npts,3) )
    do ii=1,npts
     call get_3d_indices(ii,array_size(IDIM,ZBLOCK,pid),array_size(JDIM,ZBLOCK,pid),i,j,k)
     ijk_indices_ZBLOCK(ii,1)=i
     ijk_indices_ZBLOCK(ii,2)=j
     ijk_indices_ZBLOCK(ii,3)=k
    enddo
   endif


   if(myid==0) write(1,*) pid
   if(myid==0) write(1,*) '  1st index  ',array_size(IDIM,ZBLOCK,pid),global_z_indices(:,ZBLOCK,pid),'           z'
   if(myid==0) write(1,*) '  2nd index  ',array_size(JDIM,ZBLOCK,pid),global_x_indices(:,ZBLOCK,pid),'           x'
   if(myid==0) write(1,*) '  3rd index  ',array_size(KDIM,ZBLOCK,pid),global_y_indices(:,ZBLOCK,pid),'           y'
  enddo
  
  if(myid==0) close(1)

 call mpi_barrier(comm,ierr)
 return
end subroutine decomposition2d



