subroutine xblock_2_yblock(in,out)
!------------------------------------------------------------
!  Transpose operation from XBLOCK format to YBLOCK FORMAT
!  XBLOCK   x vals local, y divided by p1, z divided by p2
!  YBLOCK   y vals local, x divided by p1, z divided by p2
!  ==> if( p1==1 )   
!    both xblock_2_yblock and yblock_2_xblock
!    operations are LOCAL to each of the p2 processors
!------------------------------------------------------------
 use decomposition_params
 use independent_variables, only: nx,ny,nz
 use mpi_params,            only: myid
 implicit none
 integer                       :: i,j,k
 integer,parameter             :: INBLOCK  = XBLOCK
 integer,parameter             :: OUTBLOCK = YBLOCK
 
 real(kind=8)                  ::  in( array_size(IDIM,INBLOCK,myid),   &
                                       array_size(JDIM,INBLOCK,myid),   &
                                       array_size(KDIM,INBLOCK,myid)  )
                                       
 real(kind=8)                  :: out( array_size(IDIM,OUTBLOCK,myid),   &
                                       array_size(JDIM,OUTBLOCK,myid),   &
                                       array_size(KDIM,OUTBLOCK,myid)  )
 
 integer,save                      :: transpose_plan
 integer,dimension(3),save         :: gdims
 integer,dimension(3),save         :: dims_in, dims_out
 integer,dimension(3),save         :: layout_in, layout_out
 integer,dimension(3),save         :: memorder_in, memorder_out
 character(len=4),save             :: datatype='real'
 character(len=4),save             :: flavor='none'
 integer,parameter                 :: dir=-999      
 integer,parameter                 :: isign=-999 
 
 logical,save                      :: first_entry=.TRUE.
 

!--------------------------------------------------------------------- 
!  FIRST ENTRY INITIALIZATION
!--------------------------------------------------------------------- 
 if( first_entry ) then
   gdims=(/nx,ny,nz/)
   dims_in = array_size(:,INBLOCK,myid)
   layout_in  =  layout(:,INBLOCK)
   memorder_in = mem_order(:,INBLOCK)    
   dims_out = array_size(:,OUTBLOCK,myid)
   layout_out  =  layout(:,OUTBLOCK) 
   memorder_out = mem_order(:,OUTBLOCK)
   if( p1 > 1 ) then   

    call par_transpose2D_init(gdims,layout_in,layout_out,  &
                              memorder_in,memorder_out,    &
                              transpose_plan)
   endif   
  first_entry = .FALSE.
 endif
!--------------------------------------------------------------------- 
!  END FIRST ENTRY INITIALIZATION
!--------------------------------------------------------------------- 



!--------------------------------------------------------------------- 
!  for p1=1, XBLOCK<-->YBLOCK transposes are local to a processor, 
!  no need for MPI version of transpose
!---------------------------------------------------------------------
 if( p1==1 ) then 
 
!------------------------------------------
!  transpose => swap 1st & 2nd indices
!------------------------------------------  

  do k=1,array_size(KDIM,INBLOCK,myid)
   out(:,:,k) = TRANSPOSE( in(:,:,k) )
  enddo



!--------------------------------------------------------------------- 
!  for p1 .NE. 1, data is distributed across processors
!  need to call MPI version of transpose
!---------------------------------------------------------------------    
 elseif( p1 > 1 ) then

  if( p2==1 ) stop 'transpose 2D will not do 1d decompositions, xblock_2_yblock'
  call par_transpose2D(in,out,gdims,              &
                       memorder_in,memorder_out,  &
                       layout_in,layout_out,      &
                       transpose_plan)

 elseif( p1 <= 0 ) then
  stop 'p1 <= 0 in xblock_2_yblock'
 endif
 
 return
end subroutine xblock_2_yblock




subroutine yblock_2_xblock(in,out)
!------------------------------------------------------------
!  Transpose operation from YBLOCK format to XBLOCK FORMAT
!  XBLOCK   x vals local, y divided by p1, z divided by p2
!  YBLOCK   y vals local, x divided by p1, z divided by p2
!  ==> if( p1==1 )   
!    both xblock_2_yblock and yblock_2_xblock
!    operations are LOCAL to each of the p2 processors
!------------------------------------------------------------
 use decomposition_params
 use independent_variables, only: nx,ny,nz
 use mpi_params,            only: myid
 implicit none
 integer                       :: i,j,k
 integer,parameter             :: INBLOCK  = YBLOCK
 integer,parameter             :: OUTBLOCK = XBLOCK
 real(kind=8)                  ::  in( array_size(IDIM,INBLOCK,myid),   &
                                       array_size(JDIM,INBLOCK,myid),   &
                                       array_size(KDIM,INBLOCK,myid)  )
                                       
 real(kind=8)                  :: out( array_size(IDIM,OUTBLOCK,myid),   &
                                       array_size(JDIM,OUTBLOCK,myid),   &
                                       array_size(KDIM,OUTBLOCK,myid)  )
 
 integer,save                      :: transpose_plan
 integer,dimension(3),save         :: gdims
 integer,dimension(3),save         :: dims_in, dims_out
 integer,dimension(3),save         :: layout_in, layout_out
 integer,dimension(3),save         :: memorder_in, memorder_out
 character(len=4),save             :: datatype='real'
 character(len=4),save             :: flavor='none'
 integer,parameter                 :: dir=-999      
 integer,parameter                 :: isign=-999 
 
 logical,save                      :: first_entry=.TRUE.
 
!--------------------------------------------------------------------- 
!  FIRST ENTRY INITIALIZATION
!--------------------------------------------------------------------- 
 if( first_entry ) then 
   gdims=(/nx,ny,nz/)
   dims_in = array_size(:,INBLOCK,myid)
   layout_in  =  layout(:,INBLOCK)
   memorder_in = mem_order(:,INBLOCK)      
   dims_out = array_size(:,OUTBLOCK,myid)
   layout_out  =  layout(:,OUTBLOCK)
   memorder_out = mem_order(:,OUTBLOCK)
   
   if( p1 > 1 ) then

    call par_transpose2D_init(gdims,layout_in,layout_out,  &
                              memorder_in,memorder_out,    &
                              transpose_plan)

   endif   
  first_entry = .FALSE.
 endif
!--------------------------------------------------------------------- 
!  END FIRST ENTRY INITIALIZATION
!---------------------------------------------------------------------



!--------------------------------------------------------------------- 
!  for p1=1, XBLOCK<-->YBLOCK transposes are local to a processor, 
!  no need for MPI version of transpose
!---------------------------------------------------------------------                                        
 if( p1==1 ) then 
 
!------------------------------------------
!  transpose => swap 1st & 2nd indices
!------------------------------------------

  do k=1,array_size(KDIM,INBLOCK,myid)
   out(:,:,k) = TRANSPOSE( in(:,:,k) )
  enddo


 elseif( p1 > 1 ) then

  if( p2==1 ) stop 'transpose 2D will not do 1d decompositions, yblock_2_xblock'
  call par_transpose2D(in,out,gdims,              &
                       memorder_in,memorder_out,  &
                       layout_in,layout_out,      &
                       transpose_plan)

 elseif( p1 <= 0 ) then
  stop 'p1 <= 0 in yblock_2_xblock'
 endif
 
 return
end subroutine yblock_2_xblock




subroutine yblock_2_zblock(in,out)
!------------------------------------------------------------
!  Transpose operation from YBLOCK format to ZBLOCK FORMAT
!  y and z are both divided by p2
!  YBLOCK   y vals local, x divided by p1, z divided by p2
!  ZBLOCK   z vals local, x divided by p1, y divided by p2
!------------------------------------------------------------
 use decomposition_params
 use independent_variables, only: nx,ny,nz
 use mpi_params,            only: myid,comm,ierr
 implicit none
 integer                       :: i,j,k,npts,ii
 integer,parameter             :: INBLOCK  = YBLOCK
 integer,parameter             :: OUTBLOCK = ZBLOCK
 real(kind=8)                  ::  in( array_size(IDIM,INBLOCK,myid),   &
                                       array_size(JDIM,INBLOCK,myid),   &
                                       array_size(KDIM,INBLOCK,myid)  )
                                       
 real(kind=8)                  :: out( array_size(IDIM,OUTBLOCK,myid),   &
                                       array_size(JDIM,OUTBLOCK,myid),   &
                                       array_size(KDIM,OUTBLOCK,myid)  )
 
 integer,save                      :: transpose_plan
 integer,dimension(3),save         :: gdims
 integer,dimension(3),save         :: dims_in, dims_out
 integer,dimension(3),save         :: layout_in, layout_out
 integer,dimension(3),save         :: memorder_in, memorder_out
 character(len=4),save             :: datatype='real'
 character(len=4),save             :: flavor='none'
 integer,parameter                 :: dir=-999      
 integer,parameter                 :: isign=-999 
 
 logical,save                      :: first_entry=.TRUE.
 
!--------------------------------------------------------------------- 
!  FIRST ENTRY INITIALIZATION
!---------------------------------------------------------------------  
 if( first_entry ) then 
   gdims=(/nx,ny,nz/)
   dims_in = array_size(:,INBLOCK,myid)
   layout_in  =  layout(:,INBLOCK)
   memorder_in = mem_order(:,INBLOCK)      
   dims_out = array_size(:,OUTBLOCK,myid)
   layout_out  =  layout(:,OUTBLOCK)
   memorder_out = mem_order(:,OUTBLOCK)
   if( np > 1 ) then 

   if( p1 > 1 ) then
    call par_transpose2D_init(gdims,layout_in,layout_out,  &
                              memorder_in,memorder_out,    &
                              transpose_plan) 
   endif

   endif
  first_entry = .FALSE.
 endif 
!--------------------------------------------------------------------- 
!  END FIRST ENTRY INITIALIZATION
!--------------------------------------------------------------------- 



!--------------------------------------------------------------------- 
!  for p2=1, YBLOCK<-->ZBLOCK transposes are local to a processor, 
!  no need for MPI version of transpose
!---------------------------------------------------------------------
 if( p2==1 ) then 
!------------------------------------------
!  transpose => swap 1st & 3rd indices
!------------------------------------------

 do k=1,array_size(KDIM,INBLOCK,myid)
  do j=1,array_size(JDIM,INBLOCK,myid)
   do i=1,array_size(IDIM,INBLOCK,myid)
    out(k,j,i) = in(i,j,k)
   enddo
  enddo
 enddo

  
 elseif( p2 > 1 ) then

  if( p1==1 ) then
   call yblock_2_zblock_1d(in,out)
  else
   call par_transpose2D(in,out,gdims,              &
                        memorder_in,memorder_out,  &
                        layout_in,layout_out,      &
                        transpose_plan)
  endif

 elseif( p2 <= 0 ) then
  stop 'p2 <= 0 in yblock_2_zblock'
 endif
 
 return
end subroutine yblock_2_zblock



subroutine zblock_2_yblock(in,out)
!------------------------------------------------------------
!  Transpose operation from YBLOCK format to ZBLOCK FORMAT
!  y and z are both divided by p2
!  YBLOCK   y vals local, x divided by p1, z divided by p2
!  ZBLOCK   z vals local, x divided by p1, y divided by p2
!------------------------------------------------------------
 use decomposition_params
 use independent_variables, only: nx,ny,nz
 use mpi_params,            only: myid
 implicit none
 integer                       :: i,j,k,npts,ii
 integer,parameter             :: INBLOCK  = ZBLOCK
 integer,parameter             :: OUTBLOCK = YBLOCK
 real(kind=8)                  ::  in( array_size(IDIM,INBLOCK,myid),   &
                                       array_size(JDIM,INBLOCK,myid),   &
                                       array_size(KDIM,INBLOCK,myid)  )
                                       
 real(kind=8)                  :: out( array_size(IDIM,OUTBLOCK,myid),   &
                                       array_size(JDIM,OUTBLOCK,myid),   &
                                       array_size(KDIM,OUTBLOCK,myid)  )
 
 integer,save                      :: transpose_plan
 integer,dimension(3),save         :: gdims
 integer,dimension(3),save         :: dims_in, dims_out
 integer,dimension(3),save         :: layout_in, layout_out
 integer,dimension(3),save         :: memorder_in, memorder_out
 character(len=4),save             :: datatype='real'
 character(len=4),save             :: flavor='none'
 integer,parameter                 :: dir=-999      
 integer,parameter                 :: isign=-999 
 
 logical,save                      :: first_entry=.TRUE.
 
!--------------------------------------------------------------------- 
!  FIRST ENTRY INITIALIZATION
!--------------------------------------------------------------------- 
 if( first_entry ) then 
   gdims=(/nx,ny,nz/)
   dims_in = array_size(:,INBLOCK,myid)
   layout_in  =  layout(:,INBLOCK)
   memorder_in = mem_order(:,INBLOCK)      
   dims_out = array_size(:,OUTBLOCK,myid)
   layout_out  =  layout(:,OUTBLOCK)
   memorder_out = mem_order(:,OUTBLOCK)
   if( np > 1 ) then

   if( p1 > 1 ) then
    call par_transpose2D_init(gdims,layout_in,layout_out,  &
                              memorder_in,memorder_out,    &
                              transpose_plan)
   endif

   endif
  first_entry = .FALSE.
 endif  
!--------------------------------------------------------------------- 
!  END FIRST ENTRY INITIALIZATION
!---------------------------------------------------------------------


!--------------------------------------------------------------------- 
!  for p2=1, YBLOCK<-->ZBLOCK transposes are local to a processor, 
!  no need for MPI version of transpose
!---------------------------------------------------------------------
 if( p2==1 ) then
 
!------------------------------------------
!  transpose => swap 1st & 3rd indices
!------------------------------------------

 do k=1,array_size(KDIM,INBLOCK,myid)
  do j=1,array_size(JDIM,INBLOCK,myid)
   do i=1,array_size(IDIM,INBLOCK,myid)
    out(k,j,i) = in(i,j,k)
   enddo
  enddo
 enddo

    
 elseif( p2 > 1 ) then


  if( p1==1 ) then
   call zblock_2_yblock_1d(in,out)
  else
   call par_transpose2D(in,out,gdims,              &
                        memorder_in,memorder_out,  &
                        layout_in,layout_out,      &
                        transpose_plan)
  endif

 elseif( p2 <= 0 ) then
  stop 'p2 <= 0 in zblock_2_yblock'
 endif
 
 return
end subroutine zblock_2_yblock





subroutine yblock_2_zblock_1d(in,out)
!----------------------------------------------------------------
!   YBLOCK to ZBLOCK, assuming 1d decomposiiton p1=1
!----------------------------------------------------------------
  use mpi_params,               only: myid,comm,ierr
  use decomposition_params
  
  implicit none
  include 'mpif.h'
  integer,parameter                :: INBLOCK  = YBLOCK
  integer,parameter                :: OUTBLOCK = ZBLOCK
  integer                          :: pid,i,j,k,position 
  integer                          :: sum_send,sum_recv
  
  real(kind=8)                     ::  in( array_size(IDIM,INBLOCK,myid),   &
                                           array_size(JDIM,INBLOCK,myid),   &
                                           array_size(KDIM,INBLOCK,myid)  )
                                       
  real(kind=8)                     :: out( array_size(IDIM,OUTBLOCK,myid),   &
                                           array_size(JDIM,OUTBLOCK,myid),   &
                                           array_size(KDIM,OUTBLOCK,myid)  )
  
  real(kind=8),allocatable,save         :: sendbuf(:),recvbuf(:)
  integer,allocatable,dimension(:),save :: ystart,yend,zstart,zend
  integer,allocatable,dimension(:),save :: size_send_chunk,size_recv_chunk
  integer,allocatable,dimension(:),save :: send_displs,recv_displs
 
  logical,save                          :: first_entry=.true.
  logical,save                          :: verbose=.false.

!---------------------------------------------------------------------------- 
!      FIRST ENTRY BLOCK
!----------------------------------------------------------------------------
  if( first_entry ) then 
  
   allocate( sendbuf(PRODUCT(array_size(:,INBLOCK,myid))) )
   allocate( recvbuf(PRODUCT(array_size(:,OUTBLOCK,myid))) )
   allocate( ystart(0:np-1), yend(0:np-1) )
   allocate( zstart(0:np-1), zend(0:np-1) )
   allocate( size_send_chunk(0:np-1), size_recv_chunk(0:np-1) )
   allocate( send_displs(0:np-1), recv_displs(0:np-1) )
   
   sum_send=0
   sum_recv=0
   
   do pid=0,np-1
   
    !------------------------------------------------
    ! size of chunk to be sent to processor pid
    !------------------------------------------------
    size_send_chunk(pid) = array_size(KDIM,INBLOCK,myid)     &   ! all local z indices
                          *array_size(JDIM,INBLOCK,myid)     &   ! all local x indices
                          *array_size(KDIM,OUTBLOCK,pid)         ! some of the local y indices
    sum_send = sum_send + size_send_chunk(pid)
    
    !------------------------------------------------
    ! size of chunk to be recvd from processor pid
    !------------------------------------------------
    size_recv_chunk(pid) = array_size(KDIM,OUTBLOCK,myid)    &   ! all local y indices
                          *array_size(JDIM,OUTBLOCK,myid)    &   ! all local x indices
                          *array_size(KDIM,INBLOCK,pid)          ! some of the local z indices                          
    sum_recv = sum_recv + size_recv_chunk(pid)
    
    !------------------------------------------------
    ! y indices of chunk to be sent to processor pid
    !------------------------------------------------
    ystart(pid) = global_y_indices(START,OUTBLOCK,pid)
    yend(pid)   = global_y_indices(END,OUTBLOCK,pid)
    
    !------------------------------------------------
    ! z indices of chunk to be recvd from processor pid
    !------------------------------------------------
    zstart(pid) = global_z_indices(START,INBLOCK,pid)
    zend(pid)   = global_z_indices(END,INBLOCK,pid)
    
    !------------------------------------------------
    ! displacements for chunk sent to processor pid
    !  (uses convention that only last chunk
    !  is allowed to be odd sized)
    !------------------------------------------------
    send_displs(pid) = size_send_chunk(0)*pid
                      
    
    !------------------------------------------------
    ! displacements for chunk recvd from processor pid
    !  (uses convention that only last chunk
    !  is allowed to be odd sized)
    !------------------------------------------------
    recv_displs(pid) = size_recv_chunk(0)*pid
                      
   enddo
   
   if( sum_send .ne. PRODUCT(array_size(:,INBLOCK,myid)) ) then
    write(0,*) 'send chunks dont add up ',myid,sum_send,PRODUCT(array_size(:,INBLOCK,myid))
    write(0,*) 'array size YBLOCK ',myid,array_size(:,INBLOCK,myid)
    write(0,*) 'size_send_chunks: ',myid,size_send_chunk(:)
    stop
   endif
    
   if( sum_recv .ne. PRODUCT(array_size(:,OUTBLOCK,myid)) ) then
    write(0,*) 'recv chunks dont add up ',myid,sum_recv,PRODUCT(array_size(:,OUTBLOCK,myid))
    write(0,*) 'array size ZBLOCK ',myid,array_size(:,OUTBLOCK,myid)
    write(0,*) 'size_recv_chunks: ',myid,size_recv_chunk(:)
    stop
   endif
  
   if( verbose ) then
    write(0,*) myid,' size of my INBLOCK ',array_size(:,INBLOCK,myid)
    write(0,*) myid,' size of my sendbuf ',SIZE(sendbuf) 
    write(0,*) myid,' size of sendchunks ',size_send_chunk(:)
    write(0,*) myid,' send displs ',send_displs(:)
   
    write(0,*) myid,' size of my OUTBLOCK ',array_size(:,OUTBLOCK,myid) 
    write(0,*) myid,' size of my recvbuf ',SIZE(recvbuf)
    write(0,*) myid,' size of recvchunks ',size_recv_chunk(:)
    write(0,*) myid,' recv displs ',recv_displs(:)
   endif
   
   
   first_entry=.false.
  endif
!---------------------------------------------------------------------------- 
!      END FIRST ENTRY BLOCK
!----------------------------------------------------------------------------

    if( p1 .ne. 1 ) stop 'p1 must be 1 to use yblock_2_zblock_1d'
    if( p2 == 1 ) stop 'no need for mpi 1d transpose when p1=p2=1: yblock_2_zblock_1d' 

    !------------------------------------------------
    ! pack up the send buffer
    ! convention: z, then x, then y vary fastest w/ position
    !------------------------------------------------
    position = 1
    do pid = 0,np-1                          ! dest pid           
     do j=ystart(pid),yend(pid)              ! y coord
      do i=1,array_size(JDIM,INBLOCK,myid)   ! x coord
       do k=1,array_size(KDIM,INBLOCK,myid)  ! z coord
      
        sendbuf(position) = IN(j,i,k)        ! IN is YBLOCK
        position = position + 1
       
       enddo
      enddo
     enddo
    enddo
    
    
    !------------------------------------------------
    ! do the mpi communications
    !------------------------------------------------
    call MPI_ALLTOALLV(sendbuf,size_send_chunk,send_displs,MPI_DOUBLE_PRECISION,   &
                       recvbuf,size_recv_chunk,recv_displs,MPI_DOUBLE_PRECISION,   &
                       comm,ierr)
    if( ierr .ne. 0 ) stop 'mpi_alltoallv error yblock_2_zblock_1d'
    
    
    !------------------------------------------------
    ! unpack up the recv buffer
    ! ===> must keep convention that
    ! convention: z, then x, then y vary fastest w/ position
    !------------------------------------------------
    position = 1
    do pid = 0,np-1                           ! source pid 
     do j=1,array_size(KDIM,OUTBLOCK,myid)    ! y coord   
      do i=1,array_size(JDIM,OUTBLOCK,myid)   ! x coord 
       do k=zstart(pid),zend(pid)             ! z coord             
      
        OUT(k,i,j) = recvbuf(position)        ! OUT is ZBLOCK
        position = position + 1
       
       enddo
      enddo
     enddo
    enddo


return
end subroutine yblock_2_zblock_1d


















subroutine zblock_2_yblock_1d(in,out)
!----------------------------------------------------------------
!   ZBLOCK to YBLOCK, assuming 1d decomposiiton p1=1
!----------------------------------------------------------------
  use mpi_params,               only: myid,comm,ierr
  use decomposition_params
  
  implicit none
  include 'mpif.h'
  integer,parameter                :: INBLOCK  = ZBLOCK
  integer,parameter                :: OUTBLOCK = YBLOCK
  integer                          :: pid,i,j,k,position
  integer                          :: sum_send,sum_recv
  
  real(kind=8)                     ::  in( array_size(IDIM,INBLOCK,myid),   &
                                           array_size(JDIM,INBLOCK,myid),   &
                                           array_size(KDIM,INBLOCK,myid)  )
                                       
  real(kind=8)                     :: out( array_size(IDIM,OUTBLOCK,myid),   &
                                           array_size(JDIM,OUTBLOCK,myid),   &
                                           array_size(KDIM,OUTBLOCK,myid)  )
  
  real(kind=8),allocatable,save         :: sendbuf(:),recvbuf(:)
  integer,allocatable,dimension(:),save :: ystart,yend,zstart,zend
  integer,allocatable,dimension(:),save :: size_send_chunk,size_recv_chunk
  integer,allocatable,dimension(:),save :: send_displs,recv_displs
 
  logical,save                          :: first_entry=.true.
  logical,save                          :: verbose=.false.

!---------------------------------------------------------------------------- 
!      FIRST ENTRY BLOCK
!----------------------------------------------------------------------------
  if( first_entry ) then 
  
   allocate( sendbuf(PRODUCT(array_size(:,INBLOCK,myid))) )
   allocate( recvbuf(PRODUCT(array_size(:,OUTBLOCK,myid))) )
   allocate( ystart(0:np-1), yend(0:np-1) )
   allocate( zstart(0:np-1), zend(0:np-1) )
   allocate( size_send_chunk(0:np-1), size_recv_chunk(0:np-1) )
   allocate( send_displs(0:np-1), recv_displs(0:np-1) )
   
   sum_send=0
   sum_recv=0
       
   do pid=0,np-1
   
    !------------------------------------------------
    ! size of chunk to be sent to processor pid
    !------------------------------------------------
    size_send_chunk(pid) = array_size(KDIM,INBLOCK,myid)     &
                          *array_size(JDIM,INBLOCK,myid)     &
                          *array_size(KDIM,OUTBLOCK,pid)
    sum_send = sum_send + size_send_chunk(pid)
    
    !------------------------------------------------
    ! size of chunk to be recvd from processor pid
    !------------------------------------------------
    size_recv_chunk(pid) = array_size(KDIM,OUTBLOCK,myid)     &
                          *array_size(JDIM,OUTBLOCK,myid)     &
                          *array_size(KDIM,INBLOCK,pid)
    sum_recv = sum_recv + size_recv_chunk(pid)
    
    !------------------------------------------------
    ! z indices of chunk to be sent to processor pid
    !------------------------------------------------
    zstart(pid) = global_z_indices(START,OUTBLOCK,pid)
    zend(pid)   = global_z_indices(END,OUTBLOCK,pid)
    
    !------------------------------------------------
    ! y indices of chunk to be recvd from processor pid
    !------------------------------------------------
    ystart(pid) = global_y_indices(START,INBLOCK,pid)
    yend(pid)   = global_y_indices(END,INBLOCK,pid)
    
    !------------------------------------------------
    ! displacements for chunk sent to processor pid
    !  (uses convention that only last chunk
    !  is allowed to be odd sized)
    !------------------------------------------------
    send_displs(pid) = array_size(KDIM,INBLOCK,myid)      &
                      *array_size(JDIM,INBLOCK,myid)      &
                      *array_size(KDIM,OUTBLOCK,0)*pid
                      
    
    !------------------------------------------------
    ! displacements for chunk recvd from processor pid
    !  (uses convention that only last chunk
    !  is allowed to be odd sized)
    !------------------------------------------------
    recv_displs(pid) = array_size(KDIM,OUTBLOCK,myid)      &
                      *array_size(JDIM,OUTBLOCK,myid)      &
                      *array_size(KDIM,INBLOCK,0)*pid
                      
   enddo
   
      if( sum_send .ne. PRODUCT(array_size(:,INBLOCK,myid)) ) then
    write(0,*) 'send chunks dont add up ',myid,sum_send,PRODUCT(array_size(:,INBLOCK,myid))
    write(0,*) 'array size YBLOCK ',myid,array_size(:,INBLOCK,myid)
    write(0,*) 'size_send_chunks: ',myid,size_send_chunk(:)
    stop
   endif
    
   if( sum_recv .ne. PRODUCT(array_size(:,OUTBLOCK,myid)) ) then
    write(0,*) 'recv chunks dont add up ',myid,sum_recv,PRODUCT(array_size(:,OUTBLOCK,myid))
    write(0,*) 'array size ZBLOCK ',myid,array_size(:,OUTBLOCK,myid)
    write(0,*) 'size_recv_chunks: ',myid,size_recv_chunk(:)
    stop
   endif
   
   if( verbose ) then
    write(0,*) myid,' size of my INBLOCK ',array_size(:,INBLOCK,myid)
    write(0,*) myid,' size of my sendbuf ',SIZE(sendbuf) 
    write(0,*) myid,' size of sendchunks ',size_send_chunk(:)
    write(0,*) myid,' send displs ',send_displs(:)
   
    write(0,*) myid,' size of my OUTBLOCK ',array_size(:,OUTBLOCK,myid) 
    write(0,*) myid,' size of my recvbuf ',SIZE(recvbuf)
    write(0,*) myid,' size of recvchunks ',size_recv_chunk(:)
    write(0,*) myid,' recv displs ',recv_displs(:)
   endif
   
   first_entry=.false.
  endif
!---------------------------------------------------------------------------- 
!      END FIRST ENTRY BLOCK
!----------------------------------------------------------------------------

    if( p1 .ne. 1 ) stop 'p1 must be 1 to use yblock_2_zblock_1d'
    if( p2 == 1 ) stop 'no need for mpi 1d transpose when p1=p2=1: yblock_2_zblock_1d' 

    !------------------------------------------------
    ! pack up the send buffer
    ! convention: y,x,z vary fastest w/ position
    !------------------------------------------------
    position = 1
    do pid = 0,np-1                          ! dest pid
     do k=zstart(pid),zend(pid)              ! z coord
      do i=1,array_size(JDIM,INBLOCK,myid)   ! x coord
       do j=1,array_size(KDIM,INBLOCK,myid)  ! y coord
      
        sendbuf(position) = IN(k,i,j)        ! IN is ZBLOCK
        position = position + 1
       
       enddo
      enddo
     enddo
    enddo
    
    
    
    !------------------------------------------------
    ! do the mpi communications
    !------------------------------------------------
    call MPI_ALLTOALLV(sendbuf,size_send_chunk,send_displs,MPI_DOUBLE_PRECISION,   &
                       recvbuf,size_recv_chunk,recv_displs,MPI_DOUBLE_PRECISION,   &
                       comm,ierr)
    if( ierr .ne. 0 ) stop 'mpi_alltoallv error zblock_2_yblock_1d'
    
    !------------------------------------------------
    ! unpack up the recv buffer
    ! ===> must keep convention that
    ! y,x,z coordinates vary fastest w/ position
    !------------------------------------------------
    position = 1
    do pid = 0,np-1                           ! source pid      
     do k=1,array_size(KDIM,OUTBLOCK,myid)    ! z coord
      do i=1,array_size(JDIM,OUTBLOCK,myid)   ! x coord
       do j=ystart(pid),yend(pid)             ! y coord
      
        OUT(j,i,k) = recvbuf(position)       ! OUT is YBLOCK
        position = position + 1
       
       enddo
      enddo
     enddo
    enddo



return
end subroutine zblock_2_yblock_1d
