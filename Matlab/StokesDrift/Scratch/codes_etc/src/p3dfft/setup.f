C ------------------------------------------------------------
      subroutine mpisetup_2D(p1,p2,dims,s)

      use comm

      integer p1,p2,pm,dims(3),s,i,iisize,ijsize,jisize,jjsize,kisize,kjsize,gran
      integer p_id(2),cartid(2),ierr,mpi_comm_cart
      integer taskid,pg(2),size(3,2)
      logical periodic(2),remain_dims(2)


      do s=1,ns
         if(gar_dim(1,s) .ne. dims(1) .or. gar_dim(2,s) .ne. dims(2) .or. gar_dim(3,s) .ne. dims(3)) then
            goto 99
         endif
         if(proc_grid(1,s) .ne. p1 .or. proc_grid(2,s) .ne. p2) then
            goto 99
         endif
         return
 99      continue
      enddo

      s = ns+1
      if(s .gt. MAX_NS) then
         if(ls .gt. MAX_NS) then
            ls = 1
         endif
         s = ls
         ls = ls + 1
         deallocate (st(s)%R)
         deallocate (sz(s)%R)
         deallocate (en(s)%R)
         do i=1,3
            do j=i+1,3
               deallocate(Exch(i,j,s)%RSndCnts)
               deallocate(Exch(i,j,s)%RRcvCnts)
               deallocate(Exch(i,j,s)%RSndStrt)
               deallocate(Exch(i,j,s)%RRcvStrt)
               deallocate(Exch(j,i,s)%RSndCnts)
               deallocate(Exch(j,i,s)%RRcvCnts)
               deallocate(Exch(j,i,s)%RSndStrt)
               deallocate(Exch(j,i,s)%RRcvStrt)
            enddo
         enddo
      else
         ns = ns+1
      endif

      mpi_set(s) = .true.
      gar_dim(1,s) = dims(1)
      gar_dim(2,s) = dims(2)
      gar_dim(3,s) = dims(3)
      proc_grid(1,s) = p1
      proc_grid(2,s) = p2
      proc_grid(3,s) = 1
c      nxh = nx(s)/2
c      nxhp = nxh+1
      pg(1) = p2
      pg(2) = p1

      call MPI_COMM_SIZE (MPI_COMM_WORLD,numtasks,ierr)
      if(p1 * p2 .ne. numtasks) then
         print *,'Seeing ',numtasks,' tasks instead of ',p1*p2
      endif
      call MPI_COMM_RANK (MPI_COMM_WORLD,taskid,ierr)
      
      periodic(1) = .true.
      periodic(2) = .true.

c      print *,'Processor grid: ',proc_grid(1:2,s)
! creating cartesian processor grid
      call MPI_Cart_create(MPI_COMM_WORLD,2,pg,periodic,
     &     .false.,mpi_comm_cart,ierr)
! Obtaining process ids with in the cartesian grid
      call MPI_Cart_coords(mpi_comm_cart,taskid,2,cartid,ierr)

      p_id(1) = cartid(2)
      p_id(2) = cartid(1)
      gpid(1,s) = p_id(1)
      gpid(2,s) = p_id(2)
c      ip_id(s) = ipid
c      jp_id(s) = jpid

c      print *,'ipid,jpid: ',p_id(1),p_id(2)
      
      remain_dims(1) = .true.
      remain_dims(2) = .false.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col(s),ierr)
! using cart comworld create northc-south(column) sub comworld
      remain_dims(1) = .false.
      remain_dims(2) = .true.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row(s),ierr)
      
      if(p1 .gt. p2) then
         pm = p1
      else
         pm = p2
      endif

      allocate (st(s)%R(0:pm-1,3,2))
      allocate (sz(s)%R(0:pm-1,3,2))
      allocate (en(s)%R(0:pm-1,3,2))

      do i=1,3
         do j=1,2
            call MapDataToProc(dims(i),proc_grid(j,s),st(s)%R(0,i,j),
     &            en(s)%R(0,i,j), sz(s)%R(0,i,j))
         enddo
c         print *,'s=',s,', start,size 1:', st(s)%R(:,i,1),sz(s)%R(:,i,1)
c         print *,'start,size 2:', st(s)%R(:,i,2),sz(s)%R(:,i,2)
      enddo

      
c      iistart = iist(ipid,s)
c      jjstart = jjst(jpid,s)
c      ijstart = ijst(jpid,s)
c      kistart = kist(ipid,s)
c      kjstart = kjst(jpid,s)

      do j=1,2
         do i=1,3
            size(i,j) = sz(s)%R(p_id(j),i,j)
         enddo
      enddo

c      iisize = sz(s)%R(ipid,1,1)
c      ijsize = sz(s)%R(jpid,1,2)
c      jisize = sz(s)%R(ipid,2,1)
c      jjsize = sz(s)%R(jpid,2,2)
c      kisize = sz(s)%R(ipid,3,1)
c      kjsize = sz(s)%R(jpid,3,2)

c      print *,'iisize,jisize: ',iisize,jisize
c      iiend = iien(ipid,s)
c      jjend = jjen(jpid,s)
c      ijend = ijen(jpid,s)
c      kiend = kien(ipid,s)
c      kjend = kjen(jpid,s)

      do i=1,3
         do j=i+1,3
            allocate(Exch(i,j,s)%RSndCnts(0:pm-1,2))
            allocate(Exch(j,i,s)%RSndCnts(0:pm-1,2))
            allocate(Exch(i,j,s)%RRcvCnts(0:pm-1,2))
            allocate(Exch(j,i,s)%RRcvCnts(0:pm-1,2))
            allocate(Exch(i,j,s)%RSndStrt(0:pm-1,2))
            allocate(Exch(j,i,s)%RSndStrt(0:pm-1,2))
            allocate(Exch(i,j,s)%RRcvStrt(0:pm-1,2))
            allocate(Exch(j,i,s)%RRcvStrt(0:pm-1,2))
         enddo
      enddo

c      allocate (SndCnts_12(s)%R(0:pm-1,2))
c      allocate (SndStrt_12(s)%R(0:pm-1,2))
c      allocate (SndCnts_21(s)%R(0:pm-1,2))
c      allocate (SndStrt_21(s)%R(0:pm-1,2))
c      allocate (SndCnts_13(s)%R(0:pm-1,2))
c      allocate (SndStrt_13(s)%R(0:pm-1,2))
c      allocate (SndCnts_31(s)%R(0:pm-1,2))
c      allocate (SndStrt_31(s)%R(0:pm-1,2))
c      allocate (SndCnts_23(s)%R(0:pm-1,2))
c      allocate (SndStrt_23(s)%R(0:pm-1,2))
c      allocate (SndCnts_32(s)%R(0:pm-1,2))
c      allocate (SndStrt_32(s)%R(0:pm-1,2))

c      allocate (RcvCnts_12(s)%R(0:pm-1,2))
c      allocate (RcvStrt_12(s)%R(0:pm-1,2))
c      allocate (RcvCnts_21(s)%R(0:pm-1,2))
c      allocate (RcvStrt_21(s)%R(0:pm-1,2))
c      allocate (RcvCnts_13(s)%R(0:pm-1,2))
c      allocate (RcvStrt_13(s)%R(0:pm-1,2))
c      allocate (RcvCnts_31(s)%R(0:pm-1,2))
c      allocate (RcvStrt_31(s)%R(0:pm-1,2))
c      allocate (RcvCnts_23(s)%R(0:pm-1,2))
c      allocate (RcvStrt_23(s)%R(0:pm-1,2))
c      allocate (RcvCnts_32(s)%R(0:pm-1,2))
c      allocate (RcvStrt_32(s)%R(0:pm-1,2))

c      gran = kjsize / gblock
c      if(gran .lt. 1) then
c         gran = 1
c      endif


      do i=1,3
         do j=i+1,3
            
c            gran = size(xor(i,j),2) / gblock
c            if(gran .lt. 1) then
c               gran = 1
c            endif

            do ii=0,p1-1

         Exch(i,j,s)%RSndCnts(ii,1) = sz(s)%R(ii,i,1) * size(j,1) *mytype 
         Exch(i,j,s)%RSndStrt(ii,1) = (st(s)%R(ii,i,1)-1) * size(j,1) *mytype 

         Exch(i,j,s)%RRcvCnts(ii,1) = sz(s)%R(ii,j,1) * size(i,1) *mytype 
         Exch(i,j,s)%RRcvStrt(ii,1) = (st(s)%R(ii,j,1)-1) * size(i,1) *mytype
         Exch(j,i,s)%RSndCnts(ii,1) = sz(s)%R(ii,j,1) * size(i,1) *mytype 
         Exch(j,i,s)%RSndStrt(ii,1) = (st(s)%R(ii,j,1)-1) * size(i,1) *mytype

         Exch(j,i,s)%RRcvCnts(ii,1) = sz(s)%R(ii,i,1) * size(j,1) *mytype 
         Exch(j,i,s)%RRcvStrt(ii,1) = (st(s)%R(ii,i,1)-1) * size(j,1) *mytype 
      enddo
      enddo
      enddo


      do i=1,3
         do j=i+1,3
         
c            gran = size(xor(i,j),1) / gblock
c            if(gran .lt. 1) then
c               gran = 1
c            endif

            do ii=0,p2-1

         Exch(i,j,s)%RSndCnts(ii,2) = sz(s)%R(ii,i,2) * size(j,2) *mytype 
         Exch(i,j,s)%RSndStrt(ii,2) = (st(s)%R(ii,i,2)-1) * size(j,2) *mytype

         Exch(i,j,s)%RRcvCnts(ii,2) = sz(s)%R(ii,j,2) * size(i,2) *mytype 
         Exch(i,j,s)%RRcvStrt(ii,2) = (st(s)%R(ii,j,2)-1) * size(i,2) *mytype
         Exch(j,i,s)%RSndCnts(ii,2) = sz(s)%R(ii,j,2) * size(i,2) *mytype 
         Exch(j,i,s)%RSndStrt(ii,2) = (st(s)%R(ii,j,2)-1) * size(i,2) *mytype

         Exch(j,i,s)%RRcvCnts(ii,2) = sz(s)%R(ii,i,2) * size(j,2) *mytype 
         Exch(j,i,s)%RRcvStrt(ii,2) = (st(s)%R(ii,i,2)-1) * size(j,2) *mytype
      enddo
      enddo
      enddo


c         SndCnts_12(s)%R(i,1) = sz(s)%R(i,1,1) * jisize * mytype 
c         SndStrt_12(s)%R(i,1) = (st(s)%R(i,1,1) -1)*jisize*mytype 

c         RcvCnts_12(s)%R(i,1) = sz(s)%R(i,2,1) * iisize * mytype 
c         RcvStrt_12(s)%R(i,1) = (st(s)%R(i,2,1) -1)*iisize*mytype
 
c         SndCnts_21(s)%R(i,1) = sz(s)%R(i,2,1) * iisize * mytype 
c         SndStrt_21(s)%R(i,1) = (st(s)%R(i,2,1) -1)*iisize*mytype

c         RcvCnts_21(s)%R(i,1) = sz(s)%R(i,1,1) * jisize * mytype 
c         RcvStrt_21(s)%R(i,1) = (st(s)%R(i,1,1) -1)*jisize*mytype


      return
      end subroutine
C ------------------------------------------------------------

!==================================================================       
      subroutine MapDataToProc (data,proc,st,en,sz)
c========================================================
!    
       implicit none
       integer data,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
       integer i,size,nl,nu

       size=data/proc
       nu = data - size * proc
       nl = proc - nu
       st(0) = 1
       sz(0) = size
       en(0) = size
       do i=1,nl-1
         st(i) = st(i-1) + size
         sz(i) = size
         en(i) = en(i-1) + size
      enddo
      size = size + 1
      do i=nl,proc-1
         st(i) = en(i-1) + 1
         sz(i) = size
         en(i) = en(i-1) + size
      enddo
      en(proc-1)= data 
      sz(proc-1)= data-st(proc-1)+1

      end subroutine

C 
C ------------------------------------------------------------
      subroutine get_dims(layout,mem_order, ext,start,dp,s)
C ------------------------------------------------------------
c Evaluate and return local dimensions for input and output for a given plan 
c (scheme) s
c Also evaluate processor geometry dp

      use comm

      integer layout(3),mem_order(3),ext(3),dp(3),s
      integer i,start(3)
      integer proc_ndim,ip0,ip1,ip2,nd,gp0,gp1,gp2

      do i=1,3
         ext(mem_order(i)) = gar_dim(i,s) / layout(i)
      enddo
      
      dp(:) = proc_grid(:,s)

      nd = proc_ndim(layout)
c      print *,'layout=',layout,', nd=',nd
      if(nd .eq. 1) then
         if(layout(1) .ne. 1) then
            ip1 = mem_order(1)
         else if(layout(2) .ne. 1) then
            ip1 = mem_order(2)
         else if(layout(3) .ne. 1) then
            ip1 = mem_order(3)
         endif
         start(ip1) = st(s)%R(gpid(1,s),ip1,1)
         ext(ip1) = sz(s)%R(gpid(1,s),ip1,1)
         i = mod(ip1+1,3)
         start(i) = 1
         i = mod(i+1,3)
         start(i) = 1

      else if(nd .eq. 2) then
         if(layout(1) .eq. 1) then
            gp0 = 1
            gp1 = 2
            gp2 = 3
            ip0 = mem_order(1)
            ip1 = mem_order(2)
            ip2 = mem_order(3)
         else if(layout(2) .eq. 1) then
            gp0 = 2
            gp1 = 1
            gp2 = 3
            ip0 = mem_order(2)
            ip1 = mem_order(1)
            ip2 = mem_order(3)
         else
            gp0 = 3
            gp1 = 1
            gp2 = 2
            ip0 = mem_order(3)
            ip1 = mem_order(1)
            ip2 = mem_order(2)
         endif

c         print *,'ip0,ip1,ip2:',ip0,ip1,ip2
c         print *,'ipid,jpid: ',gpid(1,s),gpid(2,s)
         start(ip1) = st(s)%R(gpid(1,s),gp1,1)
         start(ip2) = st(s)%R(gpid(2,s),gp2,2)
         start(ip0) = 1
         ext(ip1) = sz(s)%R(gpid(1,s),gp1,1)
         ext(ip2) = sz(s)%R(gpid(2,s),gp2,2)
      else
         print *,'Currently cannot handle more than 2-dim. proc. grids'
      endif


      return
      end

      function proc_ndim(dom)

      integer proc_ndim,dom(3),a,i

c      PRINT *,'DOM=',DOM
      a = 0
      do i=1,3
         if(dom(i) .ne. 1) then
            a = a + 1
         endif
      enddo

c      PRINT *,'a=',a
      if(a .eq. 0) then
         a = 1
      endif

      proc_ndim = a
      return 

      end

