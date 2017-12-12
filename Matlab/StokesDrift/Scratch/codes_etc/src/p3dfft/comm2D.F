C By convention, the local dimension on output has stride-1 storage

C mem_order(1) = the fortran array index corresponding to x
C mem_order(2) = the fortran array index corresponding to y
C mem_order(3) = the fortran array index corresponding to z

      subroutine par_transpose_init(gdims,layout_in,layout_out,mem_order,new_order, plan)

      use comm
      implicit none

      integer p1,p2,gdims(3),mem_order(3),new_order(3),layout_in(3),layout_out(3),plan
      integer procs(3),ip1,ip2,inv_order(3),inv_new(3),i

      if(layout_in(1) .lt. layout_in(2)) then
         procs(1) = layout_in(2)
         procs(2) = layout_in(1)
      else
         procs(1) = layout_in(1)
         procs(2) = layout_in(2)
      endif
      if(procs(2) .lt. layout_in(3)) then
         procs(3) = procs(2)
         if(procs(1) .lt. layout_in(3)) then
            procs(2) = procs(1)
            procs(1) = layout_in(3)
         else
            procs(2) = layout_in(3)
         endif
      else
         procs(3) = layout_in(3)
      endif

C      call mpisetup(procs,gdims,plan)

      return
      end



      subroutine par_transpose2D_init(gdims,layout_in,layout_out,mem_order,new_order, plan)

      use comm
      implicit none

      integer p1,p2,gdims(3),mem_order(3),new_order(3),layout_in(3),layout_out(3),plan
      integer procs_in(3),procs_out(3),ip1,ip2,inv_order(3),inv_new(3),i

C Extract p1 and p2 ?

      if(layout_in(1) .eq. layout_out(1)) then
         if(layout_in(2) .eq. layout_out(2)) then
            return
         else
            p1 = layout_in(1)
            if(layout_in(2) .eq. 1) then
               p2 = layout_in(3)
               ip1 = mem_order(2)
               ip2 = mem_order(3)
            else
               p2 = layout_in(2)
               ip1 = mem_order(3)
               ip2 = mem_order(2)
            endif
         endif
      else if(layout_in(2) .eq. layout_out(2)) then
         print *,'Error: par_transpose2D cannot transpose X <-> Z'
         return
      else if(layout_in(3) .eq. layout_out(3)) then
         if(layout_in(1) .eq. 1) then
            p1 = layout_in(2)
            ip1 = mem_order(1)
            ip2 = mem_order(2)
         else 
            p1 = layout_in(1)
            ip1 = mem_order(2)
            ip2 = mem_order(1)
         endif
         p2 = layout_in(3)
      endif

      if(p1 .eq. 1) then
         p1 = p2
         p2 = 1
      endif

      do i=1,3
         inv_order(mem_order(i)) = i
      enddo

c      print *,'p1,p2: ',p1,p2

c      call mpisetup_2D(p1,p2,gdims(inv_order(:)),plan)
      call mpisetup_2D(p1,p2,gdims,plan)

      return
      end



C --------------------------------------------------------------------
      subroutine par_transpose2D(IN,OUT,gdims,mem_order,new_order,layout_in,
     &   layout_out, plan)

      use comm
      implicit none

      integer plan,init,mpcomm,IERR
      integer gdims(3),mem_order(3),new_order(3),procs_in(3),procs_out(3)
c      complex(mytype) :: IN(n_in(1),n_in(2),n_in(3)),OUT(n_out(1),n_out(2),n_out(3))
      real(mytype) :: IN(*),OUT(*)
      integer i,j,k,p,p1,p2,layout_in(3),layout_out(3),inv_order_in(3),inv_order_out(3)
      integer dp(3),nu(3),nl(3),ip1,ip2,l_in(3),l_out(3),st_in(3),st_out(3)
      type(count) mpichunk
      integer phys_coord1,phys_coord2

      p = layout_in(1) *layout_in(2)*layout_in(3)
      p1 = layout_out(1)*layout_out(2)*layout_out(3)
      if(p .ne. p1) then
         print *,'Error: different number of CPUs for input and output',p,p1
         return
      endif

c      print *,'Start x 1: ',st(plan)%R(:,1,1)
      
c      procs_in(mem_order(:)) = layout_in(:)
c      procs_out(new_order(:)) = layout_out(:)

      call get_dims(layout_in,mem_order,l_in,st_in,dp,plan)
      call get_dims(layout_out,new_order,l_out,st_out,dp,plan)

C Check to see that the setup fits the problem
c      nu = ubound(IN)
c      nl = lbound(IN)
c      if(nu-nl .ne. l_in -1) then
c         print *,'Error in input array dimensions in par_transpose2D: expected ',l_in,' got ',nu-nl+1
c      endif
c      nu = ubound(OUT)
c      nl = lbound(OUT)
c      if(nu-nl .ne. l_in -1) then
c         print *,'Error in output array dimensions in par_transpose2D: expected ',l_out,' got ',nu-nl+1
c      endif

      if(layout_in(1) .eq. layout_out(1)) then
         if(layout_in(2) .eq. layout_out(2)) then
            return
         else
            p1 = layout_in(1)
            mpcomm = mpi_comm_col(plan)
            if(layout_in(2) .eq. 1) then
               p2 = layout_in(3)
               ip1 = mem_order(2)
               ip2 = mem_order(3)
c               mpichunk = Exch(2,3,plan)
               phys_coord1 = 2
               phys_coord2 = 3
            else
               p2 = layout_in(2)
               ip1 = mem_order(3)
               ip2 = mem_order(2)
c                mpichunk = Exch(3,2,plan)
               phys_coord1 = 3
               phys_coord2 = 2
            endif
         endif
      else if(layout_in(2) .eq. layout_out(2)) then
         print *,'Error: par_transpose2D cannot transpose X <-> Z'
         return
      else if(layout_in(3) .eq. layout_out(3)) then
         if(layout_in(1) .eq. 1) then
            p1 = layout_in(2)
            ip1 = mem_order(1)
            ip2 = mem_order(2)
c             mpichunk = Exch(1,2,plan)
               phys_coord1 = 1
               phys_coord2 = 2
         else 
            p1 = layout_in(1)
            ip1 = mem_order(2)
            ip2 = mem_order(1)
c            mpichunk = Exch(2,1,plan)
               phys_coord1 = 2
               phys_coord2 = 1
         endif
         p2 = layout_in(3)
         mpcomm = mpi_comm_row(plan)
      endif

      if(dp(1) .ne. p1 .or. dp(2) .ne. p2) then
         print *,'Error in par_transpose2D: incorrect (p1,p2): expected ',dp,' got ',p1,p2
      endif

      do i=1,3
         inv_order_in(mem_order(i)) = i
         inv_order_out(new_order(i)) = i
      enddo
      

      if(ip1 .eq. 1) then
         if(ip2 .eq. 2) then
            ! Transpose in the plane of first and second dimensions
            call comm_12(IN,OUT,l_in,l_out,p1,p2,new_order(inv_order_in(:)),
     &          phys_coord1,phys_coord2,mpcomm,plan)
         else
            ! Transpose in the plane of first and third dimensions
            call comm_13(IN,OUT,l_in,l_out,p1,p2,new_order(inv_order_in(:)),
     &          phys_coord1,phys_coord2,mpcomm,plan)
         endif
      else if(ip1 .eq. 2) then
         if(ip2 .eq. 1) then
            ! Transpose in the plane of first and second dimensions
            call comm_21(IN,OUT,l_in,l_out,p1,p2,new_order(inv_order_in(:)),
     &          phys_coord1,phys_coord2,mpcomm,plan)
         else
            ! Transpose in the plane of second and third dimensions
            call comm_23(IN,OUT,l_in,l_out,p1,p2,new_order(inv_order_in(:)),
     &          phys_coord1,phys_coord2,mpcomm,plan)
         endif
      else 
         if(ip2 .eq. 1) then
            ! Transpose in the plane of first and third dimensions
            call comm_31(IN,OUT,l_in,l_out,p1,p2,new_order(inv_order_in(:)),
     &          phys_coord1,phys_coord2,mpcomm,plan)
         else
            ! Transpose in the plane of second and third dimensions
            call comm_32(IN,OUT,l_in,l_out,p1,p2,new_order(inv_order_in(:)),
     &          phys_coord1,phys_coord2,mpcomm,plan)
         endif
      endif

      
      return
      end


C ------------------------------------------------------------
      subroutine comm_12(IN,OUT,l_in,l_out,p1,p2,reorder,phys_coord1,phys_coord2,mpcomm,s)

      use comm
      implicit none

      integer p1,p2,pt,ind,l_in(3),l_out(3),s,reorder(3), mpcomm
      integer phys_coord1,phys_coord2,s1,s2
      real(mytype) IN(l_in(1),l_in(2),l_in(3)),OUT(l_out(1),l_out(2),l_out(3))
      integer i,j,x,y,z,position,gran,zp,zen,ierr
      real(mytype), allocatable :: sendbuf(:),recvbuf(:)
      type(count) mpichunk

      if( .not. mpi_set(s)) then
         print *,'comm_12 error: mpi communicators not set'
      endif

      gran = l_in(3) / gblock
      if(gran .lt. 1) then
         gran = 1
      endif

      if(l_in(3) .ne. l_out(reorder(3))) then
c        .or. l_in(2) * l_in(1) .ne. l_out(reorder(1)) * l_out(reorder(2))) then
         print *,'comm_12: Error in in/out array dimensions: l_in=',l_in,
     &             'l_out=',l_out,'reorder=',reorder
      endif

      mpichunk = Exch(phys_coord1,phys_coord2,s)

      allocate(sendbuf(l_in(1)*l_in(2)*gran))
      allocate(recvbuf(l_out(reorder(1))*l_out(reorder(2))*gran))

#ifdef VERBOSE
      print *,'First dimension (phys. ',phys_coord1,
     &       ') is local; exchanging it with second (phys. ',phys_coord2,')'
      print *,'Dimensions: ',l_in,' on input, ',l_out,' on output'
      print *,'With storage order swap ',reorder
      print *,'gran=',gran
#endif

      if(mpcomm .eq. mpi_comm_row(s)) then
         pt = p1
         ind = 1
#ifdef VERBOSE
         print *,p1,' processors in the row (',p2,' orthogonal)'
#endif
      else if(mpcomm .eq. mpi_comm_col(s)) then
         pt = p2
         ind = 2
#ifdef VERBOSE
         print *,p2,' processors in the column (',p1,' orthogonal)'
#endif
      else
         print *,'Error: wrong value of mpcomm =',mpcomm
         return
      endif

      do zp=1,l_in(3),gran

         zen = zp+ gran -1
         if(zen .gt. l_in(3) - gran) then
            zen = l_in(3)
         endif

         position = 1
#ifdef VERBOSE
c         print *,'Packing; zp=',zp, '; s=',s,', ind=',ind,', x limits: (',
c     &          st(s)%R(0,1,ind),en(s)%R(0,1,ind),') to (',
c     &         st(s)%R(pt-1,1,ind),en(s)%R(pt-1,1,ind),')'
#endif
         do i=0,pt-1
            do j=0,gran-1
               z = zp + j
               do y=1,l_in(2)
                  do x=st(s)%R(i,phys_coord1,ind),en(s)%R(i,phys_coord1,ind)
                     sendbuf(position) = IN(x,y,z)
                     position = position+1
                  enddo
               enddo
            enddo
         enddo
         
c Exchange y-z buffers
         
#ifdef VERBOSE
c         print *,'Sending; zp=',zp,', SndCnts=',mpichunk%RSndCnts(:,ind),', RcvCnts=',mpichunk%RRcvCnts(:,ind)
#endif

#ifdef CHECK
         s1 = 0
         s2 = 0
         do i=0,pt-1
            s1 = mpichunk%RSndCnts(i,ind) + s1
            s2 = mpichunk%RRcvCnts(i,ind) + s2
         enddo
         if(s1 .ne. l_in(2) * l_in(1)*mytype) then
            print *,'comm_12: error in dimensions, ', mpichunk%RSndCnts(0:pt-1,ind)
         endif
         if(s2 .ne. l_out(reorder(2)) * l_out(reorder(1))*mytype) then
            print *,'comm_12: error in dimensions, ', mpichunk%RRcvCnts(0:pt-1,ind)
         endif
#endif

         call mpi_alltoallv(sendbuf,mpichunk%RSndCnts(:,ind)*(zen -zp+1), 
     &        mpichunk%RSndStrt(:,ind)*(zen -zp+1),mpi_byte,
     &        recvbuf,mpichunk%RRcvCnts(:,ind)*(zen -zp+1), 
     &        mpichunk%RRcvStrt(:,ind)*(zen -zp+1),mpi_byte,mpcomm,ierr)
         
         position = 1
         if(reorder(3) .eq. 3) then
            if(reorder(2) .eq. 1) then
               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=1,l_out(2)
                           OUT(y,x,z) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo

            else 
               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=1,l_out(1)
                           OUT(x,y,z) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo

            endif
         else if(reorder(3) .eq. 2) then
            if(reorder(2) .eq. 1) then
               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=1,l_out(3)
                           OUT(y,z,x) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo

            else
               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=1,l_out(1)
                           OUT(x,z,y) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo

            endif
         else
            if(reorder(2) .eq. 2) then
               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=1,l_out(3)
                           OUT(z,y,x) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo

            else
               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=1,l_out(2)
                           OUT(z,x,y) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo
            endif
         endif

      enddo

      deallocate(sendbuf)
      deallocate(recvbuf)
      
      return
      end

C ------------------------------------------------------------
      subroutine comm_21(IN,OUT,l_in,l_out,p1,p2,reorder,phys_coord1,phys_coord2,mpcomm,s)

      use comm
      implicit none

      integer p1,p2,pt,ind,l_in(3),l_out(3),reorder(3),mpcomm,s
      integer phys_coord1,phys_coord2,s1,s2
      real(mytype) IN(l_in(1),l_in(2),l_in(3)),OUT(l_out(1),l_out(2),l_out(3))
      integer i,j,x,y,z,position,gran,zp,zen,ierr
      real(mytype), allocatable :: sendbuf(:),recvbuf(:)
      type(count) mpichunk

      if( .not. mpi_set(s)) then
         print *,'comm_21 error: mpi communicators not set'
      endif

      mpichunk = Exch(phys_coord1,phys_coord2,s)

      gran = l_in(3) / gblock
      if(gran .lt. 1) then
         gran = 1
      endif

      if(l_in(3) .ne. l_out(reorder(3))) then
c          .or. l_in(2) * l_in(1) .ne. l_out(reorder(2)) * l_out(reorder(1))) then
         print *,'comm_21: Error in in/out array dimensions'
      endif

      allocate(sendbuf(l_in(1)*l_in(2)*gran))
      allocate(recvbuf(l_out(reorder(1))*l_out(reorder(2))*gran))

#ifdef VERBOSE
      print *,'Second dimension is local (phys. ',phys_coord1,
     &       '); exchanging it with first (phys. ',phys_coord2,')'
      print *,'Dimensions: ',l_in,' on input, ',l_out,' on output'
      print *,'With storage order swap ',reorder
#endif

      if(mpcomm .eq. mpi_comm_row(s)) then
         pt = p1
         ind = 1
#ifdef VERBOSE
         print *,p1,' processors in the row (',p2,' orthogonal)'
#endif
      else if(mpcomm .eq. mpi_comm_col(s)) then
         pt = p2
         ind = 2
#ifdef VERBOSE
         print *,p2,' processors in the column (',p1,' orthogonal)'
#endif
      else
         print *,'Error: wrong value of mpcomm=',mpcomm
      endif

      do zp=1,l_in(3),gran
         zen = zp+ gran -1
         if(zen .gt. l_in(3) - gran) then
            zen = l_in(3)
         endif
         position = 1
         do i=0,pt-1
            do j=0,gran-1
               z = zp + j
               do y=st(s)%R(i,phys_coord1,ind),en(s)%R(i,phys_coord1,ind)
                  do x=1,l_in(1)
                     sendbuf(position) = IN(x,y,z)
                     position = position + 1
                  enddo
               enddo
            enddo
         enddo
         
#ifdef CHECK
         s1 = 0
         s2 = 0
         do i=0,pt-1
            s1 = mpichunk%RSndCnts(i,ind) + s1
            s2 = mpichunk%RRcvCnts(i,ind) + s2
         enddo
         if(s1 .ne. l_in(2) * l_in(1)*mytype) then
            print *,'comm_21: error in dimensions, ', mpichunk%RSndCnts(0:pt-1,ind)
         endif
         if(s2 .ne. l_out(reorder(2)) * l_out(reorder(1))*mytype) then
            print *,'comm_21: error in dimensions, ', mpichunk%RRcvCnts(0:pt-1,ind)
         endif
#endif

         call mpi_alltoallv(sendbuf,mpichunk%RSndCnts(:,ind)*(zen -zp+1), 
     &        mpichunk%RSndStrt(:,ind)*(zen -zp+1),mpi_byte,
     &        recvbuf,mpichunk%RRcvCnts(:,ind)*(zen -zp+1), 
     &        mpichunk%RRcvStrt(:,ind)*(zen -zp+1),mpi_byte,mpcomm,ierr)
         
         position = 1
         if(reorder(3) .eq. 3) then
            if(reorder(2) .eq. 2) then

               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=1,l_out(2)
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(x,y,z) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            else 
               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=1,l_out(1)
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(y,x,z) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            endif

         else if(reorder(3) .eq. 2) then
            if(reorder(2) .eq. 3) then

               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=1,l_out(3)
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(x,z,y) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            else 
               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=1,l_out(1)
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(y,z,x) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            endif

         else 
            if(reorder(2) .eq. 3) then

               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=1,l_out(3)
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(z,x,y) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            else 
               do i=0,pt-1
                  do j=0,gran-1
                     z = zp + j
                     do y=1,l_out(2)
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(z,y,x) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            endif
         endif

      enddo

      deallocate(sendbuf)
      deallocate(recvbuf)

      return 
      end

C ------------------------------------------------------------
      subroutine comm_13(IN,OUT,l_in,l_out,p1,p2,reorder,phys_coord1,phys_coord2,mpcomm,s)

      use comm
      implicit none

      integer p1,p2,pt,ind,l_in(3),l_out(3),s,reorder(3),mpcomm
      integer phys_coord1,phys_coord2,s1,s2
      real(mytype) IN(l_in(1),l_in(2),l_in(3)),OUT(l_out(1),l_out(2),l_out(3))
      integer i,x,y,z,position,iy,yst,yen,gran,j,ierr
      real(mytype), allocatable :: sendbuf(:),recvbuf(:)
      type(count) mpichunk
                  
      if( .not. mpi_set(s)) then
         print *,'comm_13 error: mpi communicators not set'
      endif

      if(l_in(2) .ne. l_out(reorder(2))) then
c         .or. l_in(1) * l_in(3) .ne. l_out(reorder(1)) * l_out(reorder(3))) then
         print *,'comm_13: Error in in/out array dimensions'
      endif

      mpichunk = Exch(phys_coord1,phys_coord2,s)

      gran = l_in(2) / gblock
      if(gran .lt. 1) then
         gran = 1
      endif

      allocate(sendbuf(l_in(1)*l_in(3)*gran))
      allocate(recvbuf(l_out(reorder(1))*l_out(reorder(3))*gran))

#ifdef VERBOSE
      print *,'First dimension is local (phys. ',phys_coord1,
     &   '); exchanging it with third (phys. ',phys_coord2,')'
      print *,'Dimensions: ',l_in,' on input, ',l_out,' on output'
      print *,'With storage order swap ',reorder
#endif

      if(mpcomm .eq. mpi_comm_row(s)) then
         pt = p1
         ind = 1
#ifdef VERBOSE
         print *,p1,' processors in the row (',p2,' orthogonal)'
#endif
      else if(mpcomm .eq. mpi_comm_col(s)) then
         pt = p2
         ind = 2
#ifdef VERBOSE
         print *,p2,' processors in the column (',p1,' orthogonal)'
#endif
      else
         print *,'Error: wrong value of mpcomm=',mpcomm
      endif

      do iy = 0,gblock
         if(iy .eq. gblock) then
            if(yen .ge. l_in(2)) goto 100
            yst = yen+1
            yen = l_in(2)
         else
            yst = 1 + iy * gran
            yen = yst + gran -1
         endif
         
         position = 1
         do i=0,pt-1
            do z=1,l_in(3)
               do y=yst,yen
                  do x=st(s)%R(i,phys_coord1,ind),en(s)%R(i,phys_coord1,ind)
                     sendbuf(position) = IN(x,y,z)
                     position = position+1
                  enddo
               enddo
            enddo
         enddo

#ifdef CHECK
         s1 = 0
         s2 = 0
         do i=0,pt-1
            s1 = mpichunk%RSndCnts(i,ind) + s1
            s2 = mpichunk%RRcvCnts(i,ind) + s2
         enddo
         if(s1 .ne. l_in(3) * l_in(1)*mytype) then
            print *,'comm_13: error in dimensions, ', mpichunk%RSndCnts(0:pt-1,ind)
         endif
         if(s2 .ne. l_out(reorder(3)) * l_out(reorder(1))*mytype) then
            print *,'comm_13: error in dimensions, ', mpichunk%RRcvCnts(0:pt-1,ind)
         endif
#endif

         call mpi_alltoallv(sendbuf,mpichunk%RSndCnts(:,ind)*(yen-yst+1), 
     &       mpichunk%RSndStrt(:,ind)*(yen-yst+1),mpi_byte,
     &          recvbuf,mpichunk%RRcvCnts(:,ind)*(yen-yst+1),
     &          mpichunk%RRcvStrt(:,ind)*(yen-yst+1),mpi_byte,mpcomm,ierr)

         position = 1
         if(reorder(2) .eq. 1) then
            if(reorder(1) .eq. 2) then
               do i=0,pt-1
                  do z=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                     do y=yst,yen
                        do x=1,l_out(2)
                           OUT(y,x,z) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo               
            else if(reorder(1) .eq. 3) then
               do i=0,pt-1
                  do z=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                     do y=yst,yen
                        do x=1,l_out(3)
                           OUT(y,z,x) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo
            endif
         else  if(reorder(2) .eq. 2) then
            if(reorder(1) .eq. 1) then
               do i=0,pt-1
                  do z=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                     do y=yst,yen
                        do x=1,l_out(1)
                           OUT(x,y,z) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo
            else if(reorder(1) .eq. 3) then
               do i=0,pt-1
                  do z=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                     do y=yst,yen
                        do x=1,l_out(3)
                           OUT(z,y,x) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo
            endif 
         else  
            if(reorder(1) .eq. 1) then
               do i=0,pt-1
                  do z=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                     do y=yst,yen
                        do x=1,l_out(1)
                           OUT(x,z,y) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo               
            else if(reorder(1) .eq. 2) then
               do i=0,pt-1
                  do z=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                     do y=yst,yen
                        do x=1,l_out(2)
                           OUT(z,x,y) = recvbuf(position)
                           position = position + 1
                        enddo
                     enddo
                  enddo
               enddo
            endif
         endif
         
 100     continue
      enddo
   
      deallocate(sendbuf)
      deallocate(recvbuf)
      
      return
      end

C ------------------------------------------------------------
      subroutine comm_31(IN,OUT,l_in,l_out,p1,p2,reorder,phys_coord1,phys_coord2,mpcomm,s)

      use comm
      implicit none

      integer p1,p2,pt,ind,l_in(3),l_out(3),s,reorder(3),mpcomm
      integer phys_coord1,phys_coord2,s1,s2
      real(mytype) IN(l_in(1),l_in(2),l_in(3)),OUT(l_out(1),l_out(2),l_out(3))
      integer i,x,y,z,position,j,yst,yen,gran,ierr,iy
      real(mytype), allocatable :: sendbuf(:),recvbuf(:)
      type(count) mpichunk
    

      if( .not. mpi_set(s)) then
         print *,'comm_31 error: mpi communicators not set'
      endif

      if(l_in(2) .ne. l_out(reorder(2))) then
c         .or. l_in(1) * l_in(3) .ne. l_out(reorder(1)) * l_out(reorder(3))) then
         print *,'comm_31: Error in in/out array dimensions'
      endif

      mpichunk = Exch(phys_coord1,phys_coord2,s)

      gran = l_in(2) / gblock
      if(gran .lt. 1) then
         gran = 1
      endif

      allocate(sendbuf(l_in(1)*l_in(3)*gran))
      allocate(recvbuf(l_out(reorder(1))*l_out(reorder(3))*gran))

#ifdef VERBOSE
      print *,'Third dimension is local (phys. ',phys_coord1,
     &   '); exchanging it with first (phys. ',phys_coord2,')'
      print *,'Dimensions: ',l_in,' on input, ',l_out,' on output'
      print *,'With storage order swap ',reorder
#endif

      if(mpcomm .eq. mpi_comm_row(s)) then
         pt = p1
         ind = 1
#ifdef VERBOSE
         print *,p1,' processors in the row (',p2,' orthogonal)'
#endif
      else if(mpcomm .eq. mpi_comm_col(s)) then
         pt = p2
         ind = 2
#ifdef VERBOSE
         print *,p2,' processors in the column (',p1,' orthogonal)'
#endif
      else
         print *,'Error: wrong value of mpcomm=',mpcomm
      endif
      
      do iy = 0,gblock
         if(iy .eq. gblock) then
            if(yen .ge. l_in(2)) goto 101
            yst = yen+1
            yen = l_in(2)
         else
            yst = 1 + iy * gran
            yen = yst + gran -1
         endif
         
         position = 1
         do i=0,pt-1
            do z=st(s)%R(i,phys_coord1,ind),en(s)%R(i,phys_coord1,ind)
               do y=yst,yen
                  do x=1,l_in(1)
                     sendbuf(position) = IN(x,y,z)
                     position = position + 1
                  enddo
               enddo
            enddo
         enddo
         
#ifdef CHECK
         s1 = 0
         s2 = 0
         do i=0,pt-1
            s1 = mpichunk%RSndCnts(i,ind) + s1
            s2 = mpichunk%RRcvCnts(i,ind) + s2
         enddo
         if(s1 .ne. l_in(3) * l_in(1)*mytype) then
            print *,'comm_31: error in dimensions, ', mpichunk%RSndCnts(0:pt-1,ind)
         endif
         if(s2 .ne. l_out(reorder(3)) * l_out(reorder(1))*mytype) then
            print *,'comm_31: error in dimensions, ', mpichunk%RRcvCnts(0:pt-1,ind)
         endif
#endif

         call mpi_alltoallv(sendbuf,mpichunk%RSndCnts(:,ind)*(yen-yst+1), 
     &        mpichunk%RSndStrt(:,ind)*(yen-yst+1),mpi_byte,
     &        recvbuf,mpichunk%RRcvCnts(:,ind)*(yen-yst+1), 
     &        mpichunk%RRcvStrt(:,ind)*(yen-yst+1),mpi_byte,mpcomm,ierr)
         
         position = 1
         if(reorder(2) .eq. 1) then
            if(reorder(1) .eq. 2) then
               do i=0,pt-1
                  do z=1,l_out(3)
                     do y=yst,yen
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(y,x,z) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            else if(reorder(1) .eq. 3) then
               do i=0,pt-1
                  do z=1,l_out(2)
                     do y=yst,yen
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(y,z,x) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            endif
         else         if(reorder(2) .eq. 2) then
            if(reorder(1) .eq. 1) then
               do i=0,pt-1
                  do z=1,l_out(3)
                     do y=yst,yen
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(x,y,z) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            else if(reorder(1) .eq. 3) then
               do i=0,pt-1
                  do z=1,l_out(1)
                     do y=yst,yen
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(z,y,x) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            endif
         else 
            if(reorder(1) .eq. 2) then
               do i=0,pt-1
                  do z=1,l_out(1)
                     do y=yst,yen
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(z,x,y) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            else if(reorder(1) .eq. 1) then
               do i=0,pt-1
                  do z=1,l_out(2)
                     do y=yst,yen
                        do x=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                           OUT(x,z,y) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            endif
         endif

 101     continue

      enddo
      
      deallocate(sendbuf)
      deallocate(recvbuf)

      return
      end

C ------------------------------------------------------------
      subroutine comm_23(IN,OUT,l_in,l_out,p1,p2,reorder,phys_coord1,phys_coord2,mpcomm,s)

      use comm
      implicit none

      integer p1,p2,pt,ind,l_in(3),l_out(3),s,reorder(3),mpcomm
      integer phys_coord1,phys_coord2,s1,s2
      real(mytype) IN(l_in(1),l_in(2),l_in(3)),OUT(l_out(1),l_out(2),l_out(3))
      integer i,x,y,z,position,j,xst,xen,gran,ierr,ix
      real(mytype), allocatable :: sendbuf(:),recvbuf(:)
      type(count) mpichunk

      if( .not. mpi_set(s)) then
         print *,'comm_23 error: mpi communicators not set'
      endif

      if(l_in(1) .ne. l_out(reorder(1))) then
c          .or. l_in(2) * l_in(3) .ne. l_out(reorder(2)) * l_out(reorder(3))) then
         print *,'comm_23: Error in in/out array dimensions'
      endif

      mpichunk = Exch(phys_coord1,phys_coord2,s)

      gran = l_in(1) / gblock
      if(gran .lt. 1) then
         gran = 1
      endif
      

#ifdef VERBOSE
      print *,'Second dimension is local (phys. ',phys_coord1,
     &   '); exchanging it with third (phys. ',phys_coord2,')'
      print *,'Dimensions: ',l_in,' on input, ',l_out,' on output'
      print *,'With storage order swap ',reorder
#endif

      allocate(sendbuf(gran*l_in(2)*l_in(3)))
      allocate(recvbuf(gran*l_out(reorder(2))*l_out(reorder(3))))

      if(mpcomm .eq. mpi_comm_row(s)) then
         pt = p1
         ind = 1
#ifdef VERBOSE
         print *,p1,' processors in the row (',p2,' orthogonal)'
#endif
      else if(mpcomm .eq. mpi_comm_col(s)) then
         pt = p2
         ind = 2
#ifdef VERBOSE
         print *,p2,' processors in the column (',p1,' orthogonal)'
#endif
      else
         print *,'Error: wrong value of mpcomm=',mpcomm
      endif

      do ix = 0,gblock
         if(ix .eq. gblock) then
            if(xen .ge. l_in(1)) goto 102
            xst = xen+1
            xen = l_in(1)
         else
            xst = 1 + ix * gran
            xen = xst + gran -1
         endif
         
         position =1
         do i=0,pt-1
            do z=1,l_in(3)
               do y=st(s)%R(i,phys_coord1,ind),en(s)%R(i,phys_coord1,ind)
                  do x=xst,xen
                     sendbuf(position) = IN(x,y,z)
                     position = position+1
                  enddo
               enddo
            enddo
         enddo
         
#ifdef CHECK
         s1 = 0
         s2 = 0
         do i=0,pt-1
            s1 = mpichunk%RSndCnts(i,ind) + s1
            s2 = mpichunk%RRcvCnts(i,ind) + s2
         enddo
         if(s1 .ne. l_in(2) * l_in(3) * mytype) then
            print *,'comm_23: error in IN dimensions, ',mpichunk%RSndCnts(0:pt-1,ind)
         endif
         if(s2 .ne. l_out(reorder(2)) * l_out(reorder(3)) *mytype) then
            print *,'comm_23: error in OUT dimensions,', mpichunk%RRcvCnts(0:pt-1,ind)
         endif
#endif

         call mpi_alltoallv(sendbuf,mpichunk%RSndCnts(:,ind)*(xen-xst+1), 
     &        mpichunk%RSndStrt(:,ind)*(xen-xst+1),mpi_byte,
     &        recvbuf,mpichunk%RRcvCnts(:,ind)*(xen-xst+1), 
     &        mpichunk%RRcvStrt(:,ind)*(xen-xst+1),mpi_byte,mpcomm,ierr)

         position = 1
         if(reorder(1) .eq. 1) then
            if(reorder(2) .eq. 2) then
               do z=1,l_out(3)
                  do y=1,l_out(2)
                     do x=xst,xen
                        OUT(x,y,z) = recvbuf(position)
                        position = position + 1
                     enddo
                  enddo
               enddo
            else if(reorder(2) .eq. 3) then
               do z=1,l_out(2)
                  do y=1,l_out(3)
                     do x=xst,xen
                        OUT(x,z,y) = recvbuf(position)
                        position = position + 1
                     enddo
                  enddo
               enddo
            endif
         else          if(reorder(1) .eq. 2) then
            if(reorder(2) .eq. 1) then
               do z=1,l_out(3)
                  do y=1,l_out(1)
                     do x=xst,xen
                        OUT(y,x,z) = recvbuf(position)
                        position = position + 1
                     enddo
                  enddo
               enddo
            else if(reorder(2) .eq. 3) then
               do z=1,l_out(1)
                  do y=1,l_out(3)
                     do x=xst,xen
                        OUT(z,x,y) = recvbuf(position)
                        position = position + 1
                     enddo
                  enddo
               enddo
            endif
         else
            if(reorder(2) .eq. 2) then
               do z=1,l_out(1)
                  do y=1,l_out(2)
                     do x=xst,xen
                        OUT(z,y,x) = recvbuf(position)
                        position = position + 1
                     enddo
                  enddo
               enddo
            else if(reorder(2) .eq. 1) then
               do z=1,l_out(2)
                  do y=1,l_out(1)
                     do x=xst,xen
                        OUT(y,z,x) = recvbuf(position)
                        position = position + 1
                     enddo
                  enddo
               enddo
            endif
         endif

 102     continue
      enddo

      deallocate(sendbuf)
      deallocate(recvbuf)
         
      return
      end

C ------------------------------------------------------------
      subroutine comm_32(IN,OUT,l_in,l_out,p1,p2,reorder,phys_coord1,phys_coord2,mpcomm,s)

      use comm
      implicit none

      integer p1,p2,pt,ind,l_in(3),l_out(3),s,reorder(3),mpcomm
      integer phys_coord1,phys_coord2,s1,s2
      real(mytype) IN(l_in(1),l_in(2),l_in(3)),OUT(l_out(1),l_out(2),l_out(3))
      integer i,x,y,z,position,j,xst,xen,gran,ierr,ix
      real(mytype), allocatable :: sendbuf(:),recvbuf(:)
      type(count) mpichunk

      if( .not. mpi_set(s)) then
         print *,'comm_32 error: mpi communicators not set'
      endif

      if(l_in(1) .ne. l_out(reorder(1))) then 
c           .or. l_in(2) * l_in(3) .ne. l_out(reorder(2)) * l_out(reorder(3))) then
         print *,'comm_32: Error in in/out array dimensions'
      endif

      mpichunk = Exch(phys_coord1,phys_coord2,s)

      gran = l_in(1) / gblock
      if(gran .lt. 1) then
         gran = 1
      endif
      

C third dimension is local, to be exchanged with second

#ifdef VERBOSE
      print *,'Third dimension is local (phys. ',phys_coord1,
     &   '); exchanging it with second (phys. ',phys_coord2,')'
      print *,'Dimensions: ',l_in,l_out
      print *,'With storage order swap ',reorder
#endif

      allocate(sendbuf(gran*l_in(2)*l_in(3)))
      allocate(recvbuf(gran*l_out(reorder(2))*l_out(reorder(3))))

      if(mpcomm .eq. mpi_comm_row(s)) then
         pt = p1
         ind = 1
#ifdef VERBOSE
         print *,p1,' processors in the row (',p2,' orthogonal)'
#endif
      else if(mpcomm .eq. mpi_comm_col(s)) then
         pt = p2
         ind = 2
#ifdef VERBOSE
         print *,p2,' processors in the column (',p1,' orthogonal)'
#endif
      else
         print *,'Error: wrong value of mpcomm=',mpcomm
      endif

      do ix = 0,gblock
         if(ix .eq. gblock) then
            if(xen .ge. l_in(1)) goto 103
            xst = xen+1
            xen = l_in(1)
         else
            xst = 1 + ix * gran
            xen = xst + gran -1
         endif
         
         position = 1
         do z=1,l_in(3)
            do y=1,l_in(2)
               do x=xst,xen
                  sendbuf(position) = IN(x,y,z)
                  position = position + 1
               enddo
            enddo
         enddo
         
#ifdef CHECK
         s1 = 0
         s2 = 0
         do i=0,pt-1
            s1 = mpichunk%RSndCnts(i,ind) + s1
            s2 = mpichunk%RRcvCnts(i,ind) + s2
         enddo
         if(s1 .ne. l_in(2) * l_in(3)*mytype) then
            print *,'comm_32: error in dimensions, ', mpichunk%RSndCnts(0:pt-1,ind)
         endif
         if(s2 .ne. l_out(reorder(2)) * l_out(reorder(3))*mytype) then
            print *,'comm_32: error in dimensions, ', mpichunk%RRcvCnts(0:pt-1,ind)
         endif
#endif

         call mpi_alltoallv(sendbuf,mpichunk%RSndCnts(:,ind)*(xen-xst+1), 
     &        mpichunk%RSndStrt(:,ind)*(xen-xst+1),mpi_byte,
     &        recvbuf,mpichunk%RRcvCnts(:,ind)*(xen-xst+1), 
     &        mpichunk%RRcvStrt(:,ind)*(xen-xst+1),mpi_byte,mpcomm,ierr)

c Unpack receive buffers into XZgY

         position=1
         if(reorder(1) .eq. 1) then
            if(reorder(2) .eq. 2) then
               do i=0,pt-1
                  do z=1,l_out(3)
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=xst,xen
                           OUT(x,y,z) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            else if(reorder(2) .eq. 3) then
               do i=0,pt-1
                  do z=1,l_out(2)
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=xst,xen
                           OUT(x,z,y) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            endif
         else if(reorder(1) .eq. 2) then
            if(reorder(2) .eq. 1) then
               do i=0,pt-1
                  do z=1,l_out(3)
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=xst,xen
                           OUT(y,x,z) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            else if(reorder(2) .eq. 3) then
               do i=0,pt-1
                  do z=1,l_out(1)
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=xst,xen
                           OUT(z,x,y) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            endif
         else
            if(reorder(2) .eq. 2) then
               do i=0,pt-1
                  do z=1,l_out(1)
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=xst,xen
                           OUT(z,y,x) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            else if(reorder(2) .eq. 1) then
               do i=0,pt-1
                  do z=1,l_out(2)
                     do y=st(s)%R(i,phys_coord2,ind),en(s)%R(i,phys_coord2,ind)
                        do x=xst,xen
                           OUT(y,z,x) = recvbuf(position)
                           position = position+1
                        enddo
                     enddo
                  enddo
               enddo
            endif
         endif

 103     continue
      enddo

      deallocate(sendbuf)
      deallocate(recvbuf)
         
      return
      end

