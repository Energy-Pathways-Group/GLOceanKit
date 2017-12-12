      program test1
      implicit none
      include 'mpif.h'

      integer, parameter    :: Nx=1025, Ny=8, Nz=513, p1=2
      logical, parameter    :: YZ=.FALSE., XY=.TRUE.   ! only 1 true

#ifdef DOUBLE_PREC
      real(kind=8), allocatable :: A(:,:,:),B(:,:,:), C(:,:,:)
#else
      real(kind=4), allocatable :: A(:,:,:),B(:,:,:), C(:,:,:)
#endif
      
      integer                   :: i,j,k,ii,plan,myid,ierr
      integer                   :: Np,p2, gdims(3)
      integer                   :: plan1,plan2
      integer                   :: order1(3),order2(3)
      integer                   :: layout1(3),layout2(3)
      integer                   :: procs(3)
      integer                   :: ld_in(3),ld_out(3)
      integer                   :: st_in(3),st_out(3)
      logical                   :: cor

      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD,Np,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
      p2 = Np / p1
      
      if(myid==0) print *,'Running on ',Np,' tasks, ',p1,' x ',p2
      

c------------------------------------------------------      
c     YBLOCK to ZBLOCK  (and back)
c------------------------------------------------------
      if( YZ ) then
       layout1(1) = p1
       layout1(2) = 1
       layout1(3) = p2

       layout2(1) = p1
       layout2(2) = p2
       layout2(3) = 1

       order1(1) = 2
       order1(2) = 1
       order1(3) = 3

       order2(1) = 2
       order2(2) = 3
       order2(3) = 1

c------------------------------------------------------      
c     XBLOCK to YBLOCK  (and back)
c------------------------------------------------------ 
      elseif( XY ) then
       layout1(1) = 1
       layout1(2) = p1
       layout1(3) = p2

       layout2(1) = p1
       layout2(2) = 1
       layout2(3) = p2

       order1(1) = 1
       order1(2) = 2
       order1(3) = 3

       order2(1) = 2
       order2(2) = 1
       order2(3) = 3
      endif

c------------------------------      
c     global dims, both cases
c------------------------------
      gdims(1) = Nx
      gdims(2) = Ny
      gdims(3) = Nz
 
      if(myid==0) print *,'Initialize the library; get the plan number'
      call par_transpose2D_init(gdims,layout1,layout2,order1,order2,plan1)
      if(myid==0) print *,'Get local array dimensions for input and output'
      call get_dims(layout1,order1,ld_in,st_in,procs,plan1)
      call get_dims(layout2,order2,ld_out,st_out,procs,plan1)
      print *,'Processor ',myid,': local input array starts at ',st_in,', dimensions ',ld_in
      print *,'Processor ',myid,': local output array starts at ',st_out,', dimensions ',ld_out
      if(myid==0) print *,'Allocate arrays'
      allocate(A(0:ld_in(1)-1,0:ld_in(2)-1,0:ld_in(3)-1))
      allocate(B(0:ld_out(1)-1,0:ld_out(2)-1,0:ld_out(3)-1))
      allocate(C(0:ld_in(1)-1,0:ld_in(2)-1,0:ld_in(3)-1))
      if(myid==0) print *,'Initialize input array      '
      do k=0,ld_in(3)-1
         do j=0,ld_in(2)-1
            do i=0,ld_in(1)-1
               A(i,j,k) = 10000.0 * (k+st_in(3)) + 100.0 * (j+st_in(2)) +i+st_in(1)
            enddo
         enddo
      enddo
      
      if(myid==0) print *,'Call transpose; plan=',plan1 
      call par_transpose2D(A,B,gdims,order1,order2,layout1,layout2,plan1)

#ifdef VERBOSE_MAIN
      print *,'Processor ',myid,' Input array: '
      do k=0,ld_in(3)-1
         do j=0,ld_in(2)-1
            do i=0,ld_in(1)-1,8
               print '(8F9.0)',(A(ii,j,k),ii=i,i+7)
            enddo
            print *,'N'
         enddo
      enddo      

      print '(((32F9.0)))',(((A(i,j,k),i=0,ld_in(1)-1),j=0,ld_in(2)-1),k=0,ld_in(3)-1)
      print *,'Processor ',myid,' Output array: '
      do k=0,ld_out(3)-1
         do j=0,ld_out(2)-1
            do i=0,ld_out(1)-1,8
               print '(8F9.0)',(B(ii,j,k),ii=i,i+7)
            enddo
            print *,'N'
         enddo
      enddo

      print '(((32F9.0)))',(((B(i,j,k),i=0,ld_out(1)-1),j=0,ld_out(2)-1),k=0,ld_out(3)-1)
#endif


      if(myid==0) print *,'Transposing back to original shape'
      print *,'Layout in:',layout2,'Layout_out:',layout1
      call par_transpose2D_init(gdims,layout2,layout1,order2,order1,plan2)
      call par_transpose2D(B,C,gdims,order2,order1,layout2,layout1,plan2)

      print *,'Checking for correctness'
      cor = .true.
      do k=0,ld_in(3)-1
         do j=0,ld_in(2)-1
            do i=0,ld_in(1)-1
               if(abs(A(i,j,k) - C(i,j,k)) .gt. 0.00001) then
                  cor = .false. 
                  goto 11
               endif
            enddo
         enddo
      enddo

11    continue

      if(cor) then
         print *,'Correctness test passed',myid
      else
         print *,'Correctness test NOT passed'
#ifdef VERBOSE_MAIN
         do k=0,ld_in(3)-1
            do j=0,ld_in(2)-1
               do i=0,ld_in(1)-1,8
                  print '(8F9.0)',(C(ii,j,k),ii=i,i+7)
               enddo
               print *,'N'
            enddo
         enddo
#endif
      endif

      print *,'Finished, myid:',myid
      call mpi_finalize(ierr)

      end
