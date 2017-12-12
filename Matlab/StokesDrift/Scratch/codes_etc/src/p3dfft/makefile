include Make.inc
FF = ${MPI_ROOT}/bin/mpif90 

ifeq ($(VENDOR), gcc)
 FFLAGS = -g -O3 -ffast-math -ffixed-line-length-132
endif

INCL= 
LDFLAGS=
DFLAGS = "-DDOUBLE_PREC"


P3DFFT_ROOT = .
P3DFFT_LIB = ${P3DFFT_ROOT}/libp3dfft.3.a


FFT3DLIB = module.o setup.o comm2D.o

 

#all: lib    don't actually make a library...
all: $(FFT3DLIB)

lib: $(FFT3DLIB)
	$(AR) $(ARFLAGS) $(P3DFFT_LIB) $(FFT3DLIB)
	
test: $(FFT3DLIB) test1_2D.o
	$(FF) test1_2D.o -L. -lp3dfft.3 $(LDFLAGS) 

setup.o: setup.f module.o
comm2D.o: comm2D.F module.o

.c.o: 
		$(CC) -c $(CFLAGS) *.c
.F.o:
		$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.f.o: 
		$(FF) -c $(DFLAGS) $(FFLAGS) $(INCL) $<
clean:
	/bin/rm -f libp3dfft.3.a $(FFT3DLIB) *.mod *.i a.out test1_2D.o *~
	
