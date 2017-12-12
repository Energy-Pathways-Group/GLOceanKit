SHELL=/bin/sh
VENDOR=gcc



#--------------------------------------------------------------------
#    Library locations ... OSX 10.6.6 MacPorts gfortran from gcc4.4
#--------------------------------------------------------------------
#MPI_ROOT= /usr/mpi/gcc/mvapich2-1.0.3
#MPI_ROOT = /usr
MPI_ROOT = /usr/lib64/openmpi
NETCDF_ROOT= /usr
FFTW3_ROOT= /usr

#VENDOR=portland
#MPI_ROOT= /usr/local/pgi/mpich2
#NETCDF_ROOT= /usr/local/pgi/netcdf
#FFTW3_ROOT= /usr/local/pgi/fftw3
LD_LIBRARY_PATH =  # only search for explicitly indicated libraries


#--------------------------------------------------------------------
#    Run and output base directory locations
#--------------------------------------------------------------------
BASE=${PWD}


#--------------------------------------------------------------------
#    Compiler wrappers for MPI
#--------------------------------------------------------------------
CC = ${MPI_ROOT}/bin/mpicc
F77 = ${MPI_ROOT}/bin/mpif77 
F90 = ${MPI_ROOT}/bin/mpif90 
FLINKER = ${F90}

 
#--------------------------------------------------------------------
#    Vendor specific compilation parameters
#--------------------------------------------------------------------

ifeq ($(VENDOR), gcc)  
    COMPILERS = /usr/bin
	LDFLAGS = 
	FOPTFLAGS   = -g -O2  -ffree-form -ffree-line-length-none
	F77OPTFLAGS = -g -O2  -ffixed-form 
	BLAS_LIB=  /usr/lib64/libblas.a
	LAPACK_LIB = /usr/lib64/liblapack.a
endif


ifeq ($(VENDOR), portland)
    COMPILERS = /opt/pgi/osx86-64/7.1-4/bin
	LDFLAGS = 
	FOPTFLAGS   = -Mpreprocess -g -O2 -Mvect=sse -tp p7-64 
	F77OPTFLAGS = -Mpreprocess -g -O2 -Mvect=sse -tp p7-64 
	BLAS_LIB	= -L/opt/pgi/osx86-64/7.0/lib -lblas
	LAPACK_LIB	= -L/opt/pgi/osx86-64/7.0/lib -llapack
endif

ifeq ($(VENDOR), intel) 
    FCOMPILERS = /opt/intel/fce/10.1.006/bin
    CCOMPILERS = /opt/intel/cce/10.1.006/bin
    LDFLAGS = -openmp
	LD_LIBRARY_PATH = /opt/intel/fce/10.1.006/lib:/opt/intel/cce/10.1.006/lib
 	DYLD_LIBRARY_PATH = /opt/intel/fce/10.1.006/lib:/opt/intel/cce/10.1.006/lib
 	FOPTFLAGS =  -O2 -g -fpp -xP -free 
 	F77OPTFLAGS =  -O2 -g -fpp -xP -fixed
 	BLAS_LIB=
 	LAPACK_LIB=	-framework VecLib
endif


AR = /usr/bin/ar 	 
AR_FLAGS =	cru

#--------------------------------------------------------------------
#  External libraries
#  (no explicit mpi  ==> using mpif77, mpif90)
#--------------------------------------------------------------------  
NETCDF_INC=	-I${NETCDF_ROOT}/include
NETCDF_LIB=	-L${NETCDF_ROOT}/lib64 -lnetcdff
ifeq ($(VENDOR), gcc)
	NETCDF_LIB=	-L${NETCDF_ROOT}/lib64 -lnetcdff
endif
FFTW_INC=	    -I${FFTW3_ROOT}/include
FFTW_LIB=	    -L${FFTW3_ROOT}/lib64 -lfftw3

ALL_EXTERNAL_LIBS=	$(NETCDF_LIB) $(LAPACK_LIB) $(BLAS_LIB) \
 					$(FFTW_LIB)
 					
ALL_INCLUDES=  	$(NETCDF_INC) $(FFTW_INC)	

 
 
 		 
#--------------------------------------------------------------------
showconfig:
	@echo
	-@echo  "host:                 " `hostname -s`
	-@echo  "                      " `uname -sr`
	-@echo  "user:                 " `whoami`
	-@echo  "compile date:         " `date -u +"%Y-%m-%dT%H:%M:%SZ"`
	-@echo  "VENDOR:               " ${VENDOR}
	-@echo  "LD_LIBRARY_PATH:      " ${LD_LIBRARY_PATH}
	-@echo  "DYLD_LIBRARY_PATH:    " ${DYLD_LIBRARY_PATH}
	-@echo  "BASE:                 " ${BASE}
	-@echo  "OUTPUT_ROOT:          " ${BASE}/output
	-@echo  "F90:                  " ${F90}
	-@echo  "F77:                  " ${F77}
	-@echo  "DFLAGS:               " ${DFLAGS}
	-@echo  "FOPTFLAGS:            " ${FOPTFLAGS}
	-@echo  "F77OPTFLAGS:          " ${F77OPTFLAGS}
	-@echo  "FLINKER:              " ${FLINKER}
	-@echo  "NETCDF_INC:           " ${NETCDF_INC}
	-@echo  "NETCDF_LIB:           " ${NETCDF_LIB}
	-@echo  "FFTW_INC:             " ${FFTW_INC}
	-@echo  "FFTW_LIB:             " ${FFTW_LIB}
	-@echo  "LAPACK_LIB:           " ${LAPACK_LIB}
	-@echo  "BLAS_LIB:             " ${BLAS_LIB}
	-@echo  "TRANSPOSE_LIB:        " ${TRANSPOSE_LIB}
	-@echo  "ALL_EXTERNAL_LIBS:    " ${ALL_EXTERNAL_LIBS}
	-@echo  "ALL_INCLUDES:         " ${ALL_INCLUDES}
	-@echo  "XTRA_LIBS:            " ${XTRA_LIBS}
