#=====================================================	
#  FlowSolve Makefile for src directory
#=====================================================	
SHELL=/bin/sh
include Make.inc

LIBS= ${ALL_EXTERNAL_LIBS}

FS_SRCS=data_modules.f90 user_params_module.f90 flow_solve.f90 initialize.f90 mpi_utils.f90 read_userdata.f90 \
		flow_solve_utils.f90 decomposition2D.f90 preliminary_tasks.f90 setup_domain.f90 \
		flow_solve_user.f90 setup_derivs.f90 deriv_toolbox.f90 differentiate.f90 \
		transpose_routines.f90 allocate_memory.f90 setup_diffusion.f90 explicit_rhs.f90 \
		explicit_time_step.f90 eos.f90 initial_conditions.f90 pressure_projection.f90 \
		transform_xy.f90 z_poisson_solver.f90 deboor.f machine.f90 \
		write_results.f90  interpolation_routines.f90 particle_routines.f90 \
		setup_immersed_boundary.f90 immersed_boundary.f90 apply_forcing.f90 newton.f \
		

FS_OBJS=data_modules.o user_params_module.o flow_solve.o initialize.o mpi_utils.o read_userdata.o \
		flow_solve_utils.o decomposition2D.o preliminary_tasks.o setup_domain.o \
		flow_solve_user.o setup_derivs.o deriv_toolbox.o differentiate.o \
		transpose_routines.o allocate_memory.o setup_diffusion.o explicit_rhs.o \
		explicit_time_step.o eos.o initial_conditions.o pressure_projection.o \
		transform_xy.o z_poisson_solver.o  \
		deboor.o machine.o interpolation_routines.o particle_routines.o write_results.o \
		setup_immersed_boundary.o immersed_boundary.o apply_forcing.o newton.o \
		p3dfft/module.o p3dfft/setup.o p3dfft/comm2D.o 
		

DRIVER_SRCS=test_driver.f90 data_modules.f90 
			
DRIVER_OBJS=test_driver.o data_modules.o 


#==================
# Targets
#==================

all:
	@echo
	make flow.x
	mv   flow.x ${BASE}	
	
flow.x:	$(FS_SRCS) $(FS_OBJS)
	${FLINKER} ${LDFLAGS} -o $@ ${FS_OBJS} $(LIBS)

test_driver.x: ${DRIVER_OBJS}
	 ${FLINKER} ${LDFLAGS} -o $@ ${DRIVER_OBJS} $(LIBS)
	        
clean:
	@echo
	@echo "src:  removing  *.o, *.mod,  Make.inc, flow_solve_user.f90 "
	rm -f *.o *.mod  Make.inc flow_solve_user.f90 *.x fort*

#=============================================================
# Dependencies  
#=============================================================
flow.x:	$(FS_OBJS)





#=============================================================
# Suffix Rules
# Gfortran creates a file name.mod for each module name. 
# This confuses the gnu make, which thinks that .mod files
# are Modula2 source code files. Override w/ %.o : %.mod
#=============================================================

.SUFFIXES: $(SUFFIXES) .f90  .F90 .f .mod

.f.o:
	${F77} $(DFLAGS) $(F77OPTFLAGS) -c $<

.f90.o:
	${F90} $(DFLAGS) $(FOPTFLAGS) $(ALL_INCLUDES) -c $<
	

.f90.mod:
	${F90} $(FOPTFLAGS) $(ALL_INCLUDES) -c $<
	

%.o : %.mod
