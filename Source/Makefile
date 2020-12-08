FCFLAGS= -O1  -funroll-loops -p -g 
BUILD_DIR ?= ./../build
current_dir = '"'$(shell pwd)'"'

ifeq ($(BUILD),fast)
OPT= '"FAST"'
ifeq ($(F90),gfortran)
FCFLAGS= -fconvert=big-endian -fno-realloc-lhs -fopenmp -fPIC -O3 -funroll-loops -fno-signed-zeros -g -fbacktrace -march=native -mno-avx -fbackslash
endif
ifeq ($(F90),ifort)
FCFLAGS=-convert big_endian -assume norealloc_lhs -diag-disable 8290 -diag-disable 8291 -qopenmp -fPIC -O3 -debug minimal -traceback -xHost -fbackslash
endif

endif

ifeq ($(BUILD),debug)
OPT= '"DEBUG"'
ifeq ($(F90),gfortran)
FCFLAGS= -static-libgfortran -Ddebug -O0 -g -fbounds-check -fbacktrace -Wall -Waliasing -Wsurprising -Wline-truncation -Wno-tabs -Wno-uninitialized -Wno-unused-dummy-argument -Wno-unused -Wno-character-truncation -Wl,-no_pie -fbackslash
endif
ifeq ($(F90),ifort)
FCFLAGS= -convert big_endian -assume norealloc_lhs -diag-disable 8290 -diag-disable 8291 -qopenmp -fPIC -Ddebug -O0 -g -debug extended -traceback -C -stand f03 -warn all -diag-disable 7025 -diag-disable 7712 -diag-disable 5268 -diag-disable 7847 -diag-disable 7893 -diag-disable 6439 -fbackslash
endif
endif

#FCFLAGS= -O1  -funroll-loops -p -g

ifeq ($(COMMS_ARCH),mpi)

n_body.mpi: n_body.o io.o trace.o comms_mpi.o pot.o diff.o
	$(MPI_F90) -o n_body.mpi  n_body.o io.o trace.o comms_mpi.o pot.o diff.o $(FCFLAGS)
# Here are the compile steps

n_body.o:n_body.f90 io.o trace.o comms_mpi.o pot.o diff.o
	$(MPI_F90) -c  n_body.f90 $(FCFLAGS)

io.o:./io.f90 trace.o comms_mpi.o
	@echo $(OPT)
	$(MPI_F90) -c  ./io.f90 -cpp -D life_dir=$(current_dir) -D opt_strat=$(OPT) $(FCFLAGS)

trace.o:./trace.f90
	$(MPI_F90) -c  ./trace.f90 $(FCFLAGS)

comms_mpi.o:./comms_mpi.f90 trace.o
	$(MPI_F90) -c  ./comms_mpi.f90 $(FCFLAGS)

pot.o:./pot.f90 comms_mpi.o io.o trace.o
	$(MPI_F90) -c ./pot.f90 $(FCFLAGS)
diff.o:./diff.f90 comms_mpi.o io.o trace.o pot.o
	$(MPI_F90) -c ./diff.f90 $(FCFLAGS)



#life.o:./life.f90 comms_mpi.o io.o trace.o
#	$(MPI_F90) -c  ./life.f90 $(FCFLAGS)


endif


ifeq ($(COMMS_ARCH),serial)

n_body.serial: n_body.o io.o trace.o comms_serial.o pot.o diff.o
	$(SERIAL_F90) -o n_body.serial  n_body.o io.o trace.o comms_serial.o pot.o diff.o $(FCFLAGS)
# Here are the coserialle steps

n_body.o:n_body.f90 io.o trace.o comms_serial.o pot.o diff.o
	$(SERIAL_F90) -c  n_body.f90 $(FCFLAGS)

io.o:./io.f90 trace.o comms_serial.o
	@echo $(OPT)
	$(SERIAL_F90) -c  ./io.f90 -cpp -D life_dir=$(current_dir) -D opt_strat=$(OPT) $(FCFLAGS)

trace.o:./trace.f90
	$(SERIAL_F90) -c  ./trace.f90 $(FCFLAGS)

comms_serial.o:./comms_serial.f90 trace.o
	$(SERIAL_F90) -c  ./comms_serial.f90 $(FCFLAGS)

pot.o:./pot.f90 comms_serial.o io.o trace.o
	$(SERIAL_F90) -c ./pot.f90 $(FCFLAGS)
diff.o:./diff.f90 comms_serial.o io.o trace.o pot.o
	$(SERIAL_F90) -c ./diff.f90 $(FCFLAGS)



#mc_pop.serial: pop.o io.o trace.o comms_serial.o life.o
#	$(SERIAL_F90) -o mc_pop.serial  pop.o io.o trace.o comms_serial.o life.o  $(FCFLAGS)
# Here are the compile steps
#
#pop.o:pop.f90 io.o trace.o comms_serial.o life.o   
#	$(SERIAL_F90) -c  pop.f90 $(FCFLAGS)
#
#io.o:./io.f90 trace.o comms_serial.o
#	$(SERIAL_F90) -c  ./io.f90 -cpp -D life_dir=$(current_dir) -D opt_strat=$(OPT) $(FCFLAGS)
#
#trace.o:./trace.f90
#	$(SERIAL_F90) -c  ./trace.f90 $(FCFLAGS)
#
#comms_serial.o:./comms_serial.f90 trace.o
#	$(SERIAL_F90) -c  ./comms_serial.f90 $(FCFLAGS)
#
#life.o:./life.f90 comms_serial.o io.o trace.o
#	$(SERIAL_F90) -c  ./life.f90 $(FCFLAGS)


endif

