#########################################
#   TOP LEVEL MAKEFILE                  #
#########################################

SOURCE = ./Source
BUILD_DIR = ./Build

#########################################
# User editable options

# Communications architectiure: mpi,serial
COMMS_ARCH:=mpi

# Fortran compiler, allowed: ifort, gfortran
F90:=gfortran

# Speed of the build, allowed: fast, debug
BUILD= fast


########################################


MPI_F90:=mpif90
SERIAL_F90=$(F90)
export F90
export COMMS_ARCH
export MPI_F90
export SERIAL_F90
export BUILD

ifeq ($(COMMS_ARCH),mpi)





subsystem:
	export $(COMMS_ARCH)
	$(MAKE) -C $(SOURCE)
	mkdir  -p Build/$(COMMS_ARCH)_$(F90)
	install -m 557 $(SOURCE)/n_body.mpi $(BUILD_DIR)/$(COMMS_ARCH)_$(F90)
	mv $(SOURCE)/*.o $(BUILD_DIR)/$(COMMS_ARCH)_$(F90)
	mv $(SOURCE)/*.mod $(BUILD_DIR)/$(COMMS_ARCH)_$(F90)
	rm -f $(SOURCE)/n_body.mpi	

.phony: install

clean:
	rm -f $(BUILD_DIR)/$(COMMS_ARCH)_$(F90)/*



endif 


ifeq ($(COMMS_ARCH),serial)


subsystem:
	export $(COMMS_ARCH)
	$(MAKE) -C $(SOURCE)
	mkdir  -p Build/$(COMMS_ARCH)_$(F90)
	install -m 557 $(SOURCE)/n_body.serial $(BUILD_DIR)/$(COMMS_ARCH)_$(F90)/n_body.serial
	mv $(SOURCE)/*.o $(BUILD_DIR)/$(COMMS_ARCH)_$(F90)
	mv $(SOURCE)/*.mod $(BUILD_DIR)/$(COMMS_ARCH)_$(F90)


	rm -f $(SOURCE)/n_body.serial
.phony: install



clean:
	rm -f $(BUILD_DIR)/$(COMMS_ARCH)_$(F90)/*


endif
clean_all:
	rm -f -r $(BUILD_DIR)/*/ 

dist:
	tar  --exclude="./.git" --exclude="./Test" --exclude="./*/*.mpi" --exclude="./*/*.serial" --exclude="./Source/*.o" --exclude="./Source/*.mod" -cvf N_BODY.tar .

