
PRECISION  = DOUBLE

DEBUG      = FALSE
#TINY_PROFILE = TRUE
#PROFILE    = TRUE

#USE_ASSERTION = FALSE
#TEST          = TRUE

DIM        = 2

COMP	   = gnu

#FFLAGS    += -ffree-line-length-none
#fFLAGS    += -ffree-line-length-none

DEFINES += -DDO_PROBLEM_POST_INIT

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_GRAV   = TRUE
USE_GR     = FALSE
USE_POINTMASS = TRUE

USE_REACT  = FALSE

USE_RATES     = FALSE
USE_SCREENING = FALSE
USE_NEUTRINOS = FALSE

USE_MODELPARSER  = TRUE

USE_HDF5   = TRUE

USE_FORTRAN_INTERFACE = TRUE
USE_F_INTERFACES = TRUE

libraries += -lhdf5 -lhdf5_fortran -lz
LDFLAGS += -L$(HDF5_DIR)/lib
VPATH_LOCATIONS += $(HDF5_DIR)/include
FINCLUDE_LOCATIONS += $(HDF5_DIR)/include


ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz
#EOS_DIR     := stellarcollapse

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
#Network_dir := aprox21
#Network_dir := anp56
#Network_dir := general_null
NETWORK_DIR := aprox13
GENERAL_NET_INPUTS = $(MICROPHYSICS_HOME)/networks/$(NETWORK_DIR)/anp56.net

INTEGRATOR_DIR := VODE

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
