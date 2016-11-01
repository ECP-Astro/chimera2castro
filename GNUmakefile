CASTRO_HOME ?= ../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

#USE_ASSERTION = FALSE
#TEST          = TRUE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

FFLAGS    += -ffree-line-length-none
fFLAGS    += -ffree-line-length-none

DEFINES += -DDO_PROBLEM_POST_INIT

USE_MPI    = TRUE
USE_OMP    = TRUE

USE_GRAV   = TRUE
USE_GR     = FALSE
USE_POINTMASS = TRUE

USE_REACT  = TRUE

USE_RATES     = TRUE
USE_SCREENING = TRUE
USE_NEUTRINOS = TRUE

USE_MODELPARSER  = TRUE

USE_HDF5   = TRUE

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_dir     := helmholtz
#EOS_dir     := stellarcollapse

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
#Network_dir := aprox21
#Network_dir := anp56
Network_dir := general_null
GENERAL_NET_INPUTS = $(MICROPHYSICS_HOME)/networks/$(Network_dir)/anp56.net

INTEGRATOR_DIR := BS

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
