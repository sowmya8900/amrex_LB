# # AMREX_HOME defines the directory in which we will find all the AMReX code.
# # If you set AMREX_HOME as an environment variable, this line will be ignored
# AMREX_HOME ?= ../../../amrex

# DEBUG        = FALSE
# USE_CCACHE   = TRUE
# USE_OMP      = FALSE
# TINY_PROFILE = TRUE
# USE_MPI      = FALSE 
# USE_CUDA     = FALSE
# BL_NO_FORT   = TRUE
# COMP         = gnu
# DIM          = 3

# # Add additional include paths if needed
# INCLUDE_LOCATIONS += ../../../amrex-install/include
# INCLUDE_LOCATIONS += .

# include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# # include ./Make.package
# include Make.package
# include $(AMREX_HOME)/Src/Base/Make.package
# include $(AMREX_HOME)/Tools/GNUMake/Make.rules


# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex
# AMREX_HOME ?= ../../amrex

# AMREX_HOME ?= ../amrex

#DEBUG        = TRUE
#USE_CCACHE   = TRUE
USE_OMP      = FALSE
TINY_PROFILE = TRUE
USE_MPI      = FALSE 
USE_CUDA     = FALSE
BL_NO_FORT   = TRUE
COMP         = gnu
DIM          = 3
ARCH         = x86-milan

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# include ./Make.package
include Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.rules