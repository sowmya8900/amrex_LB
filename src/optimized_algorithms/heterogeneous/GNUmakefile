AMREX_HOME ?= ../../../../amrex

DEBUG	= FALSE
DIM	= 3
COMP    = gnu
PRECISION = DOUBLE

USE_MPI   = FALSE
USE_OMP   = TRUE
USE_CUDA  = FALSE
USE_HIP   = FALSE

BL_NO_FORT = TRUE

EBASE = main

# GPU settings
CUDA_ARCH = 80
CUDA_MAX_THREADS = 256

# Add source directories
VPATH_LOCATIONS += ..
INCLUDE_LOCATIONS += ..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

# Add our source files
CEXE_sources += main.cpp
CEXE_sources += example.cpp
CEXE_sources += HeterogeneousLB.cpp
CEXE_sources += test_heterogeneous.cpp
CEXE_sources += Knapsack.cpp
CEXE_sources += SFC.cpp
CEXE_sources += SFC_knapsack.cpp
CEXE_sources += LeastUsed.cpp

include $(AMREX_HOME)/Tools/GNUMake/Make.rules 