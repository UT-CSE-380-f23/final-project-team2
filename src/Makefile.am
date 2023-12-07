AM_CXXFLAGS   = -O2
AM_CXXFLAGS  += $(GRVY_CFLAGS) $(GSL_CFLAGS) $(HDF5_CFLAGS)
if PETSC_ENABLED
    include $(PETSC_DIR)/lib/petsc/conf/variables
    AM_CXXFLAGS += -DINCLUDE_PETSC $(PETSC_CC_INCLUDES)
endif
bin_PROGRAMS  = run
run_SOURCES = main.cpp fd_solver.cpp oned_solver.cpp parsing.cpp twod_solver.cpp
run_LDADD   = $(GRVY_LIBS) $(GSL_LIBS) $(HDF5_LIBS)
if PETSC_ENABLED
    run_LDADD += ${PETSC_LIB}
endif