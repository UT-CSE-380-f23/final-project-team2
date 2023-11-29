#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([cse380_jen_ken_rod], [0.0.1], [kmeyer2@utexas.edu])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADER(src/config.h)
AM_INIT_AUTOMAKE
# Checks for programs.
AC_PROG_CXX
# Checks for libraries.

#AC_CHECK_LIB([m],[cos])
#AC_CHECK_LIB([gslcblas],[cblas_dgemm])
#AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AX_PATH_GRVY([0.32],[yes])
AX_PATH_GSL(2.7,[],AC_MSG_ERROR([Unable to locate GSL or meet minimum version requirement]))
AX_PATH_HDF5([1.14.0],[yes])
AX_PATH_PETSC([3.19],[no])
# Checks for code coverage.
AX_CODE_COVERAGE
# Checks for header files.
# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
# Checks for library functions.
#AC_OUTPUT([Makefile src/Makefile test/Makefile])
AC_OUTPUT([Makefile src/Makefile])