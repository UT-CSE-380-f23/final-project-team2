#!/bin/bash

# run our code on ls6!
module purge
module load TACC
module unload XALT
module load boost
module load gsl
module load petsc
module load hdf5
module list 

autoreconf -f -i
echo $TACC_PHDF5_DIR
echo $TACC_HDF5_DIR

# could do a better job of searching for this, but, run the configure file
./configure CC=mpicc CXX="mpicxx -std=c++17 -O3" FC=mpif90 --with-grvy=/work/03453/villa13/ls6/CSE380/PUBLIC/grvy-intel --with-hdf5=$TACC_PHDF5_DIR --with-petsc=$TACC_PETSC_DIR
