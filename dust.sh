#!/bin/bash

# Copy latest netcdf file here
#cp ../lpi_tests/gp.nc .

# I/O fnames
ifiles="faerie.f90"
ofile="faerie"

# Compiler and flags
fcom=gfortran
fflags="-g -fcheck=all -fbacktrace -Wall -Wextra -O0 -pedantic"
nfflags=`nf-config --fflags`
nfflibs=`nf-config --flibs`
lflibs="-llapack -lblas"

# Compile command
$fcom $fflags $nfflags $ifiles $nfflibs $lflibs -o $ofile

# Optionally run
if [[ $1 == 'run' ]]
then
  ./$ofile
fi
