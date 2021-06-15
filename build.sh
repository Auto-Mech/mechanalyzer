#!/usr/bin/env bash

export TOPDIR=$PWD
export FC=$(which gfortran)

# Build and install the Fortran Fit Codes
cd $TOPDIR/external/dsarrfit
bash build.sh
cd $TOPDIR/external/troefit
bash build.sh

# Install Python Interfaces to RateFits and Other Functionality
cd $TOPDIR
python setup.py install
