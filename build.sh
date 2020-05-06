#!/usr/bin/env bash

TOPDIR=$PWD

# Build and install the Fortran Fit Codes
cd $TOPDIR/dsarrfit
bash build.sh
# cd $TOPDIR/ # Need to fix TroeFit
# bash build.sh

# Install Python Interfaces to RateFits and Other Functionality
cd $TOPDIR
$PYTHON setup.py install
