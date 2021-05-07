#!/usr/bin/env bash

TOPDIR=$PWD

# Build and install the Fortran Fit Codes
cd $TOPDIR/external/dsarrfit
bash build.sh
cd $TOPDIR/external/troefit
bash build.sh

# Install Python Interfaces to RateFits and Other Functionality
cd $TOPDIR
python setup.py install
