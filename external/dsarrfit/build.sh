#!/usr/bin/env bash

TOPDIR=$PWD

mkdir -p build
cd build

cmake -DCMAKE_INSTALL_PREFIX=$TOPDIR -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_Fortran_FLAGS="${FFLAGS}" ..

make VERBOSE=1
make install
