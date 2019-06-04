#!/bin/bash
#
# Purpose:  load environment for nsCouette on fsmcluster at ZARM
# Usage:    source nsCouetteAtFsm.sh 
# Author:   Daniel Feldmann
# Date:     07th March 2017
# Modified: 09th October 2018

#--- Intel compiler suit and other stuff
module avail
module purge
module load Paraview/5.5.2
module load Intel/PSXE2018
module load Python/3.5.1
module load Anaconda3/5.1.0
module list
which python

#--- local libs
export PATH=$PATH:/home/feldmann/hdf5/hdf5-1.10.0-patch1/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/feldmann/hdf5/hdf5-1.10.0-patch1/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/feldmann/zlib/zlib-1.2.11/lib
# echo $PATH $LD_LIBRARY_PATH | egrep 'curl|hdf5|netcdf|zlib|fftw'

cd $HOME/nsCouette
echo "Changed to nsCouette case directory:"
pwd
