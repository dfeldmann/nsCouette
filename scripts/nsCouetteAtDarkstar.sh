#!/bin/bash
#
# Purpose:  load environment for nsCouette on my local desktop machine
# Usage:    source nsCouetteAtDarkstar.sh 
# Author:   Daniel Feldmann
# Date:     07th March 2017
# Modified: 09th October 2018

#--- local libs
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/feldmann/fftw/fftw-3.3.6-pl2/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/feldmann/hdf5/hdf5-1.10.0-patch1/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/feldmann/zlib/zlib-1.2.11/lib
export PATH=$PATH:/home/feldmann/fftw/fftw-3.3.6-pl2/bin
export PATH=$PATH:/home/feldmann/hdf5/hdf5-1.10.0-patch1/bin

#--- visualisation
export PATH=$PATH:$HOME/ParaView/ParaView-5.5.2-Qt5-MPI-Linux-64bit/bin/
export PATH=$PATH:$HOME/VisIt/visit2_7_3.linux-x86_64/bin/

cd $HOME/nsCouette
echo "Changed to nsCouette case directory:"
pwd
