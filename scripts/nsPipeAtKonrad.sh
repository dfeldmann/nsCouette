#!/bin/bash
#
# Purpose:  load environment for nsPipe on Konrad at HLRN
# Usage:    source nsPipeAtKonrad.sh 
# Author:   Daniel Feldmann
# Date:     07th March 2017
# Modified: 24th January 2019

# https://www.hlrn.de/home/view/System3/CrayBuildingPrograms
# echo "Set development environment for the Konrad HPC system"
# source /etc/profile.d/modules.sh
module avail -P
# module help PrgEnv-intel
module swap PrgEnv-cray PrgEnv-intel
#
# module avail 2>&1 | grep netcdf
# module help cray-netcdf
# module load cray-netcdf
#
# module avail 2>&1 | grep fftw
# module help cray-fftw
# module load cray-fftw 
#
module load gnuplot
module list

echo "Change to nsPipe working directory"
cd $WORK/nsPipe
