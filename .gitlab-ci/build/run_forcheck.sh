#!/bin/bash

. /etc/profile.d/modules.sh

module load hdf5-mpi intel/16.0 mkl/11.3 impi/5.0.3 forcheck/14.6.26

echo "running forcheck"

make HDF5IO=no forcheck
  
