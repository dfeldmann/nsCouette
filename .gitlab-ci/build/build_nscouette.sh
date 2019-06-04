#!/bin/bash

set -e
function exittrap {
    echo "$0 failed"
    exit 1
}
trap exittrap ERR INT TERM


CODE_VARIANT=$1
ARCHITECTURE=$2
DEBUG=${3:-no}

. /etc/profile.d/modules.sh

module load git

 
case "$ARCHITECTURE" in
intel-mkl)
  module load hdf5-mpi intel/16.0 mkl/11.3 impi/5.0.3
  ;;
gcc-mkl)
  module load gcc/4.9 mkl impi/5.0.3  hdf5-mpi-gcc
  ;;
gcc-openblas)
  module load gcc/4.9 fftw/gcc/3.3.4 impi/5.1.1
  export HDF5_HOME=/afs/@cell/common/soft/hdf5/1.8.15-patch1/amd64_sles11/gcc/4.9/impi-5.1.1/
  export OPENBLAS_HOME=/afs/rzg/.cs/openblas/0.2.13/gcc-4.8
  ;;
pgi-mkl)
  module load pgi mkl impi/5.0.3  hdf5-mpi-gcc
  ;;
*)
  echo "ARCHITECTURE $ARCHITECTURE not defined for build"
  exit 1
  ;;
esac

module list

echo "building nscouette with ARCH=$ARCHITECTURE, CODE=$CODE_VARIANT"

make ARCH=$ARCHITECTURE HDF5IO=yes PROFLIB=FTIM CODE=$CODE_VARIANT DEBUG=$DEBUG


  
