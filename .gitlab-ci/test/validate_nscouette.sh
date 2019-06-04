#!/bin/bash

# compare the numerically computed wave speed with an experimentally
# determined value 0.34432 (King et al., J Fluid Mech, 1984)

set -e 
function exittrap {     
    echo "$0 failed"
    exit 1 
} 
trap exittrap ERR INT TERM


ARCHITECTURE=$1
NUMPROCS=${2:-4}
NUMTHRDS=${3:-2}
. /etc/profile.d/modules.sh

#code duplicated from build_nscouette.sh
case "$ARCHITECTURE" in
intel-mkl)
  module load hdf5-mpi intel/16.0 mkl/11.3 impi
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


#build the code and postprocessing program
make HDF5IO=no CODE=STD_CODE PROFLIB=FTIM ARCH=$ARCHITECTURE
make HDF5IO=no CODE=STD_CODE PROFLIB=FTIM ARCH=$ARCHITECTURE postproc

#run the code
rm -f coeff_*.*

cat /proc/cpuinfo > environment 
env >> environment

ln -fs .gitlab-ci/test/validation_run_1/coeff_WSPEED.00000000
cp .gitlab-ci/test/validation_run_1/coeff_WSPEED.00000000.info restart

export OMP_NUM_THREADS=$NUMTHRDS
export OMP_PLACES=cores
export OMP_STACKSIZE=128M

mpiexec -n $NUMPROCS $ARCHITECTURE/nsCouette.x < .gitlab-ci/test/validation_run_1/input_nsCouette > nscouette.out


#compute wave speed
$ARCHITECTURE/waveSpeed.x < .gitlab-ci/test/validation_run_1/input_waveSpeed > wavespeed.out

#validate results
wspeed=`grep '# non-dimensional wavespeed' wavespeed.out | awk -F'=' '{print $2}'`

wspeed_max=0.35
wspeed_min=0.33
echo "comparing computed wavespeed with experimentally determined value 0.34432"
if [ $wspeed \> $wspeed_min ] && [ $wspeed \< $wspeed_max ]; then
    echo "SUCCESS: wavespeed=$wspeed within bounds [$wspeed_min,$wspeed_max]"
    exit 0
else
    echo "FAILURE: wavespeed=$wspeed out of bounds [$wspeed_min,$wspeed_max]"
    exit 1
fi

