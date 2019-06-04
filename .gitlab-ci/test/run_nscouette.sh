#!/bin/bash

set -e
function exittrap {
    echo "$0 failed"
    exit 1
}
trap exittrap ERR INT TERM

ARCH=$1
NUMPROCS=${2:-4}
NUMTHRDS=${3:-2}
. /etc/profile.d/modules.sh

module load impi itac

export OMP_NUM_THREADS=$NUMTHRDS
export OMP_PLACES=cores
export OMP_STACKSIZE=128M

mpiexec -n $NUMPROCS -check_mpi $ARCH/nsCouette.x < input_nsCouette.test_short 


  
