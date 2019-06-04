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
REFCASE=${4:reference_run_1}
. /etc/profile.d/modules.sh

#run the code
module load impi

export OMP_NUM_THREADS=$NUMTHRDS
export OMP_PLACES=cores
export OMP_STACKSIZE=128M

mpiexec -n $NUMPROCS $ARCH/nsCouette.x < .gitlab-ci/test/$REFCASE/input_nsCouette > nscouette.out

cat /proc/cpuinfo > environment 
env >> environment

#compare results
NUMDIFF=/home/rzg/soft/numdiff/bin/numdiff

RES=.
REF=.gitlab-ci/test/$REFCASE

#informational output
FILES="ke_mode  ke_th  ke_total  ke_z  torque  Nusselt"
for file in $FILES; do
    if [ -e $REF/$file ]; then
        echo "@comparing $RES/$file $REF/$file" 
        $NUMDIFF --quiet --statistics $RES/$file $REF/$file || echo "info output" 
    fi
done

#verification
FILES="torque  Nusselt"
for file in $FILES; do
    if [ -e $REF/$file ]; then
        echo "@verifying $RES/$file $REF/$file" 
        $NUMDIFF --relative=2e-8 $RES/$file $REF/$file
    fi
done
