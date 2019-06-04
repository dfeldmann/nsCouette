#!/bin/bash

#run the same setup with different number of MPI tasks and demand identity of 
#all output files
#this can/should be used to challenge the uneven distribution of fourier modes
#and/or different parallelizations

set -e
function exittrap {
    echo "$0 failed"
    exit 1
}
trap exittrap ERR INT TERM

ARCH=$1
NUMPROCS1=${2:-4}
NUMTHRDS1=${3:-1}
NUMPROCS2=${4:-8}
NUMTHRDS2=${5:-1}
. /etc/profile.d/modules.sh

module load impi itac intel hdf5-serial

export OMP_PLACES=cores
export OMP_STACKSIZE=128M

DIR1=tmpCI_1
DIR2=tmpCI_2
mkdir $DIR1  
mkdir $DIR2  


cd $DIR1
export OMP_NUM_THREADS=$NUMTHRDS1
mpiexec -n $NUMPROCS1 -check_mpi ../$ARCH/nsCouette.x < ../input_nsCouette.test_short 
cd ..

cd $DIR2
export OMP_NUM_THREADS=$NUMTHRDS2
mpiexec -n $NUMPROCS2 -check_mpi ../$ARCH/nsCouette.x < ../input_nsCouette.test_short 
cd ..

FILES="vel_mid ke_mode temp_mid ke_th ke_z ke_total torque Temp_energy Nusselt restart `ls fields_*.xmf coeff_*|xargs`"
for file in $FILES; do
    if [ -e $DIR1/$file ] && [ -e $DIR2/$file ]; then
        echo "@cmp $DIR1/$file $DIR2/$file"
        cmp $DIR1/$file $DIR2/$file
    fi
done

FILES="`ls fields_*.h5|xargs`"
for file in $FILES; do
    if [ -e $DIR1/$file ] && [ -e $DIR2/$file ]; then
        echo "@h5diff $DIR1/$file $DIR2/$file"
        h5diff $DIR1/$file $DIR2/$file
    fi
done

rm -rf $DIR1 $DIR2


  
