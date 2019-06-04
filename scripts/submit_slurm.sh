#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -J ns_SMALL-0032-12
#SBATCH -o ns_SMALL-0032-12.%j.out
#SBATCH -e ns_SMALL-0032-12.%j.err
#SBATCH -p test
#SBATCH -N 16

# SLURM submit script for CRAY XC30 machine sisu.csc.fi
#  (24 cores per node, CRAY MPI)
# using 16 nodes, 32 MPI tasks, 12 OpenMP threads per task
# Code was compiled with Intel compilers

module switch PrgEnv-cray PrgEnv-intel

export OMP_NUM_THREADS=12

export OMP_SCHEDULE=static
export OMP_WAIT_POLICY=active
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_STACKSIZE=256M

#compute the number of MPI tasks as 16*24/12
TASKS=$(($SLURM_NNODES*$SLURM_CPUS_ON_NODE/$OMP_NUM_THREADS))
#2 MPI tasks per node
TPN=$(($TASKS/$SLURM_NNODES))

echo "NODES: $SLURM_NNODES, TASKS: $TASKS, TPN: $TPN, THREADS: $OMP_NUM_THREADS"

aprun -n $TASKS -d $OMP_NUM_THREADS -N $TPN -S $(($TPN/2)) -ss  -cc numa_node ./nsCouette.x < input_nsCouette

