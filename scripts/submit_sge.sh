#$ -S /bin/bash
#$ -l h_rt=00:10:00
#$ -N ns_SMALL-0032-12
#$ -j n
#$ -cwd
#$ -m n
#$ -pe impi_hydra 256 

# SGE submit script for SandyBridge Linux cluster @RZG
#  (16 cores per node, Intel MPI)
# using 16 nodes, 32 MPI tasks, 8 OpenMP threads per task

CORES_PER_NODE=16

export OMP_NUM_THREADS=8

export OMP_SCHEDULE=static
export OMP_WAIT_POLICY=active
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_STACKSIZE=256M

export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_CELL=core

#compute the number of MPI tasks as 256/8
TASKS=$(($NSLOTS/$OMP_NUM_THREADS))
#2 MPI tasks per node
TPN=$(($CORES_PER_NODE/$OMP_NUM_THREADS))

echo "NSLOTS: $NSLOTS, TPN: $TPN, TASKS: $TASKS, THREADS: $OMP_NUM_THREADS"

#with Intel MPI it might be worth trying a different alltoall algorithm
#export I_MPI_ADJUST_ALLTOALL=3

mpiexec -n $TASKS -perhost $TPN ./nsCouette.x < input_nsCouette

