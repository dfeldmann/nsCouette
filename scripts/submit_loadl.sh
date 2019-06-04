# @ shell=/bin/bash
# @ job_name = ns_SMALL-0032-10
# @ error   = $(job_name).$(jobid).err
# @ output  = $(job_name).$(jobid).out
# @ job_type = parallel
# @ node_usage= not_shared
# @ node = 16
# @ tasks_per_node = 2
# @ resources = ConsumableCpus(10)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 00:10:00
# @ notification = never
# @ queue

# LoadLeveler submit script for MPG's iDataPlex HPC cluster Hydra @MPCDF
#  (20 cores per node, IBM PE)
# using 16 nodes, 32 MPI tasks, 10 OpenMP threads per task

export OMP_NUM_THREADS=10

export OMP_SCHEDULE=static
export OMP_WAIT_POLICY=active
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_STACKSIZE=256M

#get the number of MPI tasks from hostfile as 16*2
TASKS=`wc -l $LOADL_HOSTFILE| awk '{print $1}'`

echo "TASKS: $TASKS, THREADS: $OMP_NUM_THREADS"

poe ./nsCouette.x < input_nsCouette

