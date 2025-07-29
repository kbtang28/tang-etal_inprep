#!/bin/bash
#SBATCH --job-name=base_scenarios_independence
#SBATCH --nodes=1
#SBATCH --ntasks=34
#SBATCH --cpus-per-task=2
#SBATCH --array=1-3
#SBATCH --output=slurm/out-%A_%a.txt
#SBATCH --error=slurm/err-%A_%a.txt
#SBATCH --mail-user=kbt28@cornell.edu
#SBATCH --mail-type=END

# print properties of job as submitted
echo "SLURM_JOB_ID = $SLURM_JOB_ID"
echo "SLURM_ARRAY_JOB_ID = $SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID"
echo "SLURM_JOB_NAME = $SLURM_JOB_NAME"
echo "SLURM_NTASKS = $SLURM_NTASKS"
echo "SLURM_NTASKS_PER_NODE = $SLURM_NTASKS_PER_NODE"
echo "SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"
echo "SLURM_JOB_NUM_NODES = $SLURM_JOB_NUM_NODES"

# print properties of job as scheduled by Slurm
echo "SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo "SLURM_TASKS_PER_NODE = $SLURM_TASKS_PER_NODE"
echo "SLURM_JOB_CPUS_PER_NODE = $SLURM_JOB_CPUS_PER_NODE"
echo "SLURM_CPUS_ON_NODE = $SLURM_CPUS_ON_NODE"

# run script
julia +1.10.0 workflow/base_scenarios/test_independence.jl

exit 0