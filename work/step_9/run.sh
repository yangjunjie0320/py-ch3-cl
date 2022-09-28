#!/bin/bash                                                                          
#SBATCH --partition=debug
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --job-name=TDDFT
#SBATCH --mem=0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;
export PYSCF_MAX_MEMORY=$SLURM_MEM_PER_NODE;
echo SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK;
echo SLURM_MEM_PER_NODE  = $SLURM_MEM_PER_NODE;

export PYSCF_TMPDIR=/scratch/global/yangjunjie/;
export PYTHONPATH=/home/yangjunjie/packages/pyscf/pyscf-main

python main.py py-ch3-cl.xyz
