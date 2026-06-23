#!/bin/bash

echo "rodando no host $(hostname)"
echo "rodando no host $(hostname)" > saida-slurm-exp1-1n-8ppn-1thpp.txt

echo "SLURM_JOB_NAME: $SLURM_JOB_NAME" | tee -a saida-slurm-exp1-1thpp.txt
echo "SLURM_NODELIST: $SLURM_NODELIST" | tee -a saida-slurm-exp1-1thpp.txt
echo "SLURM_JOB_CPUS_PER_NODE: $SLURM_JOB_CPUS_PER_NODE" | tee -a saida-slurm-exp1-1thpp.txt

echo "==========================================" | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt
echo "Experiencia 1" | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt
echo "np=8  nth=1" | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt
echo "nelements=32000000" | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt
echo "npivots=32000" | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt
echo "nbins=512" | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt
echo "nr=10" | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt
echo "==========================================" | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt

export OMP_NUM_THREADS=1

mpirun -np 8 ./MPI+Pthreads-Multipartition \
    32000000 \
    32000 \
    512 \
    1 \
    10 \
    | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt

echo "Tempo total da shell: $SECONDS segundos" \
    | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt

squeue -j $SLURM_JOBID \
    | tee -a saida-slurm-exp1-1n-8ppn-1thpp.txt
