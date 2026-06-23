#!/bin/bash

SAIDA="saida-slurm-exp2-2n-4ppn-2thpp.txt"

echo "rodando no host $(hostname)"
echo "rodando no host $(hostname)" > $SAIDA

echo "SLURM_JOB_NAME: $SLURM_JOB_NAME" | tee -a $SAIDA
echo "SLURM_NODELIST: $SLURM_NODELIST" | tee -a $SAIDA
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST" | tee -a $SAIDA
echo "SLURM_JOB_CPUS_PER_NODE: $SLURM_JOB_CPUS_PER_NODE" | tee -a $SAIDA

echo "==========================================" | tee -a $SAIDA
echo "Experiencia 2" | tee -a $SAIDA
echo "np=4  nth=2" | tee -a $SAIDA
echo "nelements=32000000" | tee -a $SAIDA
echo "npivots=32000" | tee -a $SAIDA
echo "nbins=512" | tee -a $SAIDA
echo "nr=10" | tee -a $SAIDA
echo "==========================================" | tee -a $SAIDA

export OMP_NUM_THREADS=1

mpirun -np 4 ./MPI+Pthreads-Multipartition \
    32000000 \
    32000 \
    512 \
    2 \
    10 \
    | tee -a $SAIDA

echo "Tempo total da shell: $SECONDS segundos" \
    | tee -a $SAIDA

squeue -j $SLURM_JOBID \
    | tee -a $SAIDA
