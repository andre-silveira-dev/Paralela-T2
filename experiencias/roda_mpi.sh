#!/bin/bash

# compile copia.c com: mpicc -O3 -g -Wall -Wextra -o teste_mpi copia.c -lpthread
# rode com: sbatch --exclusive -N 1 ./roda_mpi.sh
# veja a saida gerada no slurm

mpirun -np 4 ./teste_mpi 32 16 4 1 1
