parallel-multiPartition: main.c
	mpicc -O3 -g -Wall -Wextra -o parallel_multiPartition main.c -lpthread 
