/** 
 * Trabalho 2 - Programação Paralela
 * Multipartição paralela usando MPI
 * 
 * Discentes:
 *   Andreus Gustavo Schultz - GRR20232374
 *   André Gustavo Silvera - GRR20232345
*/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#include "chrono.c"

#define MAX_THREADS 8
#define CACHE_SIZE_MB 3
#define EVICTION_MULTIPLIER 3

#include <mpi.h>

// variaveis de entrada
long long int nelements = 0, 
    npivots = 0, 
    nbins = 0, 
    nthreads = 0, 
    nr = 0;
static int use_verify = 0;
static int use_tb2 = 0;

// flag para evitar a criação de threads
static int Histo_thread_pool_initialized = 0;
static int MP_thread_pool_initialized = 0;

// variáveis de gerenciamento de pool
pthread_barrier_t parallelHisto_barrier;
pthread_barrier_t parallelMP_barrier;
int parallelHisto_threads_id[MAX_THREADS];
int parallelMP_threads_id[MAX_THREADS];
pthread_t parallelHisto_thread[MAX_THREADS];
pthread_t parallelMP_thread[MAX_THREADS];
static int parallelHisto_pool_threads = 0;
static int parallelMP_pool_threads = 0;

// variáveis de dados
static const long long *parallelHisto_data = NULL;
static const long long *parallelHisto_limits = NULL;
static long long **parallelHisto_local_hist = NULL;
static long long **parallelMP_local_output = NULL;
static long long *parallelMP_global_pos = NULL;
static long long parallelHisto_nelements = 0;
static int parallelHisto_nbins = 0;
static int parallelHisto_active_threads = 1;
static int local_partition_verify_ok = 1;

static long long bin_size(
    int bin,
    const long long *Pos,
    long long nElements,
    int nbins)
{
    if (bin == nbins - 1)
        return nElements - Pos[bin];

    return Pos[bin + 1] - Pos[bin];
}

long long *multi_partition_mpi(
    const long long *Output,
    const long long *Pos,
    long long nElements,
    int nbins,
    int *nElementsReceived,
    MPI_Comm comm
) {
    int rank;
    int nproc;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    int bins_per_rank = nbins / nproc;

    int *sendcounts = calloc(nproc, sizeof(int));

    for (int bin = 0; bin < nbins; bin++) {
        int dest = bin / bins_per_rank;

        long long size = bin_size(bin, Pos, nElements, nbins);

        sendcounts[dest] += (int)size;
    }

    int *sdispls = calloc(nproc, sizeof(int));

    for (int i = 1; i < nproc; i++) {
        sdispls[i] = sdispls[i - 1] + sendcounts[i - 1];
    }

    long long *sendbuf = malloc(nElements * sizeof(long long));

    int *cursor = calloc(nproc, sizeof(int));

    for (int i = 0; i < nproc; i++) {
        cursor[i] = sdispls[i];
    }
    
    for (int bin = 0; bin < nbins; bin++) {
        int dest = bin / bins_per_rank;

        long long begin = Pos[bin];

        long long size = bin_size(bin, Pos, nElements, nbins);

        memcpy(
            &sendbuf[cursor[dest]],
            &Output[begin],
            size * sizeof(long long)
        );

        cursor[dest] += size;
    }

    free(cursor);

    int *recvcounts = calloc(nproc, sizeof(int));

    MPI_Alltoall(
        sendcounts,
        1,
        MPI_INT,
        recvcounts,
        1,
        MPI_INT,
        comm
    );

    int *recvdispls = calloc(nproc, sizeof(int));

    for (int i = 1; i < nproc; i++) {
        recvdispls[i] = recvdispls[i - 1] + recvcounts[i - 1];
    }

    int total_recv = 0;

    for (int i = 0; i < nproc; i++)
        total_recv += recvcounts[i];

    long long *recvbuf = malloc((size_t)total_recv * sizeof(long long));

    MPI_Alltoallv(
        sendbuf,
        sendcounts,
        sdispls,
        MPI_LONG_LONG,
        recvbuf,
        recvcounts,
        recvdispls,
        MPI_LONG_LONG,
        comm
    );

    free(sendbuf);
    free(sendcounts);
    free(sdispls);
    free(recvcounts);
    free(recvdispls);

    *nElementsReceived = total_recv;

    return recvbuf;
}

static inline unsigned long long rand63(void) {
    return (  (unsigned long long)(unsigned)rand()        )
           | ( (unsigned long long)(unsigned)rand() << 21 )
           | ( (unsigned long long)(unsigned)rand() << 42 );
}

static inline void ll_swap(long long *a, long long *b) {
    long long temp = *a;
    *a = *b;
    *b = temp;
}

static inline long long rand64(void) {
    union {
        unsigned long long u;
        long long          s;
    } x;

    x.u = ((unsigned long long)(unsigned)rand() << 49)
        | ((unsigned long long)(unsigned)rand() << 34)
        | ((unsigned long long)(unsigned)rand() << 19)
        | ((unsigned long long)(unsigned)rand() <<  4)
        | ((unsigned long long)(unsigned)rand() & 0xF);

    return x.s;
}


static void gen_test_data_balanced2(
    long long *data,
    long long nelements,
    int nbins
) {
    if (nelements <= 0 || nbins <= 0) return;

    for (long long i = 0; i < nelements; i++)
        data[i] = i;

    for (long long i = nelements - 1; i > 0; i--) {
        long long j = (long long)(rand63() % (unsigned long long)(i + 1));
        ll_swap(&data[i], &data[j]);
    }

    for (long long i = 0; i < nelements; i++)
        data[i] %= nbins;
}

static volatile char *eviction_buffer = NULL;
static size_t eviction_size = 0;

static void evictCacheInit(void) {
    size_t cache_size_bytes = (size_t)CACHE_SIZE_MB * 1024 * 1024;
    eviction_size = cache_size_bytes * EVICTION_MULTIPLIER;
    eviction_buffer = (volatile char *)malloc(eviction_size);
    if (!eviction_buffer) {
        fprintf(stderr, "Error: Could not allocate eviction buffer\n");
        exit(1);
    }
    memset((void *)eviction_buffer, 0, eviction_size);
}

static void evictCache(void) {
    volatile char result = 0;
    for (size_t i = 0; i < eviction_size; i += 64) {
        result += eviction_buffer[i];
    }
    (void)result;
}

static int cmp_ll(const void *a, const void *b) {
    long long va = *(const long long *)a;
    long long vb = *(const long long *)b;
    return (va > vb) - (va < vb);
}

static void generate_local_pivots(
    const long long *data,
    long long n_local,
    int npivots_local,
    long long *pivots_local)
{
    unsigned long long stride = n_local / npivots_local;

    for(int i = 0; i < npivots_local; i++)
    {
        long long jitter_i =
            (long long)(rand63() % stride);

        pivots_local[i] =
            data[i * stride + jitter_i];
    }
}

void build_limits_mpi(
    const long long *data,
    long long n_local,
    int npivots,
    int nbins,
    long long *limits,
    MPI_Comm comm
) {
    int rank;
    int nproc;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    if(npivots % nproc != 0) {
        if(rank == 0) {
            fprintf(
                stderr,
                "Erro: npivots (%d) deve ser multiplo de nproc (%d)\n",
                npivots,
                nproc
            );
        }

        MPI_Abort(comm, 1);
    }

    int npivots_local = npivots / nproc;

    long long *pivots_local = malloc(npivots_local * sizeof(long long));

    generate_local_pivots(
        data,
        n_local,
        npivots_local,
        pivots_local
    );

    long long *all_pivots = NULL;

    if(rank == 0) {
        all_pivots = malloc(npivots * sizeof(long long));
    }

    MPI_Gather(
        pivots_local,
        npivots_local,
        MPI_LONG_LONG,
        all_pivots,
        npivots_local,
        MPI_LONG_LONG,
        0,
        comm
    );

    if(rank == 0) {
        qsort(
            all_pivots,
            npivots,
            sizeof(long long),
            cmp_ll
        );

        limits[0] = LLONG_MIN;
        limits[nbins] = LLONG_MAX;

        int step = npivots / (nbins - 1);

        for(int k = 1; k < nbins; k++) {
            limits[k] = all_pivots[k * step];
        }

        for(int k = 1; k <= nbins - 1; k++) {
            if(limits[k] <= limits[k - 1]) {
                if(limits[k - 1] < LLONG_MAX - 1)
                    limits[k] = limits[k - 1] + 1;
                else
                    limits[k] = limits[k - 1];
            }
        }

        free(all_pivots);
    }

    MPI_Bcast(
        limits,
        nbins + 1,
        MPI_LONG_LONG,
        0,
        comm
    );

    free(pivots_local);
}


static inline int findBin(const long long *limits, int nbins, long long value) {
    int lo = 0;
    int hi = nbins;

    while (lo + 1 < hi) {
        int mid = lo + (hi - lo) / 2;
        if (value < limits[mid]) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    if (lo < 0) return 0;
    if (lo >= nbins) return nbins - 1;
    return lo;
}

void *histogram(void *ptr){
    int tid = *(int *)ptr;

    while(1) {
        pthread_barrier_wait(&parallelHisto_barrier);

        long long start = ((long long)tid * parallelHisto_nelements) / parallelHisto_active_threads;
        long long end = ((long long)(tid + 1) * parallelHisto_nelements) / parallelHisto_active_threads;
        for (long long i = start; i < end; i++) {
            int bin = findBin(parallelHisto_limits, parallelHisto_nbins, parallelHisto_data[i]);
            parallelHisto_local_hist[tid][bin]++;
        }

        pthread_barrier_wait(&parallelHisto_barrier);
        
        if(tid == 0) return NULL;
    }

    return NULL;
}

void *build_output(void *ptr){
    int tid = *(int *)ptr;
    
    long long *next_local = (long long*)malloc(sizeof(long long) * nbins);
    memset(next_local, 0, sizeof(long long) * nbins);

    while(1) {
        pthread_barrier_wait(&parallelMP_barrier);

        long long start = ((long long)tid * parallelHisto_nelements) / parallelHisto_active_threads;
        long long end = ((long long)(tid + 1) * parallelHisto_nelements) / parallelHisto_active_threads;
        for (long long i = start; i < end; i++) {
            int bin = findBin(parallelHisto_limits, parallelHisto_nbins, parallelHisto_data[i]);
            parallelMP_local_output[tid][parallelMP_global_pos[bin] + next_local[bin]] = parallelHisto_data[i];
            next_local[bin]++;
        }

        pthread_barrier_wait(&parallelMP_barrier);
        
        if(tid == 0){
            free(next_local);
            return NULL;
        }
    }

    free(next_local);
    return NULL;
}

int parallelHistogram(
    const long long  *data,
    long long         nelements,
    const long long  *limits,
    int               nbins,
    long long        *hist,
    int               nthreads
){
    if (!data || !limits || !hist || nbins <= 0 || nelements < 0) {
        return 1;
    }

    if (nthreads <= 1) {
        for (long long i = 0; i < nelements; i++) {
            int bin = findBin(limits, nbins, data[i]);
            hist[bin]++;
        }

        return 0;
    }

    if (nthreads > MAX_THREADS) {
        return 1;
    }

    if (!Histo_thread_pool_initialized) {
        pthread_barrier_init(&parallelHisto_barrier, NULL, nthreads);

        parallelHisto_pool_threads = nthreads;
        parallelHisto_threads_id[0] = 0;
        for (int i = 1; i < nthreads; i++) {
            parallelHisto_threads_id[i] = i;
            pthread_create(&parallelHisto_thread[i], NULL, histogram, &parallelHisto_threads_id[i]);
        }

        Histo_thread_pool_initialized = 1;
    } else if (nthreads != parallelHisto_pool_threads) {
        return 1;
    }

    parallelHisto_local_hist = (long long **)malloc((size_t)nthreads * sizeof(long long *));
    if (!parallelHisto_local_hist) {
        return 1;
    }

    for (int t = 0; t < nthreads; t++) {
        parallelHisto_local_hist[t] = (long long *)calloc((size_t)nbins, sizeof(long long));
        if (!parallelHisto_local_hist[t]) {
            for (int j = 0; j < t; j++) {
                free(parallelHisto_local_hist[j]);
            }
            free(parallelHisto_local_hist);
            parallelHisto_local_hist = NULL;
            return 1;
        }
    }

    parallelHisto_data = data;
    parallelHisto_limits = limits;
    parallelHisto_nelements = nelements;
    parallelHisto_nbins = nbins;
    parallelHisto_active_threads = nthreads;

    histogram(&parallelHisto_threads_id[0]);

    for (int b = 0; b < nbins; b++) {
        long long sum = 0;
        for (int t = 0; t < nthreads; t++) {
            sum += parallelHisto_local_hist[t][b];
        }
        hist[b] = sum;
    }

    return 0;
}

void prefix_sum(long long * Pos, long long * hist, int nbins){
    Pos[0] = 0;
    for(int i = 1; i < nbins; i++)
        Pos[i] = Pos[i-1] + hist[i-1];
}


int parallel_multiPartition(
    const long long  *Input,       /* input array                          */
    long long        *Output,      /* Output array (partitioned per bin)   */
    long long         nElements,   /* number of elements                   */
    const long long  *Limits,      /* bin boundaries, size nbins+1         */
    int               nbins,       /* number of bins                       */
    long long        *Pos,         /* output histogram, size nbins         */
    int               nthreads     /* number of threads (1 = serial path)  */
){
    long long * local_histogram = (long long *)malloc(nbins * sizeof(long long));
    memset(local_histogram, 0, sizeof(long long) * nbins);

    if(parallelHistogram(Input, nElements, Limits, nbins, local_histogram, nthreads)){
        free(local_histogram);
        return 1;
    }

    if(!Input || !Output || !Limits || !Pos || !local_histogram){
        free(local_histogram);
        return 1;
    }

    if(nthreads > MAX_THREADS){
        free(local_histogram);
        return 1;
    }
    
    prefix_sum(Pos, local_histogram, nbins);

    if(nthreads <= 1){
        long long *next = (long long*)malloc(sizeof(long long) * nbins);
        if(!next) return 1;
        memset(next, 0, sizeof(long long) * nbins);

        memcpy(next, Pos, sizeof(long long) * nbins);

        for(int i = 0; i < nelements; i++){
            int bin = findBin(Limits, nbins, Input[i]);
            Output[next[bin]] = Input[i];
            next[bin]++;
        }

        free(local_histogram);
        free(next);
        return 0;
    }

    if (!MP_thread_pool_initialized) {
        pthread_barrier_init(&parallelMP_barrier, NULL, nthreads);

        parallelMP_pool_threads = nthreads;
        parallelMP_threads_id[0] = 0;
        for (int i = 1; i < nthreads; i++) {
            parallelMP_threads_id[i] = i;
            pthread_create(&parallelMP_thread[i], NULL, build_output, &parallelMP_threads_id[i]);
        }

        MP_thread_pool_initialized = 1;
    } else if (nthreads != parallelMP_pool_threads) {
        return 1;
    }

    parallelMP_local_output = (long long **)malloc((size_t)nthreads * sizeof(long long *));
    if (!parallelMP_local_output) {
        return 1;
    }

    for (int t = 0; t < nthreads; t++) {
        parallelMP_local_output[t] = (long long *)calloc((size_t)nelements, sizeof(long long));
        if (!parallelMP_local_output[t]) {
            for (int j = 0; j < t; j++) {
                free(parallelMP_local_output[j]);
            }
            free(parallelMP_local_output);
            parallelMP_local_output = NULL;
            return 1;
        }
    }

    parallelHisto_data = Input;
    parallelHisto_limits = Limits;
    parallelHisto_nelements = nelements;
    parallelHisto_nbins = nbins;
    parallelHisto_active_threads = nthreads;
    parallelMP_global_pos = Pos;

    build_output(&parallelMP_threads_id[0]);

    //função que faz merge de tudo
    for(int b = 0; b < nbins; b++){
        long long count = 0;
        for(int t = 0; t < nthreads; t++){
            memcpy(
                &Output[Pos[b] + count],
                &parallelMP_local_output[t][Pos[b]],
                sizeof(long long) * parallelHisto_local_hist[t][b]
            );
            
            count += parallelHisto_local_hist[t][b];
        }
    }


    free(local_histogram);

    return 0;
}

void verifica_particoesLocais(
    const long long  *Input,       /* input array                          */
    long long        *Output,      /* Output array (partitioned per bin)   */
    long long         nElements,   /* number of elements                   */
    const long long  *Limits,      /* bin boundaries, size nbins+1         */
    int               nbins,       /* number of bins                       */
    long long        *Pos,         /* output histogram, size nbins         */
    int               nthreads     /* number of threads (1 = serial path)  */
){
    int correct = 1;

    if (nbins <= 0 || nElements < 0 || !Limits || !Pos || !Output) {
        fprintf(stderr, "      ===> particionamento local COM ERROS no rank 0\n");
        local_partition_verify_ok = 0;
        return;
    }

    if (Pos[0] != 0) {
        correct = 0;
    }

    for (int b = 1; b < nbins; b++) {
        if (Pos[b] < 0 || Pos[b] < Pos[b - 1] || Pos[b] > nElements) {
            correct = 0;
            break;
        }
    }

    if (Pos[nbins - 1] > nElements) {
        correct = 0;
    }

    if (correct) {
        for (int b = 0; b < nbins; b++) {
            int lastBin = (b == nbins - 1);
            long long start = Pos[b];
            long long end = lastBin ? nElements : Pos[b + 1];
            long long low = Limits[b];
            long long high = Limits[b + 1];

            if (start < 0 || end < start || end > nElements) {
                correct = 0;
                break;
            }

            for (long long i = start; i < end; i++) {
                long long value = Output[i];
                if (lastBin) {
                    if (value < low || value > high) {
                        correct = 0;
                        break;
                    }
                } 
                else {
                    if (value < low || value >= high) {
                        correct = 0;
                        break;
                    }
                }
            }

            if (!correct) {
                break;
            }
        }
    }

    if (correct) {
        printf("      ===> particionamento local CORRETO no rank 0\n");
        local_partition_verify_ok = 1;
    } 
    else {
        printf("      ===> particionamento local COM ERROS no rank 0\n");
        local_partition_verify_ok = 0;
    }
}

void verifica_particoesGlobais(
    const long long  *Input,       /* input array                          */
    long long        *Output,      /* Output array (partitioned per bin)   */
    long long         nElements,   /* number of elements                   */
    const long long  *Limits,      /* bin boundaries, size nbins+1         */
    int               nbins,       /* number of bins                       */
    long long        *Pos,         /* output histogram, size nbins         */
    int               nthreads     /* number of threads (1 = serial path)  */
){
    verifica_particoesLocais(
        Input,
        Output,
        nElements,
        Limits,
        nbins,
        Pos,
        nthreads
    );
}

int validateInputs(){
    if(nelements <= 0){
        printf("nelements should be bigger than 0, currently: %lld\n", nelements);
        return 1;
    }

    if(nbins < 1){
        printf("nbins should be bigger than 1, currently: %lld\n", nbins);
        return 1;
    }

    if(npivots < nbins || npivots < 2 || npivots > nelements){
        printf("npivots should be bigger than 2 and %lld, and smaller than %lld, currently: %lld\n", nbins, nelements, npivots);
        return 1;
    }

    if(nthreads <= 0 || nthreads > MAX_THREADS){
        printf("nthreads should be between 1 and %d, currently: %lld\n", MAX_THREADS, nthreads);
        return 1;
    }

    if(nr <= 0){
        printf("nr should be bigger than 0, currently: %lld\n", nr);
        return 1;
    }

    return 0;
}

int main(int argc, char* argv[]){
    if(argc < 6 || argc > 8){
        printf("usage: %s <nelements> <npivots> <nbins> <nthreads> <nr> [-Verify] [-tb2] \n", argv[0]);
        return 0;
    }

    int num_proc;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    nelements = atoll(argv[1]);
    npivots = atoi(argv[2]);
    nbins = atoi(argv[3]);
    nthreads = atoi(argv[4]);
    nr = atoi(argv[5]);

    use_verify = (argc == 7 && strcmp(argv[6], "-Verify") == 0);
    use_tb2 = (argc == 7 && strcmp(argv[6], "-tb2") == 0) || (argc == 8 && strcmp(argv[7], "-tb2") == 0);
    
    if(validateInputs()){
        return 0;
    }

    srand(2025*100*rank);

    long long nelements_local = nelements / num_proc;

    evictCacheInit();

    if(rank == 0) {
        printf("\n=== Parallel Multi-Partition — Scalability Test (Persistent Thread Pool and MPI Comunication) ===\n");
        printf("  Elements : %lld  |  Pivots : %lld  |  Bins : %lld  |  Threads : %lld  |  Process's : %i  |  Rounds : %lld  |  Verify : %s  |  Tb2: %s\n",
            nelements, npivots, nbins, nthreads, num_proc, nr,
            use_verify ? "Able" : "Disable",
            use_tb2 ? "Able" : "Disable");
        printf("  LLC size : %d MiB  |  Eviction buffer : %d MiB\n\n",
            CACHE_SIZE_MB, CACHE_SIZE_MB * EVICTION_MULTIPLIER);
    }

    double *t_1thr1proc = (double *)malloc(nr * sizeof(double));
    double *t_NthrNPproc = (double *)malloc(nr * sizeof(double));
    double *spdup  = (double *)malloc(nr * sizeof(double));
    int    *ok_arr = (int *)   malloc(nr * sizeof(int));
    int all_ok = 1;

    chronometer_t sproc_chronometer;
    chronometer_t mproc_chronometer;

    if(rank == 0) {
        printf("  Round  ;  T(part-ser) s  ;  T(np:%i x nth:%lld) s    ;  Speedup  ;  OK?\n", num_proc, nthreads);
        printf(" ------- ; --------------- ; --------------------- ; --------- ; ----\n");
    }

    for(int r = 0; r < nr; r++){
        long long *data = (long long *)malloc(sizeof(long long) * nelements * 2);
        if (!data) { 
            fprintf(stderr, "malloc failed\n"); 
            return 1; 
        }

        if (use_tb2) {
            gen_test_data_balanced2(data, nelements_local, nbins);
        } 
        else {
            for (long long i = 0; i < nelements_local; i++) {
                data[i] = rand64();
            }
        }

        long long *pivots = (long long *)malloc(sizeof(long long) * npivots);
        long long *limits = (long long *)malloc(sizeof(long long) * (nbins + 1));
        long long *Output_1 = (long long *)calloc(nelements, sizeof(long long));
        long long *Output_n = (long long *)calloc(nelements, sizeof(long long));
        long long *Pos      = (long long *)malloc(sizeof(long long) * nbins);
        long long * final_output = (long long*)malloc(sizeof(long long) * nelements);
        int nelements_final_output = 0;

        if (!pivots || !limits || !Output_1 || !Output_n || !Pos || !final_output) {
            fprintf(stderr, "Memory allocation failed\n");
            return 1;
        }

        chrono_reset(&sproc_chronometer);
        chrono_reset(&mproc_chronometer);
        
        build_limits_mpi(data, nelements, npivots, nbins, limits, MPI_COMM_WORLD);
        
        // Execução 1 thread
        evictCache();

        long long * Input_all = NULL;
        if(rank == 0) {
            Input_all = malloc(sizeof(long long) * nelements);
        }

        MPI_Gather(
            data,
            nelements_local,
            MPI_LONG_LONG,

            Input_all,
            nelements_local,
            MPI_LONG_LONG,

            0,
            MPI_COMM_WORLD
        );

        long long *Output_serial = NULL;
        long long *Pos_serial    = NULL;
        if(rank == 0) {
            Output_serial = calloc(nelements, sizeof(long long));
            Pos_serial    = calloc(nbins, sizeof(long long));

            chrono_start(&sproc_chronometer);

            parallel_multiPartition(
                Input_all,
                Output_serial,
                nelements,
                limits,
                nbins,
                Pos_serial,
                1
            );

            chrono_stop(&sproc_chronometer);

            free(Output_serial);
            free(Pos_serial);
            free(Input_all);
        }
        else{
            chrono_start(&sproc_chronometer);
            chrono_stop(&sproc_chronometer);
        }
        
        // Execução paralela nthread
        evictCache();

        chrono_start(&mproc_chronometer);
        parallel_multiPartition(data, Output_n, nelements_local, limits, nbins, Pos, nthreads);

        final_output = multi_partition_mpi(Output_n, Pos, nelements_local, nbins, &nelements_final_output, MPI_COMM_WORLD);        
        chrono_stop(&mproc_chronometer);

        local_partition_verify_ok = 1;
        if (use_verify) {
            verifica_particoesLocais(data, Output_n, nelements_local, limits, nbins, Pos, nthreads);
        }
        ok_arr[r] = local_partition_verify_ok;
        if(!ok_arr[r]) {
            all_ok = 0;
        }

        t_1thr1proc[r] = (double) chrono_gettotal(&sproc_chronometer) / (double)1e9;
        t_NthrNPproc[r] = (double) chrono_gettotal(&mproc_chronometer) / (double)1e9;

        spdup[r] = t_1thr1proc[r] / t_NthrNPproc[r];

        if(rank == 0) {
            printf("  %-5d  ;   %12.6f  ;         %12.6f  ;%9.3f  ;  %s\n",
                  r + 1, t_1thr1proc[r], t_NthrNPproc[r], spdup[r],
                  ok_arr[r] ? "OK" : "FAIL");
            verifica_particoesGlobais(data, final_output, nelements, limits, nbins, Pos, nthreads);
        }

        free(data);
        free(pivots);
        free(limits);
        free(Output_1);
        free(Output_n);
        free(final_output);
        chrono_reset(&sproc_chronometer);
        chrono_reset(&mproc_chronometer);
    }

    double avg_sproc = 0, avg_mproc = 0, avg_sp = 0;
    for (int r = 0; r < nr; r++) {
        avg_sproc += t_1thr1proc[r];
        avg_mproc += t_NthrNPproc[r];
        avg_sp += spdup[r];
    }
    avg_sproc /= nr;
    avg_mproc /= nr;
    avg_sp /= nr;

    if(rank == 0) {
        printf(" ------- ; --------------- ; --------------------- ; --------- ; ----\n");
        printf("   AVG   ;  %12.6f   ;         %12.6f  ;%9.3f  ;  %s\n\n",
               avg_sproc, avg_mproc, avg_sp, all_ok ? "OK" : "FAIL");

        double meps_sp = (nelements / 1e6) / avg_sproc;
        double meps_mp = (nelements / 1e6) / avg_mproc;
    
        printf(" ============================== Summary =============================\n");
        printf("  Avg time  (1 process 1 thread )     : %.6f ; s   ; %.2f ; MEPS\n", avg_sproc, meps_sp);
        printf("  Avg time  (%i process %lld threads)     : %.6f ; s   ; %.2f ; MEPS\n", num_proc, nthreads, avg_mproc, meps_mp);
        printf("  Avg Partition problem speedup               : %.3fx\n\n", avg_sp);
        printf("  Overall correctness                 : %s\n\n", all_ok ? "PASS" : "FAIL");
    }

    free(t_1thr1proc);
    free(t_NthrNPproc);
    free(spdup);
    free(ok_arr);
    free((void *)eviction_buffer);


    MPI_Finalize();
    return 0;
}
