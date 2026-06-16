#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#include "chrono.c"

#define MAX_THREADS 8
#define CACHE_SIZE_MB 6
#define EVICTION_MULTIPLIER 3

// variaveis de entrada
long long int nelements = 0, 
              npivots = 0, 
              nbins = 0, 
              nthreads = 0, 
              nr = 0;
static int use_verify = 0;

// flag para evitar a criação de threads
static int thread_pool_initialized = 0;

static double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// variáveis de gerenciamento de pool
pthread_barrier_t parallelHisto_barrier;
int parallelHisto_threads_id[MAX_THREADS];
pthread_t parallelHisto_thread[MAX_THREADS];
static int parallelHisto_pool_threads = 0;

// variáveis de dados
static const long long *parallelHisto_data = NULL;
static const long long *parallelHisto_limits = NULL;
static long long **parallelHisto_local_hist = NULL;
static long long parallelHisto_nelements = 0;
static int parallelHisto_nbins = 0;
static int parallelHisto_active_threads = 1;


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


static void gen_test_data_balanced2(long long *data,
                                     long long nelements,
                                     int nbins)
{
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

/* Pre-allocated eviction buffer */
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

static void build_limits_sp2_serial(const long long *data, 
                                    long long n,
                                    int npivots, 
                                    int nbins,
                                    long long *pivots, 
                                    long long *limits )
{
    unsigned long long stride = n / npivots;

    for(int i=0; i<npivots; i++){
        long long jitter_i = (long long)(rand63() % stride);
        pivots[i] = data[i*stride + jitter_i];
    }

    qsort(pivots, npivots, sizeof(long long), cmp_ll);

    limits[0] = LLONG_MIN;
    limits[nbins] = LLONG_MAX;

    int step = npivots/(nbins - 1);

    for (int k = 1; k < nbins; k++){
        limits[k] = pivots[k*step];
    }

    for (int k = 1; k <= nbins - 1; k++){
        if (limits[k] <= limits[k - 1]){
            if (limits[k - 1] < LLONG_MAX - 1)
                limits[k] = limits[k - 1] + 1;
            else
                limits[k] = limits[k - 1];
        }
    }
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

    if (!thread_pool_initialized) {
        pthread_barrier_init(&parallelHisto_barrier, NULL, nthreads);

        parallelHisto_pool_threads = nthreads;
        parallelHisto_threads_id[0] = 0;
        for (int i = 1; i < nthreads; i++) {
            parallelHisto_threads_id[i] = i;
            pthread_create(&parallelHisto_thread[i], NULL, histogram, &parallelHisto_threads_id[i]);
        }

        thread_pool_initialized = 1;
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

    for (int t = 0; t < nthreads; t++) {
        free(parallelHisto_local_hist[t]);
    }
    free(parallelHisto_local_hist);
    parallelHisto_local_hist = NULL;

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
    long long * next            = (long long *)malloc(nbins * sizeof(long long));

    parallelHistogram(Input, nElements, Limits, nbins, local_histogram, nthreads);
    prefix_sum(Pos, local_histogram, nbins);

    memcpy(next, Pos, sizeof(long long) * nbins);
    memset(next, 0, sizeof(long long) * nbins);

    for(int i = 0; i < nElements; i++){
        int bin = findBin(Limits, nbins, Input[i]);
        Output[next[bin]] = Input[i];
        next[bin]++;
    }

    free(local_histogram);
    free(next);

    return 0;
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


static int verifyHistogram(
    const long long *data,
    long long        nelements,
    const long long *limits,
    int              nbins,
    const long long *hist_1thr,
    const long long *hist_nthr)
{
    /* Stage 1: 1-thread vs N-thread */
    int s1_ok = 1;
    for (int b = 0; b < nbins; b++) {
        if (hist_1thr[b] != hist_nthr[b]) {
            fprintf(stderr,
                    "  VERIFY FAIL stage1: bin %d  1thr=%lld  Nthr=%lld\n",
                    b, hist_1thr[b], hist_nthr[b]);
            s1_ok = 0;
        }
    }
    if (!s1_ok) return 0;

    /* Stage 2: serial recount vs hist_nthr */
    long long *recount = (long long *)calloc(nbins, sizeof(long long));
    if (!recount) { perror("calloc recount"); return 0; }

    for (long long i = 0; i < nelements; i++) {
        long long v = data[i];
        int b = 0;
        while (b < nbins - 1 && v >= limits[b + 1]) b++;
        recount[b]++;
    }

    int s2_ok = 1;
    for (int b = 0; b < nbins; b++) {
        if (recount[b] != hist_nthr[b]) {
            fprintf(stderr,
                    "  VERIFY FAIL stage2: bin %d  recount=%lld  Nthr=%lld\n",
                    b, recount[b], hist_nthr[b]);
            s2_ok = 0;
        }
    }
    free(recount);
    if (!s2_ok) return 0;

    /* Stage 3: sum of bins == nelements */
    long long total = 0;
    for (int b = 0; b < nbins; b++) total += hist_nthr[b];
    if (total != nelements) {
        fprintf(stderr,
                "  VERIFY FAIL stage3: sum=%lld expected=%lld\n",
                total, nelements);
        return 0;
    }

    return 1;
}

int main(int argc, char* argv[]){
    if(argc < 6 || argc > 7){
        printf("usage: %s <nelements> <npivots> <nbins> <nthreads> <nr> [-Verify]\n", argv[0]);
        return 0;
    }

    nelements = atoll(argv[1]);
    npivots = atoi(argv[2]);
    nbins = atoi(argv[3]);
    nthreads = atoi(argv[4]);
    nr = atoi(argv[5]);

    use_verify = (argc == 7 && strcmp(argv[6], "-Verify") == 0);
    
    if(validateInputs()){
        return 0;
    }

    srand(2025 * 100 + 1); // +1 mas tem que deixar +i

    /* Init eviction buffer */
    evictCacheInit();

    /* Print header */
    printf("\n=== Parallel Multi-Partition — Scalability Test (Persistent Thread Pool and MPI Comunication) ===\n");
    printf("  Elements : %lld  |  Pivots : %lld  |  Bins : %lld  |  Threads : %lld  |  Rounds : %lld  |  Verify : %s\n",
           nelements, npivots, nbins, nthreads, nr,
           use_verify ? "Able" : "Disable");
    printf("  LLC size : %d MiB  |  Eviction buffer : %d MiB\n\n",
           CACHE_SIZE_MB, CACHE_SIZE_MB * EVICTION_MULTIPLIER);

    /* Arrays for per-round times */
    double *t_bl   = (double *)malloc(nr * sizeof(double));
    double *t_1thr = (double *)malloc(nr * sizeof(double));
    double *t_nthr = (double *)malloc(nr * sizeof(double));
    double *spdup  = (double *)malloc(nr * sizeof(double));
    int    *ok_arr = (int *)   malloc(nr * sizeof(int));
    int all_ok = 1;

    chronometer_t bl_chronometer;
    chronometer_t sthr_chronometer;
    chronometer_t nthr_chronometer;

    /* Print round table header */
    printf("  Round ;  T(bl_ser) s ;  T(1 thr) s ;   T(N thr) s ;    Speedup ; OK?\n");
    printf("  ----- ;------------ ;------------ ;------------ ;------------ ;---------- ; ----\n");

    for(int r = 0; r < nr; r++){
        /* Generate data */
        long long *data = (long long *)malloc(sizeof(long long) * nelements);
        if (!data) { fprintf(stderr, "malloc failed\n"); return 1; }

        gen_test_data_balanced2(data, nelements, nbins);

        long long *pivots = (long long *)malloc(sizeof(long long) * npivots);
        long long *limits = (long long *)malloc(sizeof(long long) * (nbins + 1));
        long long *Output_1 = (long long *)calloc(nelements, sizeof(long long));
        long long *Output_n = (long long *)calloc(nelements, sizeof(long long));
        long long *Pos      = (long long *)malloc(sizeof(long long) * nbins);

        if (!pivots || !limits || !Output_1 || !Output_n || !Pos) {
            fprintf(stderr, "Memory allocation failed\n");
            return 1;
        }

        chrono_reset(&bl_chronometer);
        chrono_reset(&sthr_chronometer);
        chrono_reset(&nthr_chronometer);

        chrono_start(&bl_chronometer);
        build_limits_sp2_serial(data, nelements, npivots, nbins, pivots, limits);
        chrono_stop(&bl_chronometer);

        /* Print 8 bins */
/*        if (r == 0) {
            int show = nbins < 8 ? nbins : 8;
            printf("\n  --- Round 1: first %d partitions ---\n", show);
            printf("   Bin  ;          Lo (inclusive)  ;          Hi (exclusive)  ;         Count\n");

            long long *tmp_hist = (long long *)calloc(nbins, sizeof(long long));
            for (long long i = 0; i < nelements; i++) {
                int bin = findBin(limits, nbins, data[i]);
                tmp_hist[bin]++;
            }
            for (int b = 0; b < show; b++) {
                const char *hi_label = (b == nbins - 1) ? "" : "";
                printf("   %3d  ;  %22lld  ;  %22lld  ;  %12lld\n",
                       b, limits[b], limits[b + 1], tmp_hist[b]);
            }
            if (nbins > 8)
                printf("  ... (%lld more bins not shown)\n", nbins - 8);
            free(tmp_hist);
            printf("\n");
            printf("  Round ;  T(bl_ser) s ;  T(1 thr) s ;   T(N thr) s ;    Speedup ; OK?\n");
            printf("  ----- ;------------ ;------------ ;------------ ;------------ ;---------- ; ----\n");
        }
*/
        /* Limpar cache antes de exec de 1-thread */
        evictCache();

        /* dados para testes */
        long long Input [14] = {8, 4, 13, 7, 11, 100, 44, 3, 7, 7, 100, 110, 46, 44};
        long long Limits_teste[5] = {LLONG_MIN, 12, 70, 90, LLONG_MAX};
        nbins = 4;

        /* Histograma 1-thread */
        chrono_start(&sthr_chronometer);
        parallel_multiPartition(data, Output_1, nelements, limits, nbins, Pos, 1);
        chrono_stop(&sthr_chronometer);

        /* saida dos testes */
/*        printf("Output: [");
        for(int i = 0; i < nelements; i++)
            printf("%lld ", Output_1[i]);
        printf("]\n");

        printf("Pos: [");
        for(int i = 0; i < nbins; i++)
            printf("%lld ", Pos[i]);
        printf("]\n");
*/
        /* Limpar cache antes de exec de N-thread */
        evictCache();

        /* Histograma N-thread*/
        chrono_start(&nthr_chronometer);
        parallel_multiPartition(data, Output_n, nelements, limits, nbins, Pos, nthreads);
        chrono_stop(&nthr_chronometer);

        /* saida dos testes */
/*        printf("Output: [");
        for(int i = 0; i < nelements; i++)
            printf("%lld ", Output_1[i]);
        printf("]\n");

        printf("Pos: [");
        for(int i = 0; i < nbins; i++)
            printf("%lld ", Pos[i]);
        printf("]\n");
*/

        /* Checar */
//        ok_arr[r] = verifyHistogram(data,Helements, limits, nbins, hist_1, hist_n);
//        if (!ok_arr[r]) all_ok = 0;

        t_bl[r]   = (double) chrono_gettotal(&bl_chronometer) / (double)1e9;
        t_1thr[r] = (double) chrono_gettotal(&sthr_chronometer) / (double)1e9;
        t_nthr[r] = (double) chrono_gettotal(&nthr_chronometer) / (double)1e9;

        spdup[r] = t_1thr[r] / t_nthr[r];

        printf("  %-5d ;  %12.6f ;  %12.6f ;  %12.6f ;  %9.3f ; %s\n",
               r + 1, t_bl[r], t_1thr[r], t_nthr[r], spdup[r],
               ok_arr[r] ? "OK" : "FAIL");

        free(data);
        free(pivots);
        free(limits);
        free(Output_1);
        free(Output_n);
    }

    /* Compute averages */
/*    double avg_bl = 0, avg_1 = 0, avg_n = 0, avg_sp = 0;
    for (int r = 0; r < nr; r++) {
        avg_bl += t_bl[r];
        avg_1  += t_1thr[r];
        avg_n  += t_nthr[r];
        avg_sp += spdup[r];
    }
    avg_bl /= nr;
    avg_1  /= nr;
    avg_n  /= nr;
    avg_sp /= nr;

    printf("  ----- ;------------ ;------------ ;------------ ;------------ ;---------- ; ----\n");
    printf("  AVG   ;  %12.6f ;  %12.6f ;  %12.6f ;  %9.3f ; %s\n\n",
           avg_bl, avg_1, avg_n, avg_sp, all_ok ? "OK" : "FAIL");

    double meps_1 = (nelements / 1e6) / avg_1;
    double meps_n = (nelements / 1e6) / avg_n;
    double pct_bl_1 = (avg_bl / avg_1) * 100.0;
    double pct_bl_n = (avg_bl / avg_n) * 100.0;
    double efficiency = (avg_sp / nthreads) * 100.0;

    printf("=== Summary ===\n");
    printf("  Avg build_limits serial   : %.6f ; s   ; %5.1f%% ; of T(1 thr) | %5.1f%% ; of T(N thr)\n\n",
           avg_bl, pct_bl_1, pct_bl_n);
    printf("  Avg time  (1 thread )     : %.6f ; s   ; %.2f ; MEPS\n", avg_1, meps_1);
    printf("  Avg time  (%lld threads)    : %.6f ; s   ; %.2f ; MEPS\n", nthreads, avg_n, meps_n);
    printf("  Avg histogram speedup     : %.3fx\n\n", avg_sp);
    printf("  Parallel efficiency:\n");
    printf("    with nthreads   (%2lld) : %5.1f%%\n\n", nthreads, efficiency);
    printf("  Overall correctness       : %s\n\n", all_ok ? "PASS" : "FAIL");

    free(t_bl);
    free(t_1thr);
    free(t_nthr);
    free(spdup);
    free(ok_arr);
    free((void *)eviction_buffer);
*/
    return 0;
}
