#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

/* ----------------------------------------------------
 * vector of int
 * ---------------------------------------------------
 */
typedef struct {
    int *data;
    size_t size;
    size_t capacity;
} vint_t;

void
vint_alloc(vint_t *self, size_t size)
{
    assert(size >= 0);
    self->size = size;
    self->capacity = size;
    if (size > 0) {
        self->data = malloc(sizeof(*self->data) * size);
    } else {
        self->data = NULL;
    }
}

void
vint_free(vint_t *self)
{
    if (self && self->data) {
        free(self->data);
        self->data = NULL;
    }
}

void
vint_initialize(vint_t *self, int init_value)
{
    for (int i = 0; i < self->size; i++) {
        self->data[i] = init_value;
    }
}

extern inline void
vint_add_element(vint_t *self, int value)
{
    // fprintf(stderr, "%ld\n", self->size);
    if (self->size == self->capacity) {
        if (self->capacity == 0) {
            self->capacity = 1;
        } else {
            self->capacity *= 2;
        }
        self->data = realloc(self->data, sizeof(*self->data) * self->capacity);
        assert(self->data != NULL);
    }
    self->size += 1;
    self->data[self->size - 1] = value;
}

/* ----------------------------------------------------
 * vector of double
 * ---------------------------------------------------
 */
typedef struct {
    double *data;
    size_t size;
    size_t capacity;
} vdouble_t;

void
vdouble_alloc(vdouble_t *self, size_t size)
{
    assert(size >= 0);
    self->size = size;
    self->capacity = size;
    if (size > 0) {
        self->data = malloc(sizeof(*self->data) * size);
    } else {
        self->data = NULL;
    }
}

void
vdouble_free(vdouble_t *self)
{
    if (self && self->data) {
        free(self->data);
        self->data = NULL;
    }
}

void
vdouble_init(vdouble_t *self, double init_value)
{
    for (int i = 0; i < self->size; i++) {
        self->data[i] = init_value;
    }
}

extern inline void
vdouble_add_element(vdouble_t *self, double value)
{
    if (self->size == self->capacity) {
        if (self->capacity == 0) {
            self->capacity = 1;
        } else {
            self->capacity *= 2;
        }
        self->data = realloc(self->data, sizeof(*self->data) * self->capacity);
        assert(self->data != NULL);
    }
    self->size += 1;
    self->data[self->size - 1] = value;
    // fprintf(stderr, "aded %ld, %g\n", self->size - 1, self->data[self->size - 1]);
}

void
vdouble_fprintf(vdouble_t *self, FILE *fp)
{
    assert(fp != NULL);
    for (size_t i = 0; i < self->size; i++) {
        fprintf(fp, "%ld\t%g\n", i, self->data[i]);
    }
}

void
vdouble_fprintf_range(vdouble_t *self, FILE *fp, size_t start, size_t end)
{
    assert(fp != NULL);
    assert(start <= end);
    assert(end < self->size);
    for (size_t i = start; i <= end; i++) {
        fprintf(fp, "%ld\t%g\n", i, self->data[i]);
    }
}

/* ----------------------------------------------------
 * struct for multithreading
 * ---------------------------------------------------
 */
typedef struct {
    int num_threads;
    pthread_mutex_t mutex;
    pthread_t *threads;
    long running_id;
    long last_id;
    void *(*thread_func)(void *);
} tpool_t;

void
tpool_alloc(tpool_t *self, int num_threads)
{
    if (num_threads < 1)
        num_threads = 1;
    self->num_threads = num_threads;
    self->threads = malloc(sizeof(*self->threads) * num_threads);
    assert(self->threads != NULL);
}

void
tpool_run(tpool_t *self, long thread_init_job_id, long thread_last_job_id,
    void *(*thread_func)(void *), void *param)
{
    self->thread_func = thread_func;
    self->running_id = thread_init_job_id;
    self->last_id = thread_last_job_id;

    fprintf(stderr, "num_thread: %d\n", self->num_threads);
    for (int i = 0; i < self->num_threads; i++) {
        pthread_create(self->threads + i, NULL, thread_func, (void *) param);
    }
    for (int i = 0; i < self->num_threads; i++) {
        pthread_join(self->threads[i], NULL);
    }
}

void
tpool_free(tpool_t *self)
{
    if (self->threads != NULL) {
        free(self->threads);
        self->threads = NULL;
        self->num_threads = 0;
    }
}

/* ----------------------------------------------------
 * struct for ibdqc
 * ---------------------------------------------------
 */

typedef struct {
    // NE
    vdouble_t ne_traj;
    // ibd
    vdouble_t ibd_cM;
    vint_t ibd_num_ends;
    vdouble_t ibd_consts;
    // genome
    vint_t chr_lens;
    size_t chr_end_range;
    // constants
    unsigned int g_star;
    unsigned int g_max;
    // multithreads
    tpool_t tpool;
    // results
    vdouble_t total_ibd_due_to_tmrca;
} ibdqc_t;

void
ibdqc_alloc(
    ibdqc_t *self, unsigned int g_max, unsigned int g_star, int refsize, int num_threads)
{
    self->g_max = g_max;
    self->g_star = g_star;

    // TODO: put this in function parameters
    // how many bp from between the ends of chromosome and ends of ibdsegment
    // are considered reaching the ends.
    self->chr_end_range = 100;

    vdouble_alloc(&self->ne_traj, g_max + 1);
    vdouble_init(&self->ne_traj, refsize);

    vdouble_alloc(&self->total_ibd_due_to_tmrca, g_max + 1);
    vdouble_init(&self->total_ibd_due_to_tmrca, 0);

    vdouble_alloc(&self->ibd_cM, 0);
    vint_alloc(&self->ibd_num_ends, 0);
    vdouble_alloc(&self->ibd_consts, 0);
    vint_alloc(&self->chr_lens, 0);

    tpool_alloc(&self->tpool, num_threads);
}

void
ibdqc_free(ibdqc_t *self)
{
    vdouble_free(&self->ne_traj);
    vdouble_free(&self->ibd_cM);
    vdouble_free(&self->ibd_consts);
    vint_free(&self->ibd_num_ends);
    vint_free(&self->chr_lens);

    vdouble_free(&self->total_ibd_due_to_tmrca);
    tpool_free(&self->tpool);
}

// calc gamma(N, l)
// run after NE is set and IBD length and ends are determined
void *
ibdqc_thread_func_calc_prob_consts(void *self_void)
{
    /* variables */
    long ibd_index;
    double sum, prod, tmp, l;
    int num_ends;
    long job_id;

    /* parameters or constant */
    ibdqc_t *self = (ibdqc_t *) self_void;
    int G = self->g_max;
    int g_star = self->g_star;
    double *N = self->ne_traj.data;
    assert(G == self->ne_traj.size - 1);
    double alpha; // = l / 50.0;
    double beta = 1 - 0.5 / N[G];
    int step = 1000;

    while (1) {
        // get job id
        pthread_mutex_lock(&self->tpool.mutex);
        job_id = self->tpool.running_id;
        if (job_id < self->tpool.last_id) {
            pthread_mutex_unlock(&self->tpool.mutex);
            return NULL;
        }
        self->tpool.running_id -= step;
        pthread_mutex_unlock(&self->tpool.mutex);

        // do the job: each job do `step` loops
        for (ibd_index = job_id;
             ibd_index > job_id - step && ibd_index >= self->tpool.last_id;
             ibd_index--) {
            l = self->ibd_cM.data[ibd_index];
            num_ends = self->ibd_num_ends.data[ibd_index];
            alpha = l / 50.0;

            /* calc const_gamma */
            sum = 0;
            for (int g1 = g_star; g1 <= G; g1++) {
                prod = 1;
                for (int i = 0; i < num_ends; i++) {
                    prod *= g1 / 50.0;
                }
                prod *= exp(-0.02 * l * g1);
                for (int g2 = g_star; g2 < g1; g2++) {
                    prod *= 1 - 0.5 / N[g2];
                }
                prod *= 0.5 / N[g1];
                sum += prod;
            }

            prod = 1;
            for (int g1 = g_star; g1 <= G; g1++) {
                prod *= 1 - 0.5 / N[g1];
            }
            prod *= (1 - beta) * exp(-alpha * (G + 1)) / 2500;
            tmp = 1 - beta * exp(-alpha);
            if (num_ends == 0) {
                prod *= G * G / tmp + (2 * G - 1) / tmp / tmp + 2 / tmp / tmp / tmp;
            } else if (num_ends == 1) {
                prod *= G / tmp + 1 / tmp / tmp;
            } else {
                // if num_ends_reached == 2) prod = prod * 1
            }
            sum += prod;

            self->ibd_consts.data[ibd_index] = sum;

            // fprintf(stderr, "l = %g, num_ends = %d, gamma = %g\n", l,  num_ends, sum);
        }
    }
    return NULL;
}

extern inline double
ibdqc_calc_prob_ibd_due_to_g(
    ibdqc_t *self, long ibd_index, unsigned int g, int num_ends_reached)
{
    /* variables */
    double sum = 0, prod = 0, tmp = 0;
    double const_gamma = 0;
    double probability = 0;

    /* parameters */
    unsigned int G = self->g_max;
    unsigned int g_star = self->g_star;
    double *N = self->ne_traj.data;
    double l = self->ibd_cM.data[ibd_index];
    assert(G == self->ne_traj.size - 1);
    assert(num_ends_reached >= 0 && num_ends_reached <= 2);
    const_gamma = self->ibd_consts.data[ibd_index];

    /* calc probability due to g given ibd length */
    prod = 1;
    prod *= 1 / const_gamma;
    for (int i = 0; i < num_ends_reached; i++) {
        prod *= g / 50.0;
    }
    prod *= exp(-l * g / 50);
    for (int g1 = g_star; g1 < g; g1++) {
        prod *= 1 - 0.5 / N[g1];
    }
    prod *= 0.5 / N[g];

    probability = prod;

    return probability;
}

void *
ibdqc_thread_func_calc_total_ibd_due_to_g(void *selfp)
{
    ibdqc_t *self = (ibdqc_t *) selfp;
    int g = -1;

    while (1) {
        pthread_mutex_lock(&self->tpool.mutex);
        g = self->tpool.running_id;
        if (g < self->tpool.last_id) {
            pthread_mutex_unlock(&self->tpool.mutex);
            return NULL;
        } else {
            self->tpool.running_id -= 1;
        }
        pthread_mutex_unlock(&self->tpool.mutex);

        self->total_ibd_due_to_tmrca.data[g] = 0;
        for (size_t i = 0; i < self->ibd_num_ends.size; i++) {
            assert(self->ibd_num_ends.size == self->ibd_cM.size);
            double l = self->ibd_cM.data[i];
            self->total_ibd_due_to_tmrca.data[g]
                += l
                   * ibdqc_calc_prob_ibd_due_to_g(
                       self, i, g, self->ibd_num_ends.data[i]);
        }
        fprintf(stderr, "g = %d, total = %g\n", g, self->total_ibd_due_to_tmrca.data[g]);
    }
}

void
ibdqc_read_chr_length(ibdqc_t *self, FILE *fp_chr_length)
{
    size_t size = 0;
    char *buff = NULL;
    while (getline(&buff, &size, fp_chr_length) > 0) {
        size_t chr_len;
        chr_len = atol(buff);
        assert(chr_len > 0);
        vint_add_element(&self->chr_lens, chr_len);
    }
}

void
ibdqc_read_ibd(ibdqc_t *self, FILE *fp_ibd)
{
    // line buffer
    size_t size = 0;
    char *buff = NULL;

    // fields and counts
    int chr;
    size_t start;
    size_t end;
    double ibd_cM;
    int num_fields;
    char *tok;
    // num_ends
    int num_ends;

    while (getline(&buff, &size, fp_ibd) > 0) {
        tok = strtok(buff, "\t\n ");
        num_fields = 1;
        for (; num_fields <= 8; num_fields++) {
            assert(tok != NULL);
            if (num_fields == 5) {
                chr = atoi(tok);
                assert(chr - 1 < self->chr_lens.size);
            } else if (num_fields == 6) {
                start = atoi(tok);
            } else if (num_fields == 7) {
                end = atoi(tok);
            } else if (num_fields == 8) {
                ibd_cM = strtod(tok, NULL);
            } else {
            }
            tok = strtok(NULL, "\t\n ");
        }
        num_ends = 0;
        if (start - 0 <= self->chr_end_range) {
            num_ends++;
        }
        if (self->chr_lens.data[chr] - end <= self->chr_end_range) {
            num_ends++;
        }
        vdouble_add_element(&self->ibd_cM, ibd_cM);
        vint_add_element(&self->ibd_num_ends, num_ends);
    }

    // alloc for consts
    vdouble_alloc(&self->ibd_consts, self->ibd_cM.size);
    vdouble_init(&self->ibd_consts, -1);
}

void
test_read_ibd_chrlen_from_file(int argc, char *argv[])
{
    time_t t;
    ibdqc_t qc;
    FILE *fp_ibd, *fp_chr_length;

    if (argc < 3) {
        perror("need to provide two file names!\n");
        exit(1);
    }

    // start time:
    time(&t);
    fprintf(stderr, "Start time: %s", ctime(&t));

    // get number of core available
    long num_processors = sysconf(_SC_NPROCESSORS_ONLN);
    fprintf(stderr, "num_processors: %ld\n", num_processors);

    ibdqc_alloc(&qc, 300, 2, 10000, num_processors);

    fp_chr_length = fopen(argv[1], "r");
    assert(fp_chr_length != NULL);

    fp_ibd = fopen(argv[2], "r");
    assert(fp_ibd != NULL);

    ibdqc_read_chr_length(&qc, fp_chr_length);
    ibdqc_read_ibd(&qc, fp_ibd);

    fclose(fp_ibd);
    fclose(fp_chr_length);

    tpool_run(&qc.tpool, qc.ibd_cM.size - 1, 0, ibdqc_thread_func_calc_prob_consts,
        (void *) &qc);

    tpool_run(
        &qc.tpool, 100, 3, ibdqc_thread_func_calc_total_ibd_due_to_g, (void *) &qc);

    // print result
    vdouble_fprintf_range(&qc.total_ibd_due_to_tmrca, stdout, 3, 100);

    // end time:
    time(&t);
    fprintf(stderr, "End time: %s", ctime(&t));

    ibdqc_free(&qc);
}

int
main(int argc, char *argv[])
{
    // test_threads();
    test_read_ibd_chrlen_from_file(argc, argv);

    return 0;
}
