#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

typedef struct {
    int *data;
    size_t size;
    size_t capacity;
} vectori_t;

void
vectori_alloc(vectori_t *self, size_t size)
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
vectori_free(vectori_t *self)
{
    if (self && self->data) {
        free(self->data);
        self->data = NULL;
    }
}

void
vectori_init(vectori_t *self, int init_value)
{
    for (int i = 0; i < self->size; i++) {
        self->data[i] = init_value;
    }
}

extern inline void
vectori_add_element(vectori_t *self, int value)
{
    // fprintf(stderr, "%ld\n", self->size);
    if (self->size == self->capacity) {
        if (self->capacity == 0) {
            self->capacity = 1;
        } else {
            self->capacity *= 2;
        }
        self->data = malloc(sizeof(*self->data) * self->capacity);
        assert(self->data != NULL);
    }
    self->size += 1;
    self->data[self->size - 1] = value;
}

typedef struct {
    double *data;
    size_t size;
    size_t capacity;
} vectord_t;

void
vectord_alloc(vectord_t *self, size_t size)
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
vectord_free(vectord_t *self)
{
    if (self && self->data) {
        free(self->data);
        self->data = NULL;
    }
}

void
vectord_init(vectord_t *self, double init_value)
{
    for (int i = 0; i < self->size; i++) {
        self->data[i] = init_value;
    }
}

extern inline void
vectord_add_element(vectord_t *self, double value)
{
    if (self->size == self->capacity) {
        if (self->capacity == 0) {
            self->capacity = 1;
        } else {
            self->capacity *= 2;
        }
        self->data = malloc(sizeof(*self->data) * self->capacity);
        assert(self->data != NULL);
    }
    self->size += 1;
    self->data[self->size - 1] = value;
}

void
vectord_fprintf(vectord_t *self, FILE *fp)
{
	assert(fp != NULL);
	for(size_t i = 0; i < self->size; i++) { 
		fprintf(fp, "%ld\t%g\n", i, self->data[i]);
	}
}

void
vectord_fprintf_range(vectord_t *self, FILE *fp, size_t start, size_t end)
{
	assert(fp != NULL);
	assert(start <= end);
	assert(end < self->size);
	for(size_t i = start; i <= end; i++) { 
		fprintf(fp, "%ld\t%g\n", i, self->data[i]);
	}
}

typedef struct {
    // NE
    vectord_t ne_traj;
    // ibd
    vectord_t ibd_cM;
    vectori_t ibd_num_ends;
    // genome
    vectori_t chr_lens;
    size_t chr_end_range;
    // constants
    unsigned int g_star;
    unsigned int g_max;
    // multithreads
    int num_threads;
    pthread_mutex_t mutex;
    int running_g;
    // results
    vectord_t total_ibd_due_to_tmrca;
} ibdqc_t;

void
ibdqc_alloc(
    ibdqc_t *self, unsigned int g_max, unsigned int g_star, int refsize, int num_threads)
{
    self->g_max = g_max;
    self->g_star = g_star;
    self->num_threads = num_threads;
    self->running_g = 100;

    // TODO: put this in function parameters
    // how many bp from between the ends of chromosome and ends of ibdsegment
    // are considered reaching the ends.
    self->chr_end_range = 100;

    vectord_alloc(&self->ne_traj, g_max + 1);
    vectord_init(&self->ne_traj, refsize);

    vectord_alloc(&self->total_ibd_due_to_tmrca, g_max + 1);
    vectord_init(&self->total_ibd_due_to_tmrca, 0);

    vectord_alloc(&self->ibd_cM, 0);
    vectori_alloc(&self->ibd_num_ends, 0);

    vectori_alloc(&self->chr_lens, 0);
}

void
ibdqc_add_ibdinfo(ibdqc_t *self, vectord_t ibd_length, vectori_t ibd_num_ends)
{
    self->ibd_cM = ibd_length;
    self->ibd_num_ends = ibd_num_ends;
}

void
ibdqc_free(ibdqc_t *self)
{
    vectord_free(&self->ne_traj);
    vectord_free(&self->ibd_cM);
    vectori_free(&self->ibd_num_ends);
    vectori_free(&self->chr_lens);

    vectord_free(&self->total_ibd_due_to_tmrca);
}

extern inline double
ibdqc_calc_prob_ibd_due_to_g(
    ibdqc_t *self, double l, unsigned int g, int num_ends_reached)
{
    /* variables */
    double sum = 0, prod = 0, tmp = 0;
    double const_gamma = 0;
    double probability = 0;

    /* parameters */
    unsigned int G = self->g_max;
    unsigned int g_star = self->g_star;
    double *N = self->ne_traj.data;
    assert(G == self->ne_traj.size - 1);
    assert(num_ends_reached >= 0 && num_ends_reached <= 2);
    double alpha = l / 50.0;
    double beta = 1 - 0.5 / N[G];

    /* calc const_gamma */
    sum = 0;
    for (int g1 = g_star; g1 <= G; g1++) {
        prod = 1;
        for (int i = 0; i < num_ends_reached; i++) {
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
    if (num_ends_reached == 0) {
        prod *= G * G / tmp + (2 * G - 1) / tmp / tmp + 2 / tmp / tmp / tmp;
    } else if (num_ends_reached == 1) {
        prod *= G / tmp + 1 / tmp / tmp;
    } else {
        // if num_ends_reached == 2) prod = prod * 1
    }
    sum += prod;

    const_gamma = sum;

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
ibdqc_thread_func(void *selfp)
{
    ibdqc_t *self = (ibdqc_t *) selfp;
    int g = -1;

    while (1) {
        pthread_mutex_lock(&self->mutex);
        g = self->running_g;
        if (g < 3) {
            pthread_mutex_unlock(&self->mutex);
            return NULL;
        } else {
            self->running_g -= 1;
        }
        pthread_mutex_unlock(&self->mutex);

        self->total_ibd_due_to_tmrca.data[g] = 0;
        for (size_t i = 0; i < self->ibd_num_ends.size; i++) {
            assert(self->ibd_num_ends.size == self->ibd_cM.size);
            double l = self->ibd_cM.data[i];
            self->total_ibd_due_to_tmrca.data[g]
                += l
                   * ibdqc_calc_prob_ibd_due_to_g(
                       self, l, g, self->ibd_num_ends.data[i]);
        }
        fprintf(stderr, "g = %d, total = %g\n", g, self->total_ibd_due_to_tmrca.data[g]);
    }
}

void
ibdqc_run_multiple_threads(ibdqc_t *self)
{
    /* threading */
    pthread_t *threads = malloc(sizeof(pthread_t) * self->num_threads);
    for (int i = 0; i < self->num_threads; i++) {
        pthread_create(threads + i, NULL, ibdqc_thread_func, (void *) self);
    }
    for (int i = 0; i < self->num_threads; i++) {
        pthread_join(threads[i], NULL);
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
        vectori_add_element(&self->chr_lens, chr_len);
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
    size_t ibd_cM;
    int num_fields;
    char *tok;
    // num_ends
    int num_ends;

    while (getline(&buff, &size, fp_ibd) > 0) {
        tok = strtok(buff, "\t\n ");
        num_fields = 1;
        num_ends = 0;
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
        if (start - 0 <= self->chr_end_range) {
            num_ends++;
        }
        if (self->chr_lens.data[chr] - end <= self->chr_end_range) {
            num_ends++;
        }
        vectord_add_element(&self->ibd_cM, ibd_cM);
        vectori_add_element(&self->ibd_num_ends, num_ends);
    }
}

void
test_threads()
{
    ibdqc_t qc;
    ibdqc_alloc(&qc, 300, 2, 10000, 6);

    size_t size = 100000;

    /* make ibd info */
    for (size_t i = 0; i < size; i++) {
        vectord_add_element(&qc.ibd_cM, rand() % 10000 / 100);
        vectori_add_element(&qc.ibd_num_ends, rand() % 2);
    }

    ibdqc_run_multiple_threads(&qc);

    ibdqc_free(&qc);
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
    fprintf(stderr, "Start time: %s",  ctime(&t));

    // get number of core available
    long num_processors = sysconf(_SC_NPROCESSORS_ONLN);
    fprintf(stderr, "num_processors: %ld\n", num_processors);

    ibdqc_alloc(&qc, 300, 2, 10000, num_processors);

    fp_chr_length = fopen(argv[1], "r");
    assert(fp_chr_length != NULL);
    // fprintf(stderr, "%s", argv[1]);
    fp_ibd = fopen(argv[2], "r");
    // fprintf(stderr, "%s", argv[2]);
    assert(fp_ibd != NULL);
    ibdqc_read_chr_length(&qc, fp_chr_length);
    ibdqc_read_ibd(&qc, fp_ibd);
    fclose(fp_ibd);
    fclose(fp_chr_length);

    ibdqc_run_multiple_threads(&qc);

    vectord_fprintf_range(&qc.total_ibd_due_to_tmrca, stdout, 3, 100);

    // end time:
    time(&t);
    fprintf(stderr, "End time: %s",  ctime(&t));

    ibdqc_free(&qc);
}

int
main(int argc, char *argv[])
{
    // test_threads();
    test_read_ibd_chrlen_from_file(argc, argv);

    return 0;
}
