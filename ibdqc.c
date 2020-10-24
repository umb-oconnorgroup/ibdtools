/*
 * ibdqc:
 *
 * 10/23/20:
 * 	`ibdqc` follows the ibdNe paper and esitmates the total amount of IBD
 * 	attributable to a TMRCA of g.
 *
 *      Ref: Browning, S., Browning, B. (2015). Accurate Non-parametric
 *      Estimation of Recent Effective Population Size from Segments of
 *      Identity by Descent The American Journal of Human Genetics  97(3),
 *      404-18. https://dx.doi.org/10.1016/j.ajhg.2015.07.012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include "tpool.h"
#include "vector.h"
#include "argtable3.h"

// to allow self access different arguments by names
typedef struct
{
	struct arg_file *file1;
	struct arg_file *file2;
	struct arg_rem * rem_stdout;
	struct arg_rem *rem_stderr;
	struct arg_lit *help;
	struct arg_end *end;
	void **argtable;
	int num_arg;
}args_t;

typedef struct {
    // arguments
    args_t args;
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
    int g_star;
    int g_max;
    // for multithreading
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
    // chr_end_range means how many bp are between the end of the chromosome and
    // the end of ibdsegment are considered reaching the ends.
    self->chr_end_range = 100;

    vdouble_alloc(&self->ne_traj, g_max + 1);
    vdouble_init(&self->ne_traj, refsize);

    vdouble_alloc(&self->total_ibd_due_to_tmrca, g_max + 1);
    vdouble_init(&self->total_ibd_due_to_tmrca, 0);

    vdouble_alloc(&self->ibd_cM, 0);
    vint_alloc(&self->ibd_num_ends, 0);
    vdouble_alloc(&self->ibd_consts, 0);
    vint_alloc(&self->chr_lens, 0);

    tpool_alloc(&self->tpool, (void *) self, num_threads);
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

    arg_freetable(self->args.argtable, self->args.num_arg);
}

/* calc gamma(N, l) for each ibd segments, with num_ends taken into consideration
 *
 * NOTE:
 * 	NE trajectory, IBD length vector and ibd num ends vector must be provided before
 * 	call this function. (This is called by tpool)
 * 	job_id represents the index to an element of IBD length vector
 */
void *
ibdqc_job_func_calc_prob_consts(void *self_void, long job_id)
{
    /* variables */
    long ibd_index;
    double sum, prod, tmp, l;
    int num_ends;

    /* parameters or constant */
    ibdqc_t *self = (ibdqc_t *) self_void;
    int G = self->g_max;
    int g_star = self->g_star;
    double *N = self->ne_traj.data;
    assert(G == self->ne_traj.size - 1);
    double alpha; // = l / 50.0;
    double beta = 1 - 0.5 / N[G];

    ibd_index = job_id;
    // do the job: each job do `step` loops
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
        // blank
    }
    sum += prod;

    self->ibd_consts.data[ibd_index] = sum;

    return NULL;
}

/*
 *  Calculate the probability of an ibd segment due to tmrca=g
 *  ibd_index: the index to an element of the ibd length vector
 */
extern inline double
ibdqc_calc_prob_ibd_due_to_g(
    ibdqc_t *self, long ibd_index, unsigned int g, int num_ends_reached)
{
    /* variables */
    double prod = 0;
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

/*
 * calculate the total ibd due to tmrca=g from all ibdsegments
 *
 * 	based on ibdqc_calc_prob_ibd_due_to_g
 * 	called by tpool
 */
void *
ibdqc_job_func_calc_total_ibd_due_to_g(void *self_void, long job_id)
{
    ibdqc_t *self = (ibdqc_t *) self_void;
    int g = job_id;

    self->total_ibd_due_to_tmrca.data[g] = 0;
    for (size_t i = 0; i < self->ibd_num_ends.size; i++) {
        assert(self->ibd_num_ends.size == self->ibd_cM.size);
        double l = self->ibd_cM.data[i];
        self->total_ibd_due_to_tmrca.data[g]
            += l * ibdqc_calc_prob_ibd_due_to_g(self, i, g, self->ibd_num_ends.data[i]);
    }
    fprintf(stderr, "g = %d, total = %g\n", g, self->total_ibd_due_to_tmrca.data[g]);
    return NULL;
}

/*
 * read in the chr length for each chr
 */
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

/*
 *  read ibd info into ibd length vector and ibd num ends vector
 *  	Note: chr_length must be determined before calling the function
 */
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
ibdqc_arg_parse(ibdqc_t *self, int argc, char *argv[])
{
    char *progname = "ibdqc";
    int ret;
    args_t *pargs = &self->args; 

    // argtable
    pargs->file1 = arg_filen(NULL, NULL, "<chr_len_file>", 1, 1,
        "chr length file: single column with integer, the nth row corresponds"
        "to the nth chr length");
    pargs->file2 = arg_filen(NULL, NULL, "<ibd_file>", 1, 1,
        "ibd file: tab/space separated 8 column files, 5th column, i.e. the"
        "chr column should be inteter and the 8tch, i.e. ibd length in cM, "
        "should be float point number");
    pargs->rem_stdout = arg_rem("1>res.txt", "stdout generates result info");
    pargs->rem_stderr = arg_rem("2>log.txt", "stderr generates log info");
    pargs->help = arg_litn("h", "help", 0, 1, "help info");
    pargs->end = arg_end(20);
    void *argtable[] = { pargs->file1, pargs->file2, pargs->help, pargs->rem_stdout, 
	    pargs->rem_stderr, pargs->end };

    pargs->num_arg = sizeof(argtable) / sizeof(argtable[0]);
    pargs->argtable = malloc(sizeof(void*) * pargs->num_arg);
    memcpy(pargs->argtable, argtable, sizeof(void*) *  pargs->num_arg);

    // check all entry successfully allocated
    assert(arg_nullcheck(argtable) == 0);

    // set default value before argparse
    pargs->file1->filename[0] = "-";
    pargs->file2->filename[0] = "-";

    // parse argument
    ret = arg_parse(argc, argv, argtable);

    // error
    if (ret > 0) {
        arg_print_errors(stdout, pargs->end, progname);
        arg_print_syntaxv(stderr, argtable, "\tIBD quality control tool\n======\n");
        arg_print_glossary(stderr, argtable, " %-25s %s\n");
    } else {
        // help message
        if (pargs->help->count == 1) {
            arg_print_syntaxv(stderr, argtable, "\tIBD quality control tool\n======\n");
            arg_print_glossary(stderr, argtable, " %-25s %s\n");

        } else {
            // correct message
            fprintf(stderr, "chr_len file: %s\n", pargs->file1->filename[0]);
            fprintf(stderr, "ibd file: %s\n", pargs->file2->filename[0]);
        }
    }
}

void
test_read_ibd_chrlen_from_file(int argc, char *argv[])
{
    time_t t;
    ibdqc_t qc;
    FILE *fp_ibd, *fp_chr_length;

    // get number of core available
    long num_processors = sysconf(_SC_NPROCESSORS_ONLN);
    fprintf(stderr, "num_processors: %ld\n", num_processors);

    ibdqc_alloc(&qc, 300, 2, 10000, num_processors);
    
    // parse argument
    ibdqc_arg_parse(&qc, argc, argv);

    // start time:
    time(&t);
    fprintf(stderr, "Start time: %s", ctime(&t));


    // open files
    fp_chr_length = fopen(qc.args.file1->filename[0], "r");
    assert(fp_chr_length != NULL);
    fp_ibd = fopen(qc.args.file2->filename[0], "r");
    assert(fp_ibd != NULL);

    // read data
    ibdqc_read_chr_length(&qc, fp_chr_length);
    ibdqc_read_ibd(&qc, fp_ibd);

    fclose(fp_ibd);
    fclose(fp_chr_length);

    // multithreading
    // 	1. calculate the gamma const for each ibd segment
    // 	2. calculate estimated total ibd due to each generation
    tpool_run(&qc.tpool, qc.ibd_cM.size - 1, 0, 1000, ibdqc_job_func_calc_prob_consts);
    tpool_run(&qc.tpool, 100, 3, 1, ibdqc_job_func_calc_total_ibd_due_to_g);

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
    test_read_ibd_chrlen_from_file(argc, argv);

    return 0;
}
