#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

typedef struct {
    unsigned int *val;
    size_t size;
} ne_t;

typedef struct {
    ne_t N;
    unsigned int g_star;
    unsigned int G;
} ibdqc_t;

void
ne_alloc(ne_t *self, size_t size)
{
    assert(size > 0);
    self->size = size;
    self->val = malloc(sizeof(*self->val) * size);
}

void
ne_free(ne_t *self)
{
    if (self && self->val) {
        free(self->val);
        self->val = NULL;
    }
}

void
ne_init(ne_t *self, unsigned int refsize)
{
    for (int i = 0; i < self->size; i++) {
        self->val[i] = refsize;
    }
}

void
ibdqc_alloc(ibdqc_t *self, unsigned int G, unsigned int g_star, unsigned int refsize)
{
    self->G = G;
    self->g_star = g_star;
    ne_alloc(&self->N, G + 1);
    ne_init(&self->N, refsize);
}

void
ibdqc_free(ibdqc_t *self)
{
    ne_free(&self->N);
}

double
ibdqc_calc_prob_ibd_due_to_g(ibdqc_t *self, double l, unsigned int g)
{
    /* variables */
    double sum = 0, prod = 0, tmp = 0;
    double const_gamma = 0;
    double probability = 0;

    /* parameters */
    unsigned int G = self->G;
    unsigned int g_star = self->g_star;
    unsigned int *N = self->N.val;
    assert(G == self->N.size - 1);
    double alpha = l / 50.0;
    double beta = 1 - 0.5 / N[G];

    /* calc const_gamma */
    sum = 0;

    for (int g1 = g_star; g1 <= G; g1++) {
        prod = 1;
        prod *= g1 * g1 / 2500.0;
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
    prod *= G * G / tmp + (2 * G - 1) / tmp / tmp + 2 / tmp / tmp / tmp;
    sum += prod;

    const_gamma = sum;

    /* calc probability due to g given ibd length */
    prod = 1;
    prod *= 1 / const_gamma;
    prod *= g * g / 2500.0 * exp(-l * g / 50);
    for (int g1 = g_star; g1 < g; g1++) {
        prod *= 1 - 0.5 / N[g1];
    }
    prod *= 0.5 / N[g];

    probability = prod;

    return probability;
}

void
test_ibdqc_calc_prob_ibd_due_to_g()
{
    ibdqc_t qc;
    ibdqc_alloc(&qc, 300, 2, 10000);
    /* make ibd length array */
    size_t size = 10000;
    double *l_arr = malloc(sizeof(*l_arr) * size);
    for (int i = 0; i < size; i++) {
        l_arr[i] = rand() % 10000 / 100;
    }

    /* calculate ibd due to generation g */
    for (int g = 3; g < 100; g++) {
        double sum = 0;
        for (int i = 0; i < size; i++) {
            double l = l_arr[i];
            sum += l * ibdqc_calc_prob_ibd_due_to_g(&qc, l, g);
        }
        printf("g = %d, total_ibd = %g\n", g, sum);
    }
    ibdqc_free(&qc);
}

int
main()
{
    test_ibdqc_calc_prob_ibd_due_to_g();
    return 0;
}
