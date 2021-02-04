/* Clustering algorithm: UPGMA.

 * Examples can be found in wikipedia UPGMA
 *
 * Data structures maintained in the process of clusterization
 * 	1. linkage matrix similiar to scipy.cluster.hierachy.linkage example
 * 	2. A distance matrix of internal nodes and unmerged leaf nodes.
 * 	3. Arrays of col min: maintaing the col min info is helpful for performance
 * 	4. Arrays of col min index
 * 	5. Arrays of number of leaves a node spans
 *
 * Notes:
 *   	1. To save memory, the distance matrix is presented as as triangular matrix.
 * 	2.  When two nodes is merged, the distance of the merged node to the rest
 * 	is calculated. And the matrix shrinked.
 * */

#define _GNU_SOURCE
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct {
    size_t id1;
    size_t id2;
    float height;
    size_t leaves;
} lnk_t;

typedef struct {
    size_t N;  /* total number of samples */
    size_t n;  /* running number of unmerged samples */
    float *cM; /* the array underlying the matrix */
    lnk_t *lnks;
    size_t nlink;
    size_t *m_sid;    /* sample ids for matrix columns */
    size_t *leaf_cnt; /* how many leveas it covers */
    size_t m_sid_last;
    float *col_min_value;
    size_t *col_min_index;
    int *merge_col_updated; /* tells if cells in merge col are updated. Used for
                               optimization */
    int verbose;
} clust_t;

typedef struct {
    int verbose;
    int to_transform;
    char *matrix_file_name;
    size_t N;
    char *usage;
    int cnt_positional_params;
} args_t;

static inline size_t
row_col_to_index(size_t row, size_t col)
{
    size_t temp;
    if (row < col) {
        temp = row;
        row = col;
        col = temp;
    }
    return row * (row - 1) / 2 + col;
}

static inline void
index_to_row_col(size_t x, size_t *prow, size_t *pcol)
{
    size_t row;
    for (row = 1; x >= row; x -= row++)
        ;
    *prow = row;
    *pcol = x;
}

void
clust_init(clust_t *self, float * cM, size_t N)
{
    self->col_min_value = calloc(N, sizeof(*self->col_min_value));
    self->col_min_index = calloc(N, sizeof(*self->col_min_index));
    self->m_sid = calloc(2 * N - 1, sizeof(*self->m_sid));
    self->lnks = calloc(N - 1, sizeof(*self->lnks));
    self->leaf_cnt = calloc((2 * N - 1), sizeof(*self->leaf_cnt));
    self->merge_col_updated = calloc(N, sizeof(*self->merge_col_updated));

    assert(self->col_min_index && self->col_min_value && self->m_sid && self->lnks
           && self->leaf_cnt && self->merge_col_updated && "Can't alloc memory");


    /* intialization */
    self->N = N;
    self->n = self->N;
    self->cM = cM;
    self->verbose = 0;

    for (size_t i = 0; i < N; i++) {
        self->leaf_cnt[i] = 1;
        self->m_sid[i] = i;
    }

    self->m_sid_last = N - 1;
    self->nlink = 0;

    fprintf(stderr, "done initiating ...\n");
    cM = NULL;
}

void
clust_free(clust_t *self)
{
    free(self->cM);
    self->cM = NULL;
    free(self->col_min_value);
    self->col_min_value = NULL;
    free(self->col_min_index);
    self->col_min_index = NULL;
    free(self->m_sid);
    self->m_sid = NULL;
    free(self->lnks);
    self->lnks = NULL;
    free(self->leaf_cnt);
    self->leaf_cnt = NULL;
    free(self->merge_col_updated);
    self->merge_col_updated = NULL;
}

void
clust_col_min_init(clust_t *self)
{
    fprintf(stderr, "finding min info ...\n");
    /* find col min and the col of the min*/
    size_t min_row;
    float min = 0, value = 0;
    float *cmin = self->col_min_value;
    size_t *cmin_ind = self->col_min_index;
    float *cM = self->cM;

    size_t index;
    for (size_t col = 0; col < self->N; col++) {
        min_row = (col == 0) ? 1 : 0;
        index = row_col_to_index(min_row, col);
        min = cM[index];
        for (size_t row = 0; row < self->N; row++) {
            if (row == col)
                continue;
            index = row_col_to_index(row, col);
            value = cM[index];
            if (value < min) {
                min = value;
                min_row = row;
            }
        }
        cmin[col] = min;
        cmin_ind[col] = min_row;
    }
}

static inline size_t
clust_find_min_cM_index(clust_t *self, size_t *prow, size_t *pcol)
{
    /* find the min of the matrix */
    float min, value;
    size_t n = self->n;
    float *cmin = self->col_min_value;
    size_t *cmin_ind = self->col_min_index;
    size_t min_col;
    min_col = 0;
    min = cmin[min_col];
    for (size_t col = 1; col < n; col++) {
        value = cmin[col];
        if (value < min) {
            min = value;
            min_col = col;
        }
    }
    if (self->verbose) {
        for (size_t col = 0; col < n; col++) {
            fprintf(stderr,
                "%c, col=%zu, cmin_ind[col]= %zu:", 'a' + (char) self->m_sid[col], col,
                cmin_ind[col]);
            for (size_t row = 0; row < n; row++) {
                if (row == col)
                    fprintf(stderr, "-,");
                else {
                    fprintf(stderr, "%0.2f,", self->cM[row_col_to_index(row, col)]);
                }
            }
            fprintf(stderr, "\n");
        }
    }
    *prow = cmin_ind[min_col];
    *pcol = min_col;
    if (*prow >= n || *prow < 0 || *pcol < 0 || *pcol >= n || *prow == *pcol) {
        fprintf(stderr, "row/col out of range: ");
        fprintf(stderr, "row=%zu, col=%zu, n=%zu ", *prow, *pcol, n);
        assert(*prow < n && *prow >= 0 && "row num of min element out of range");
        assert(*pcol < n && *pcol >= 0 && "col num of min element out of range");
        assert(*pcol != *prow && "min element is on the diagonal");
    }

    return row_col_to_index(*prow, *pcol);
}

void
clust_update_links(clust_t *self, size_t row, size_t col)
{
    size_t index; /* arry index of element in matrix */
    size_t temp;
    float value;
    size_t sid1, sid2, sid_new, *lcnt;
    sid_new = self->m_sid_last;
    lcnt = self->leaf_cnt;
    if (row < col) {
        temp = row;
        row = col;
        col = temp;
    }
    index = row_col_to_index(row, col);
    value = self->cM[index];

    sid_new += 1;
    sid1 = self->m_sid[row];
    sid2 = self->m_sid[col];
    lcnt[sid_new] = lcnt[sid1] + lcnt[sid2];

    lnk_t *lnk = self->lnks + self->nlink;
    lnk->id1 = sid1;
    lnk->id2 = sid2;
    lnk->height = value / 2;
    lnk->leaves = lcnt[sid_new];
    self->nlink += 1;

    if (self->verbose) {
        fprintf(stderr, "%zu (%zu)\t%zu (%zu) \t%.4f\t%zu\n", lnk->id1, lcnt[sid1],
            lnk->id2, lcnt[sid2], lnk->height, lnk->leaves);
    } else {
        fprintf(stdout, "%zu\t%zu\t%.4f\t%zu\n", lnk->id1, lnk->id2, lnk->height,
            lnk->leaves);
    }
}

void
clust_update_mat_and_sid(clust_t *self, size_t row, size_t col)
{
    size_t index1, index2;
    size_t weight1, weight2;
    size_t n = self->n;
    size_t last;
    float *cM = self->cM;
    size_t *lcnt = self->leaf_cnt;
    size_t *sid = self->m_sid;
    size_t sidnew;
    size_t large, small;

    if (row < col) {
        large = col;
        small = row;
    } else {
        large = row;
        small = col;
    }

    /* 1. merge sid[row] and sid[col] and put sid new in col slot */
    weight1 = lcnt[sid[large]];
    weight2 = lcnt[sid[small]];
    for (size_t i = 0; i < n; i++) {
        if (i == large || i == small)
            continue;
        /* merging sample row and sample col */
        index1 = row_col_to_index(i, large);
        index2 = row_col_to_index(i, small);
        if (cM[index1] == cM[index2])
            self->merge_col_updated[i] = 0;
        else {
            self->merge_col_updated[i] = 1;
            cM[index2] = cM[index1] * weight1 + cM[index2] * weight2;
            cM[index2] /= (weight1 + weight2);
        }
        /* indicating sample indexes in running matrix  */
    }

    /*update samples labels*/
    ++(self->m_sid_last);
    sidnew = (self->m_sid_last);
    sid[small] = sidnew;
    last = n - 1;
    sid[large] = sid[last];

    /* 2. swap row slots with last slots to condense the matrix */
    if (large != last) {
        for (size_t i = 0; i < n; i++) {
            /* copy last running sample to sample row*/
            if (i == last || i == large)
                continue;
            index1 = row_col_to_index(i, large);
            index2 = row_col_to_index(i, last);
            cM[index1] = cM[index2];
            cM[index2] = -1;
        }
    }

    /* 3. update num of running samples*/
    self->n -= 1;
}

/* Update min for those:
 * 	1. whose previous col min index is #col or #row
 * 	2. the caculated column of sidnew. i.e #row
 * 	3. whose #sidnew row has a value smaller than column min */
void
clust_update_col_min(clust_t *self, size_t row, size_t col)
{
    float *cmin = self->col_min_value;
    size_t *cmin_ind = self->col_min_index;
    size_t min_row;
    float min, value;
    size_t index = 0;
    size_t n = self->n;
    size_t last = n;
    float *cM = self->cM;
    size_t small, large;
    int need_recalc = 0;
    if (row < col) {
        large = col;
        small = row;
    } else {
        large = row;
        small = col;
    }

    cmin[large] = cmin[last];
    cmin_ind[large] = cmin_ind[last];

    /* 1. for those previous col min index is #col or #row, and
     * 2. for column #sidnew. */
    for (size_t i_col = 0; i_col < n; i_col++) {
        min_row = cmin_ind[i_col];
        need_recalc = 0;

        if (i_col == small) {
            need_recalc = 1;
            if (cM[row_col_to_index(small, large)]
                == self->lnks[self->nlink - 1].height) {
                cmin_ind[small] = large;
                /* if min_row == small, keep the same; */
                need_recalc = 0;
            }

        } else if (min_row == large || min_row == small) {
            need_recalc = 1;
            /* when the two merged cells have the same value, the merged value is still
             * min */
            if (self->merge_col_updated[i_col] == 0) {
                if (min_row == large)
                    cmin_ind[i_col] = small;
                /* if min_row == small, keep the same; */
                need_recalc = 0;
            }
        } else /* if (min_row != large && min_row != small) */ {
            /* no change only because merge of two non-min is still non-min
             * rediect from last to large*/
            if (min_row == last) {
                if (large != last) {
                    cmin_ind[i_col] = large;
                    cmin[i_col] = cmin[large];
                } else {
                    /* will be reculculated */
                }
            }
        }

        if (need_recalc == 0)
            continue;

        if (self->verbose)
            fprintf(stderr, "+");

        min_row = (i_col == 0) ? 1 : 0;
        index = row_col_to_index(min_row, i_col);
        min = cM[index];

        for (size_t j_row = 1; j_row < n; j_row++) {

            if (j_row == i_col)
                continue;

            index = row_col_to_index(j_row, i_col);
            value = cM[index];
            if (value < min) {
                min = value;
                min_row = j_row;
            }
        }
        cmin[i_col] = min;
        cmin_ind[i_col] = min_row;
    }
    if (self->verbose)
        fprintf(stderr, "\n");

    for (size_t i_col = 0; i_col < n && n >= 3; i_col++)
        if (cmin_ind[i_col] < 0 && cmin_ind[i_col] >= n) {
            fprintf(stderr, "col min index outof range: ");
            fprintf(stderr, "row=%zu, col=%zu, i_col=%zu, n=%zu, cmin_ind[i_col]=%zu\n",
                row, col, i_col, n, cmin_ind[i_col]);
            assert(cmin_ind[i_col] < 0 && cmin_ind[i_col] >= n);
        }
}

void
clust_matrix_transform(clust_t *self)
{

    fprintf(stderr, "transforming matrix ...\n");
    size_t nparis = self->N * (self->N - 1) / 2;
    float max, min, value;
    max = self->cM[0];
    min = max;
    for (size_t i = 1; i < nparis; i++) {
        value = self->cM[i];
        if (min > value)
            min = value;
        if (value > max)
            max = value;
    }
    for (size_t i = 0; i < nparis; i++) {
        self->cM[i] = max / (self->cM[i] - min + 1e-9);
    }
    fprintf(stderr, "\tmin=%g, max=%g ===>> min=%g, max=%g\n", min, max,
        max / (max - min + 1e-9), max / (min - min + 1e-9));
}

int
clust_main(clust_t *self, float *cM, size_t N, int to_transform, int verbose)
{
    clust_init(self, cM, N);
    if (to_transform != 0)
        clust_matrix_transform(self);
    clust_col_min_init(self);
    self->verbose = verbose;

    size_t i, j, index;
    while (self->n >= 2) {
        index = clust_find_min_cM_index(self, &i, &j);
        fprintf(stderr, "run_num_samples: %lu, index: %lu, i: %lu, j: %lu\n", self->n,
            index, i, j);
        assert(i != j && "col row should not be equal");
        clust_update_links(self, i, j);
        clust_update_mat_and_sid(self, i, j);
        clust_update_col_min(self, i, j);
    }
    clust_free(self);

    return 0;
}

int
main(int argc, char *argv[])
{
    args_t args;
    size_t N, num_pairs, nread;
    float *cM;
    FILE * fp;

    args.cnt_positional_params = 0;
    args.matrix_file_name = NULL;
    args.to_transform = 0;
    args.verbose = 0;
    args.N = 0;
    args.usage = "Usage: ./upgma [options]  matrix_file   num_sample\n"
                 "\t-s: transform similarity matrix to distance matrix\n"
                 "\t-v: toggle debuging info\n"
                 "\t-h/--help: help information\n"
		 "stderr: log info\n"
		 "stdout: linkage matrix\n";

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-s") == 0) {
            args.to_transform = 1;
        } else if (strcmp(argv[i], "-v") == 0) {
            args.verbose = 1;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            fprintf(stderr, "%s\n", args.usage);
            exit(1);
        } else {
            if (args.cnt_positional_params == 0)
                args.matrix_file_name = argv[i];
            else if (args.cnt_positional_params == 1)
                args.N = strtol(argv[i], NULL, 10);
            else {
                fprintf(stderr, "Error: to much positional parameters!\n");
                exit(1);
            }
            args.cnt_positional_params++;
        }
    }
    if (args.cnt_positional_params != 2) {
        fprintf(stderr, "Error: need exact 2 positional parameters!\n");
        exit(1);
    }

    N = args.N;
    num_pairs = N * (N-1) / 2;
    cM = calloc(N * N, sizeof(*cM));
    assert(cM && "can't alloc memory");

    /* read precomputed distance matrix */
    fprintf(stderr, "reading matrix ...\n");
    fp = fopen(args.matrix_file_name, "rb");
    assert(fp != NULL && "can't open matrix file");
    nread = fread(cM, sizeof(*cM), N * N, fp);
    assert(nread == N * N && "size of matrix file not consistent with sample numbers");
    fclose(fp);
    fp = NULL;

    /* matrix representation to array representation */
    fprintf(stderr, "make matrix triangle ...\n");
    for (size_t i = 1; i < N; i++)
        for (size_t j = 0; j < i; j++)
            cM[row_col_to_index(i, j)] = cM[i * N + j];
    cM = realloc(cM, num_pairs * sizeof(*cM));
    assert(cM && "failed to shrink memory of the matrix");

    clust_t myclust;
    clust_main(&myclust, cM, N, args.to_transform, args.verbose);

    return 0;
}
