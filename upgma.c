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
} clust_t;

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
clust_init(clust_t *self, const char *filename, size_t N)
{

    size_t num_pairs = N * (N - 1) / 2;

    float *cM = malloc(sizeof(*cM) * N * N);
    self->col_min_value = calloc(N, sizeof(*self->col_min_value));
    self->col_min_index = calloc(N, sizeof(*self->col_min_index));
    self->m_sid = calloc(2 * N - 1, sizeof(*self->m_sid));
    self->lnks = calloc(N - 1, sizeof(*self->lnks));
    self->leaf_cnt = malloc((2 * N - 1) * sizeof(*self->leaf_cnt));

    assert(cM && self->col_min_index && self->col_min_value && self->m_sid && self->lnks
           && self->leaf_cnt && "Can't alloc memory");

    /* read precomputed distance matrix */
    fprintf(stderr, "reading matrix ...\n");
    FILE *fp = fopen(filename, "rb");
    assert(fp != NULL && "can't open matrix file");
    size_t nread = fread(cM, sizeof(*cM), N * N, fp);
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

    /* intialization */
    self->N = N;
    self->n = self->N;
    self->cM = cM;

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
        min_row = 0;
        index = row_col_to_index(0, col);
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
    float *col_mins = self->col_min_value;
    size_t *col_min_row = self->col_min_index;
    size_t min_col;
    min = col_mins[0];
    min_col = 0;
    for (int col = 1; col < n; col++) {
        value = col_mins[col];
        if (value < min) {
            min = value;
            min_col = col;
        }
    }
    *prow = col_min_row[min_col];
    *pcol = min_col;
    if (*prow >= n || *prow < 0 || *pcol < 0 || *pcol >= n) {
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
    /* TODO: after debugging is done change stderr back to stdout*/
    fprintf(
        stderr, "%zu\t%zu\t%.4f\t%zu\n", lnk->id1, lnk->id2, lnk->height, lnk->leaves);
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
        cM[index2] = cM[index1] * weight1 + cM[index2] * weight2;
        cM[index2] /= (weight1 + weight2);
        /* indicating sample indexes in running matrix  */
    }

    /*update samples labels*/
    sidnew = ++(self->m_sid_last);
    sid[small] = sidnew;

    /* 2. swap row slots with last slots to condense the matrix */
    last = n - 1;
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
        /*update samples labels*/
        sid[large] = sid[last];
        sid[last] = -1;
        /* update col min info */
        self->col_min_value[large] = self->col_min_value[last];
        self->col_min_index[large] = self->col_min_index[last];
        /* update min_rows from last to row */
        for (size_t i = 0; i < n; i++)
            if (self->col_min_index[i] == last)
                self->col_min_index[i] = large;
    } else { /* when large = last */
        /* update min_rows from last to row */
        for (size_t i = 0; i < n; i++)
            if (self->col_min_index[i] == last)
                self->col_min_index[i] = small;
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
    float *cM = self->cM;
    size_t small = row > col ? col : row;

    /* 1. for those previous col min index is #col or #row, and
     * 2. for column #sidnew. */
    for (size_t i_col = 0; i_col < n; i_col++) {
        min_row = cmin_ind[i_col];
        if (min_row == col || min_row == row || i_col == small) {
            min_row = 0;
            index = row_col_to_index(0, i_col);
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

        /* 3. for those #sidnew row has a value smaller than column min  
        if (min_row != col && min_row != row && i_col != small)*/
       	else{
            index = row_col_to_index(small, i_col);
            value = cM[index];
            /*calcuated value < min*/
            if (value < cmin[i_col]) {
                cmin[i_col] = value;
                cmin_ind[i_col] = small;
            }
            /*calcuated value > min; no change*/
        }
    }

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
clust_main(const char *fn_matrix_binary, size_t N, int to_transform)
{
    clust_t myclust;
    clust_init(&myclust, fn_matrix_binary, N);
    if (to_transform != 0)
        clust_matrix_transform(&myclust);
    clust_col_min_init(&myclust);

    size_t i, j, index;
    while (myclust.n >= 2) {
        index = clust_find_min_cM_index(&myclust, &i, &j);
        fprintf(stderr, "run_num_samples: %lu, index: %lu, i: %lu, j: %lu\n", myclust.n,
            index, i, j);
        clust_update_links(&myclust, i, j);
        clust_update_mat_and_sid(&myclust, i, j);
        clust_update_col_min(&myclust, i, j);
    }
    clust_free(&myclust);

    return 0;
}

int
main(int argc, char *argv[])
{
    int to_transform = 0;
    if (argc < 3 || argc > 4) {
        fprintf(stderr, "\n\nUsage: ./upgma <float32_NxN_matrix_binary_file> "
                        "<N> [<-s]>  2>/dev/null \n"
                        "\twith `-s`: similarity matrix\n"
                        "\twithout `-s`: distance matrix\n\n");
        exit(1);
    }
    if (argc == 4 && strcmp(argv[3], "-s") == 0) {
        to_transform = 1;
    }

    clust_main(argv[1], atoi(argv[2]), to_transform);

    return 0;
}
