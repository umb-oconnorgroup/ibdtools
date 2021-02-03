/* Clustering algorithm: UPGMA.
 * Examples can be found in wikipedia UPGMA
 *
 * Two data structures are maintained in the process of clusterization
 * 	1. A tree building from merged nodes
 * 	2. A distance matrix of internal nodes and unmerged leaf nodes.
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

struct NODE {
    double height;
    double len_to_parent;
    int label;
    int parent;
    int left;
    int right;
};

struct CLUST {
    int num_samples;
    size_t sample_name_len;
    double *cM;
    char *samples;
    struct NODE *nodes;
    size_t num_pairs;
    int num_run_samples;
    int *run_labels;
    int *run_weights;
    int run_label_max;
    int *orders; /* order of sample index after clusterization */
    double *col_mins;
    int *col_min_row;
};

static inline int
upper_trangle_matrix_row_col_to_array_index(int id1, int id2, int num_samples)
{
    int temp;
    if (id1 > id2) {
        temp = id1;
        id1 = id2;
        id2 = temp;
    }
    if ((id1 != id2 && id1 >= 0 && id2 < num_samples) == 0)
        fprintf(stderr, "id1 %d, id2 %d\n", id1, id2);
    assert(id1 != id2 && id1 >= 0 && id2 < num_samples);
    return (2 * num_samples - id1 - 1) * id1 / 2 + id2 - id1 - 1;
}

/* convert to index to row and col num */
void
array_index_to_upper_trangle_matrix_row_col(
    int index, int num_samples, int *pid1, int *pid2)
{
    int row, col, row_len;
    for (row = 0; row < num_samples - 1; row++) {
        row_len = num_samples - 1 - row;
        if (index > row_len - 1)
            index -= row_len - 1; /*not in current row*/
        else {
            col = index + num_samples - row_len; /*find the number*/
            /*    -----   ---------------------  */
            /*     1's      0's                 */
            break;
        }
    }
    *pid1 = row;
    *pid2 = col;
}

void
clust_init_from_binary(struct CLUST *self, const char *filename)
{
    struct IBDTOTAL {
        double *cM;    // array of size equal: num_pairs = n*(n-1)/2
        char *samples; // size of n
        int num_smaple;
        size_t sample_name_len;
    };

    struct IBDTOTAL myibdtotal, *ibdtot = &myibdtotal;
    double max_cM = 0;

    // read ibdtotal structure for size information
    FILE *fp = fopen(filename, "rb");
    assert(fp != NULL);
    fread(ibdtot, sizeof(struct IBDTOTAL), 1, fp);
    ibdtot->cM = NULL;
    ibdtot->samples = NULL;

    // alloc according to size information
    int num_samples = ibdtot->num_smaple;
    int num_pairs = num_samples * (num_samples - 1) / 2;
    int num_chars = num_samples * ibdtot->sample_name_len;
    ibdtot->cM = (double *) malloc(sizeof(double) * num_pairs);
    assert(ibdtot->cM != NULL);
    ibdtot->samples = (char *) malloc(sizeof(char) * num_chars);
    assert(ibdtot->samples != NULL);
    self->col_mins = (double *) calloc(num_samples, sizeof(double));
    self->col_min_row = (int *) calloc(num_samples, sizeof(int));

    // read in arrays according to size info
    fread(ibdtot->cM, sizeof(double), num_pairs, fp);
    fread(ibdtot->samples, sizeof(char), num_chars, fp);
    fclose(fp);
    fp = NULL;

    /* Allocate memory */
    self->run_labels = (int *) calloc((2 * num_samples - 1), sizeof(int));
    self->run_weights = (int *) calloc(num_samples, sizeof(int));
    self->nodes = (struct NODE *) calloc((2 * num_samples - 1), sizeof(struct NODE));
    self->orders = (int *) calloc(num_samples, sizeof(int));

    /* intialization */
    self->num_samples = num_samples;
    self->sample_name_len = ibdtot->sample_name_len;
    self->num_pairs = num_pairs;
    self->num_run_samples = self->num_samples;
    self->samples = ibdtot->samples;
    ibdtot->samples = NULL;
    self->cM = ibdtot->cM;
    ibdtot->cM = NULL;

    for (int i = 0; i < num_samples; i++)
        self->run_labels[i] = i;
    for (int i = num_samples; i < 2 * num_samples - 1; i++)
        self->run_labels[i] = -1;
    for (int i = 0; i < num_samples; i++)
        self->run_weights[i] = 1;

    self->run_label_max = num_samples - 1;

    for (int i = 0; i < 2 * num_samples - 1; i++) {
        self->nodes[i].height = 0;
        self->nodes[i].left = -1;
        self->nodes[i].right = -1;
        self->nodes[i].parent = -1;
        self->nodes[i].label = i;
    }

    /* find max_cM */
    for (int i = 0; i < num_pairs; i++)
        if (self->cM[i] > max_cM)
            max_cM = self->cM[i];
    max_cM += 1;

    /* reverse -> similiary matrix -> dissimiliary matrix */
    for (int i = 0; i < num_pairs; i++)
        self->cM[i] = max_cM - self->cM[i];

    /* find col min and the col of the min*/
    int row_min = -1;
    double min = 0, value = 0, *cM = self->cM;
    double *col_mins = self->col_mins;
    int *col_min_rows = self->col_min_row;
    int index;
    for (int col = 0; col < num_samples; col++) {
        row_min = 0;
        min = 99999.0;
        for (int row = 0; row < num_samples; row++) {
            if (row == col)
                continue;
            index = upper_trangle_matrix_row_col_to_array_index(row, col, num_samples);
            value = cM[index];
            if (value < min) {
                min = value;
                row_min = row;
            }
        }
        col_mins[col] = min;
        col_min_rows[col] = row_min;
    }
}

void
clust_init(struct CLUST *self, const char *sample_file, const char *data_file)
{
    FILE *fp;
    // char buff[1000];
    char *p = NULL;
    size_t size = 1000;
    size_t counter = 0;
    size_t num_pairs = 0;
    size_t max_name_len = 0;
    size_t temp_len = 0;
    int ret = 0;

    fp = fopen(sample_file, "r");
    assert(fp != NULL);

    /* Read sample file 1st round,
     * get max sample name and count # of samples*/
    getline(&p, &size, fp);
    while (getline(&p, &size, fp) > 0 && *p != '\n') {
        counter++;
        temp_len = strlen(p);
        if (max_name_len < temp_len)
            max_name_len = temp_len;
    }
    self->num_samples = counter;
    self->num_run_samples = counter;
    self->sample_name_len = max_name_len + 1; /* add 1 for '\0' */

    /* Allocate memory */
    num_pairs = (counter - 1) * counter / 2;
    self->samples
        = (char *) calloc(self->num_samples * self->sample_name_len, sizeof(char));
    self->cM = (double *) calloc(num_pairs, sizeof(double));
    self->num_pairs = num_pairs;
    self->run_labels = (int *) calloc((2 * counter - 1), sizeof(int));
    self->run_weights = (int *) calloc(counter, sizeof(int));
    self->nodes = (struct NODE *) calloc((2 * counter - 1), sizeof(struct NODE));
    self->orders = (int *) calloc(counter, sizeof(int));
    self->col_mins = (double *) calloc(counter, sizeof(double));
    self->col_min_row = (int *) calloc(counter, sizeof(int));

    /* intialization */
    for (int i = 0; i < counter; i++)
        self->run_labels[i] = i;
    for (int i = counter; i < 2 * counter - 1; i++)
        self->run_labels[i] = -1;
    for (int i = 0; i < counter; i++)
        self->run_weights[i] = 1;
    self->run_label_max = counter - 1;

    for (int i = 0; i < 2 * counter - 1; i++) {
        self->nodes[i].height = 0;
        self->nodes[i].left = -1;
        self->nodes[i].right = -1;
        self->nodes[i].parent = -1;
        self->nodes[i].label = i;
    }

    // printf("counter=== %ld\n",counter);
    // printf("dsddsdd=== %d\n", self->nodes[0].left);

    /* Read sample file 2nd time
     * get sample names
     * */
    fseek(fp, 0, SEEK_SET);
    for (int i = 0; i < counter; i++) {
        p = self->samples + i * self->sample_name_len;
        size = self->sample_name_len;
        getline(&p, &size, fp);
        for (int j = 0; j < self->sample_name_len; j++) {
            if (p[j] == '\n') {
                p[j] = '\0';
                break;
            }
        }
    }
    fclose(fp);
    fp = NULL;

    /* Read data file */
    fp = fopen(data_file, "r");
    assert(fp != NULL);
    for (size_t i = 0; i < num_pairs; i++) {
        ret = fscanf(fp, "%lf", self->cM + i);
        assert(ret != EOF);
    }

    /* find col min and the col of the min*/
    int row_min = -1;
    int num_samples = self->num_samples;
    double min = 0, value = 0, *cM = self->cM;
    double *col_mins = self->col_mins;
    int *col_min_rows = self->col_min_row;
    int index;
    for (int col = 0; col < num_samples; col++) {
        row_min = 0;
        min = 99999.0;
        for (int row = 0; row < num_samples; row++) {
            if (row == col)
                continue;
            index = upper_trangle_matrix_row_col_to_array_index(row, col, num_samples);
            value = cM[index];
            if (value < min) {
                min = value;
                row_min = row;
            }
        }
        col_mins[col] = min;
        col_min_rows[col] = row_min;
    }
}

void
clust_print_info(struct CLUST *self)
{
    int index = 0;

    fprintf(stderr, "samples: ");
    for (size_t i = 0; i < self->num_samples; i++)
        fprintf(stderr, "%s, ", self->samples + i * self->sample_name_len);
    fprintf(stderr, "\n");

    fprintf(stderr, "run_labels: ");
    for (size_t i = 0; i < self->num_run_samples; i++)
        fprintf(stderr, "%d, ", self->run_labels[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "weights: ");
    for (size_t i = 0; i < self->num_run_samples; i++)
        fprintf(stderr, "%d, ", self->run_weights[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "col_min: ");
    for (size_t i = 0; i < self->num_run_samples; i++)
        fprintf(stderr, "%g, ", self->col_mins[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "col_min_row: ");
    for (size_t i = 0; i < self->num_run_samples; i++)
        fprintf(stderr, "%d, ", self->col_min_row[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "matrix: \n");
    for (size_t i = 0; i < self->num_samples; i++) {
        for (size_t j = 0; j < self->num_samples; j++) {
            if (i == j)
                fprintf(stderr, "%.1f\t", 0.0);
            else {
                index = upper_trangle_matrix_row_col_to_array_index(
                    i, j, self->num_samples);
                fprintf(stderr, "%.1lf\t", self->cM[index]);
            }
        }
        fprintf(stderr, "\n");
    }
}

int
clust_find_min_cM_index(struct CLUST *self, int *prow, int *pcol)
{
    /* find the min of the matrix */
    double min = 99999, value;
    int num_run_samples = self->num_run_samples;
    double *col_mins = self->col_mins;
    int *col_min_row = self->col_min_row;
    int min_index;
    for (int col = 0; col < num_run_samples; col++) {
        value = col_mins[col];
        if (value < min) {
            min = value;
            min_index = col;
        }
    }
    *prow = col_min_row[min_index];
    assert(*prow < num_run_samples && *prow >= 0);
    *pcol = min_index;
    assert(*pcol < num_run_samples && *pcol >= 0);
    assert(*pcol != *prow);

    return upper_trangle_matrix_row_col_to_array_index(*prow, *pcol, num_run_samples);
}

void
clust_update_nodes_before_clustering(struct CLUST *self, int row, int col)
{
    struct NODE *parent, *child1, *child2;
    int index; /* index to element of matrix */
    int temp;
    double cM;
    if (row > col) {
        temp = row;
        row = col;
        col = temp;
    }
    for (int i = 0; i < self->num_run_samples; i++)
        assert(self->run_labels[i] >= 0);
    parent = self->nodes + self->run_label_max + 1;
    child1 = self->nodes + self->run_labels[row];
    child2 = self->nodes + self->run_labels[col];

    index = upper_trangle_matrix_row_col_to_array_index(row, col, self->num_samples);

    cM = self->cM[index];
    parent->height = cM / 2;
    parent->left = child1 - self->nodes;
    parent->right = child2 - self->nodes;
    parent->label = parent - self->nodes;
    child1->parent = parent - self->nodes;
    child2->parent = parent - self->nodes;
    child1->len_to_parent = parent->height - child1->height;
    child2->len_to_parent = parent->height - child2->height;

    fprintf(stderr, "child1: %d, child2: %d -> parent: %d \n", parent->left,
        parent->right, child1->parent);
    fprintf(stderr, "child1- height: %lf, len_to_parent: %lf, left: %d, label: %d\n",
        child1->height, child1->len_to_parent, child1->left, child1->label);
    fprintf(stderr, "child2- height: %lf, len_to_parent: %lf\n", child2->height,
        child2->len_to_parent);

    assert(child2 >= self->nodes);
    assert(child1 >= self->nodes);
}

void
clust_update_matrix_after_clustering(struct CLUST *self, int row, int col)
{
    int temp;
    int index1, index2;
    int weight1, weight2;
    int num_samples = self->num_samples;
    int num_run_samples = self->num_run_samples;
    int last; // last effective sample
    double *cM = self->cM;
    if (row > col) {
        temp = row;
        row = col;
        col = temp;
    }
    weight1 = self->run_weights[row];
    weight2 = self->run_weights[col];
    for (size_t i = 0; i < num_run_samples; i++) {
        if (i == row || i == col)
            continue;
        /* merging sample row and sample col */
        index1 = upper_trangle_matrix_row_col_to_array_index(i, row, num_samples);
        index2 = upper_trangle_matrix_row_col_to_array_index(i, col, num_samples);
        cM[index1] = cM[index1] * weight1 + cM[index2] * weight2;
        cM[index1] /= (weight1 + weight2);
        /* indicating sample indexes in running matrix  */
    }
    /*update samples labels*/
    self->run_label_max++;
    self->run_labels[row] = self->run_label_max;
    /*update weights */
    self->run_weights[row] += self->run_weights[col];

    last = num_run_samples - 1;
    if (col != last) {
        for (size_t i = 0; i < num_run_samples; i++) {
            /* copy last running sample to sample col*/
            if (i == last || i == col)
                continue;
            index1 = upper_trangle_matrix_row_col_to_array_index(i, col, num_samples);
            index2 = upper_trangle_matrix_row_col_to_array_index(i, last, num_samples);
            cM[index1] = cM[index2];
            cM[index2] = -1;
        }
        /*update samples labels*/
        if (col != last) {
            self->run_labels[col] = self->run_labels[last];
            self->run_labels[last] = -1;
        }
        /*update weights */
        self->run_weights[col] = self->run_weights[last];
        self->run_weights[last] = 0;
    }

    /*update num of running samples*/
    self->num_run_samples -= 1;
    num_run_samples = self->num_run_samples;

    /* update col min info */
    int index = 0;
    self->col_mins[col] = self->col_mins[last];
    self->col_min_row[col] = self->col_min_row[last];

    // update min_row from last to col
    for (int i = 0; i < num_run_samples; i++)
        if (self->col_min_row[i] == last && last != col)
            self->col_min_row[i] = col;

    // update min for those
    // 	1. the caculated col: #row
    // 	2. whose previous min in col or row
    // 	3. whose min not in col or row but caculate col is smaller than current
    // min

    double *col_mins = self->col_mins;
    int *col_min_row = self->col_min_row;
    int min_row;
    double min, value;

    for (int i_col = 0; i_col < num_run_samples; i_col++) {
        /* min values are in rows that is merged  or i is in col #row
         * find min for each columns
         */
        index = self->col_min_row[i_col];
        if (index == col || index == row || i_col == row) {
            min_row = -1;
            min = 99999.0;
            for (int j_row = 0; j_row < num_run_samples; j_row++) {

                if (j_row == i_col)
                    continue;

                int ind = upper_trangle_matrix_row_col_to_array_index(
                    j_row, i_col, num_samples);
                value = cM[ind];
                if (value < min) {
                    min = value;
                    min_row = j_row;
                }
            }
            col_mins[i_col] = min;
            col_min_row[i_col] = min_row;
        }
    }

    for (int i_col = 0; i_col < num_run_samples; i_col++) {
        index = self->col_min_row[i_col];
        if (index != col && index != row && i_col != row) {
            if (i_col == row)
                continue;
            int ind
                = upper_trangle_matrix_row_col_to_array_index(row, i_col, num_samples);

            /*calcuated value < min*/
            if (self->cM[ind] < self->col_mins[i_col]) {
                self->col_mins[i_col] = self->cM[ind];
                self->col_min_row[i_col] = row;
            }
            /*calcuated value > min; no change*/
        }
    }

    for (int i_col = 0; i_col < num_run_samples && num_run_samples >= 3; i_col++)
        assert(col_min_row[i_col] >= 0 && col_min_row[i_col] < num_run_samples);
}

/*
 * First, generate a list of index of traversing orders with duplicated indices
 * 	if the path runs thru the node twice
 * Second, onece the list is made, based on the type of jump, print differnt
 * 	string. See comments below.
 * Third, get sample order after clusterzation.
 */
void
clust_print_newick_tree(struct CLUST *self)
{
    int *stack = (int *) calloc((2 * self->num_samples), sizeof(int));
    int last = -1;
    int *list = (int *) calloc(4 * self->num_samples, sizeof(int));
    int llast = -1;
    struct NODE *nodes = self->nodes, *pNode = NULL, *pPrevNode = NULL;

    /* add root node to stack */
    last++;
    stack[last] = self->run_label_max; /*root node index */

    /* access all node: left first */
    while (last > -1) {
        /* pop out */
        pNode = nodes + stack[last];
        last--;

        if (pPrevNode != NULL) {
            do /* go up */
            {
                pPrevNode = nodes + pPrevNode->parent;
                llast++;
                list[llast] = pPrevNode->label;
            } while (pPrevNode->label != pNode->parent);
        }

        /* go down right */
        llast++;
        list[llast] = pNode->label;

        while (pNode->left != -1) {
            /*push back the right node*/
            last++;
            stack[last] = pNode->right;

            /* go down left */
            pNode = nodes + pNode->left;
            llast++;
            list[llast] = pNode->label;
        }
        pPrevNode = pNode;
    }

    while (pNode != nodes + self->run_label_max) {
        pNode = nodes + pNode->parent;
        llast++;
        list[llast] = pNode->label;
    }

    for (int i = 0; i <= llast; i++)
        fprintf(stderr, "%d ", list[i]);
    fprintf(stderr, "\n");

    /*
     * 	parent-left : (
     * 	child-parent:  leave: label:edge
     * 		      non leave: ):edge
     * 	parent-right: ,
     * */

    for (int i = 1; i <= llast; i++) {
        pNode = nodes + list[i];
        pPrevNode = nodes + list[i - 1];
        if (pNode->label == pPrevNode->left)
            fprintf(stdout, "(");
        else if (pNode->label == pPrevNode->right)
            fprintf(stdout, ",");
        else {
            if (pPrevNode->left == -1)
                fprintf(stdout, "%d:%g", pPrevNode->label, pPrevNode->len_to_parent);
            else
                fprintf(stdout, "):%g", pPrevNode->len_to_parent);
        }
        /*
        if (i == llast)
          fprintf(stdout, "%d:%g);\n", pNode->label, pNode->len_to_parent);
          */
    }
    fprintf(stdout, ");\n");

    /* save new sample order */
    int num_samples = self->num_samples;
    int *orders = self->orders;
    for (int i = 0, j = 0; i <= llast && j < num_samples; i++) {
        if (list[i] >= num_samples)
            continue;
        orders[j] = list[i];
        j++;
    }

    for (int i = 0; i < num_samples; i++)
        fprintf(stdout, "%d ", orders[i]);
    fprintf(stdout, "\n");

    free(stack);
    free(list);
    stack = NULL;
    list = NULL;
    pNode = NULL;
}

int
test_wiki()
{
    struct CLUST myclust;
    clust_init(&myclust, "samples.txt", "data.txt");
    clust_print_info(&myclust);

    int i, j, index;
    while (myclust.num_run_samples >= 2) {
        index = clust_find_min_cM_index(&myclust, &i, &j);
        fprintf(stderr, "\n\n----------- after find min\n");
        clust_print_info(&myclust);
        fprintf(stderr, "%d, %d, %d\n", index, i, j);
        clust_update_nodes_before_clustering(&myclust, i, j);
        clust_update_matrix_after_clustering(&myclust, i, j);
        fprintf(stderr, "\n\n----------- after update matrix\n");
        clust_print_info(&myclust);
    }
    fprintf(stderr, "%d\n", myclust.run_label_max);

    clust_print_newick_tree(&myclust);
    return 0;
}

int
test_binary(const char *fn_matrix_binary)
{
    struct CLUST myclust;
    clust_init_from_binary(&myclust, fn_matrix_binary);

    int i, j, index;
    while (myclust.num_run_samples >= 2) {
        index = clust_find_min_cM_index(&myclust, &i, &j);
        fprintf(stderr, "run_num_samples: %d, index: %d, i: %d, j: %d\n",
            myclust.num_run_samples, index, i, j);
        clust_update_nodes_before_clustering(&myclust, i, j);
        clust_update_matrix_after_clustering(&myclust, i, j);
    }
    // fprintf(stderr, "%d\n", myclust.run_label_max);

    clust_print_newick_tree(&myclust);
    return 0;
}

int
main(int argc, char *argv[])
{
    // test_wiki();
    if (argc != 2) {
        fprintf(stderr, "\n\nUsage: ./ibdclust <binary_matrix_file> 2>/dev/null \n\n"
                        "stdout will output newick tree(1st line) and sample order "
                        "after clusterization (2nd line)\n\n ");
        exit(1);
    }

    test_binary(argv[1]);
    return 0;
}
