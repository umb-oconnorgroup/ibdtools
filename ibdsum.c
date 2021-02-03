#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>

#define genomeCM 3545.83

typedef struct {
    char *name;
    size_t size;
    size_t nmemb;
} sample_t;

int
str_cmp(const void *s1, const void *s2)
{
    return strcmp((char *) s1, (char *) s2);
}

void
sample_read(sample_t *self, char *sample_list_fn)
{
    char buff[1000], *p = buff, *tok = NULL;
    size_t buff_size = 1000;
    size_t ncap = 0, nmemb = 0, inc = 10000;
    size_t size = 256; /*max_sample lenth */

    FILE *fp = fopen(sample_list_fn, "r");
    assert(fp != NULL && "cannot open sample file");

    ncap += inc;
    self->name = malloc(size * sizeof(*self->name) * ncap);
    assert(self->name != NULL && "cannot alloc memory");

    while (getline(&p, &buff_size, fp) > 0 && *p != '\n') {
        if (nmemb + 1 >= ncap) {
            // alloc memeory
            ncap += inc;
            self->name = realloc(self->name, size * sizeof(*self->name) * ncap);
            assert(self->name != NULL && "cannot alloc memory");
        }
        tok = strtok(p, "\n");
        assert(tok != NULL && "sample list format error");
        strcpy(self->name + size * nmemb, p);
        nmemb++;
    }
    ncap = nmemb;
    self->name = realloc(self->name, size * sizeof(*self->name) * ncap);
    assert(self->name != NULL && "cannot alloc memory");
    self->size = size;
    self->nmemb = nmemb;

    fclose(fp);
    fp = NULL;

    qsort(self->name, self->nmemb, self->size, str_cmp);
    fprintf(stderr, "ibdsum done reading sample list\n");
}

typedef struct {
    sample_t sample;
    float *cM;
    float *cM_temp; /* for addup, and avoid realloc when multiple binary files are to be
                       added */
    size_t nsam;
} ibdsum_t;

void
ibdsum_init(ibdsum_t *self, char *sample_fn)
{
    size_t ntotal = 0;
    sample_read(&self->sample, sample_fn);
    self->nsam = self->sample.nmemb;
    fprintf(stderr, "ibdsum start init\n");

    ENTRY e, *ep;
    /* create hashtable and add sample names to the table*/
    hcreate(self->nsam * 2);
    for (size_t i = 0; i < self->nsam; i++) {
        e.key = self->sample.name + i * self->sample.size;
        e.data = (void *) i;
        ep = hsearch(e, ENTER);
        assert(ep != NULL && "hsearch enter failure!");
    }
    fprintf(stderr, "ibdsum done add sample list to hash table\n");

    /* use nsam * nam lower triagular matrix to save cM info*/
    ntotal = self->nsam * (self->nsam - 1) / 2;
    self->cM = (float *) malloc(sizeof(*self->cM) * ntotal);
    self->cM_temp = (float *) malloc(sizeof(*self->cM_temp) * ntotal);
    assert(self->cM != NULL && "can not alloc for cM matrix");
    assert(self->cM_temp != NULL && "can not alloc for cM_temp matrix");

    for (size_t i = 0; i < ntotal; i++) {
        self->cM[i] = 0;
    }
    fprintf(stderr, "ibdsum done alloc matrix of total ibd\n");
}

void
ibdsum_free(ibdsum_t *self)
{
    hdestroy();
    free(self->cM);
    free(self->sample.name);
}

void
ibdsum_summarize_sorted_ibd(ibdsum_t *self)
{
    char buff[1000], *p = buff, *tok = NULL;
    size_t buff_size = 1000;
    size_t size = 256; /*max_sample lenth */
    char samplePair[1000];
    size_t cM_index = 0;
    size_t row = 0, col = 0;
    float cM_len = 0;
    ENTRY e, *ep;
    char *s1 = NULL, *s2 = NULL;
    samplePair[0] = '\0';
    size_t count = 0;

    while (getline(&p, &buff_size, stdin) > 0 && *p != '\n') {
        /* sample pair name */
        tok = strtok(p, "\t\n");
        assert(tok != NULL && "ibd format error");
        s1 = tok;
        /* start */
        tok = strtok(NULL, "\t\n");
        assert(tok != NULL && "ibd format error");
        /* end */
        tok = strtok(NULL, "\t\n");
        assert(tok != NULL && "ibd format error");
        /* cM lenght */
        tok = strtok(NULL, "\t\n");
        assert(tok != NULL && "ibd format error");
        cM_len = strtod(tok, NULL);

        /* if same sample pair, use the same position index
         * otherwise, update the position using sample names */
        if (strcmp(s1, samplePair) != 0) {
            /* split the sample pair name into sample names */
            strtok(s1, ":");
            s2 = strtok(NULL, ":");

            e.key = s1;
            /*fprintf(stderr, "s1: %s\n", s1);*/
            ep = hsearch(e, FIND);
            assert(ep != NULL && "can not find row from hash table!");
            row = (size_t) ep->data;

            e.key = s2;
            /*fprintf(stderr, "s2: %s\n", s2);*/
            ep = hsearch(e, FIND);
            col = (size_t) ep->data;
            assert(ep != NULL && "can not find col from hash table!");

            /* lower triangle, row should be > col */
            if (row < col) {
                /*swap*/
                cM_index = col;
                col = row;
                row = cM_index;
            }
            /* TODO verify the calculation */
            cM_index = row * (row - 1) / 2 + col;
        }
        self->cM[cM_index] += cM_len;

        count++;
        if (count % 100000 == 0)
            fprintf(stderr, "ibdsum done summmarizing %lu\n", count);
    }
    fprintf(stderr, "ibdsum done summmarizing total ibd\n");
}

void
ibdsum_write_binary(ibdsum_t *self, char *outfilename)
{
    FILE *fp = fopen(outfilename, "wb");
    assert(fp != NULL && "can't open file for writing");
    size_t ntotal = self->nsam * (self->nsam - 1) / 2;
    size_t nwrite = 0;

    nwrite = fwrite(self->cM, sizeof(*self->cM), ntotal, fp);
    assert(nwrite == ntotal && "nwrite != ntotal");

    fclose(fp);
    fp = NULL;
    fprintf(stderr, "ibdsum done writing binary\n");
}

void
ibdsum_read_and_add_binary(ibdsum_t *self, char *infilename)
{
    FILE *fp = fopen(infilename, "rb");
    assert(fp != NULL && "can't open file for writing");
    size_t ntotal = self->nsam * (self->nsam - 1) / 2;
    size_t nread = 0;

    nread = fread(self->cM_temp, sizeof(*self->cM_temp), ntotal, fp);
    assert(nread ==  ntotal && "nread != ntotal");

    fclose(fp);
    fp = NULL;

    for (size_t i = 0; i < ntotal; i++) {
        self->cM[i] += self->cM_temp[i];
    }
    fprintf(stderr, "ibdsum done read and add binary\n");
}

void
ibdsum_write_txt_stdout(ibdsum_t *self)
{
    size_t ntotal = self->nsam * (self->nsam - 1) / 2;
    size_t nwrite = 0;
    char *s1 = NULL, *s2 = NULL;
    size_t nsam = self->nsam;
    size_t size = self->sample.size;
    char *name = self->sample.name;
    int ret = 0;

    for (size_t row = 0; row < nsam; row++) {
        s1 = name + row * size;
        for (size_t col = 0; col < nsam; col++) {
            s2 = name + col * size;
            ret = fprintf(stdout, "%s\t%s\t\%.4f\n", s1, s2,
                self->cM[row * (row - 1) / 2 + col] / genomeCM);
            assert(ret >= 0 && "fprintf error");
        }
    }
    fprintf(stderr, "ibdsum done writing stdout\n");
}

int
main(int argc, char *argv[])
{

    ibdsum_t ibdsum;
    char *usage
        = "Usage: zcat sorted_hapibd.gz | "
          "./ibdsum  uniq_samples.txt  out_binary_file  [in_binary_cM_file.dat ...] "
          "> out.txt";
    if (argc < 3) {
        fprintf(stderr, "%s\n", usage);
        exit(1);
    }

    /* read sample file and alloc memory */
    ibdsum_init(&ibdsum, argv[1]);

    /* read input binary file(s) if any */
    if (argc >= 4) {
        for (int i = 3; i < argc; i++) {
            ibdsum_read_and_add_binary(&ibdsum, argv[i]);
        }
    }

    /* summary ibd from stdin */
    ibdsum_summarize_sorted_ibd(&ibdsum);

    /* write summarized ibd to file */
    ibdsum_write_binary(&ibdsum, argv[2]);

    /* write summarized ibd to file */
    ibdsum_write_txt_stdout(&ibdsum);

    ibdsum_free(&ibdsum);
}
