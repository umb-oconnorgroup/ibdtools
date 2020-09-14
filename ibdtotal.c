#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct IDPAIR {
    size_t key;
    size_t id;
};

struct METAINFO {
    char *uniq_labels;
    size_t lsize; // num of chars per label names
    int n_uniq;

    struct IDPAIR *pairs; // key = sample id, id = label id;
    int npairs;
};

// those segments are non-overlapping
struct IBDTOTAL {
    double *cM;    // array of size equal: num_pairs = n*(n-1)/2
    char *samples; // size of n
    int num_sample;
    size_t sample_name_len;
};

struct RGB {
    unsigned char b;
    unsigned char g;
    unsigned char r;
};

struct RGB
RGB_from_hsv(double h, double s, double v)
{
    assert(h >= 0 && h < 360);
    assert(s >= 0 && s <= 1);
    assert(v >= 0 && v <= 1);
    fprintf(stderr, "%g, %g, %g\n", h, s, v);

    double C = v * s;
    double X = C * (1 - abs(((int) h / 60) % 2 - 1));
    double m = v - C;

    double R, G, B; // bug fixed, need to be double not int
    if (h < 60) {
        R = C;
        G = X;
        B = 0;
    } else if (h < 120) {
        R = X;
        G = C;
        B = 0;
    } else if (h < 180) {
        R = 0;
        G = C;
        B = X;
    } else if (h < 240) {
        R = 0;
        G = X;
        B = C;
    } else if (h < 300) {
        R = X;
        G = 0;
        B = C;
    } else if (h < 360) {
        R = C;
        G = 0;
        B = X;
    } else {
    };

    struct RGB color;

    color.r = (unsigned char) ((R + m) * 255);
    color.g = (unsigned char) ((G + m) * 255);
    color.b = (unsigned char) ((B + m) * 255);

    return color;
}

void
RGB_print_RGB_from_hsv(double h, double s, double v)
{
    struct RGB color = RGB_from_hsv(h, s, v);
    printf("hsv: %g\t%g\t%g\tRGB:#%02x%02x%02x\n", h, s, v, color.r, color.g, color.b);
}

void
test_RGB_from_hsv()
{
    RGB_print_RGB_from_hsv(0, 0, 0);
    RGB_print_RGB_from_hsv(0, 0, 1);
    RGB_print_RGB_from_hsv(0, 1, 1);
    RGB_print_RGB_from_hsv(120, 1, 1);
    RGB_print_RGB_from_hsv(240, 1, 1);
    RGB_print_RGB_from_hsv(60, 1, 1);
    RGB_print_RGB_from_hsv(180, 1, 1);
    RGB_print_RGB_from_hsv(300, 1, 1);
    RGB_print_RGB_from_hsv(0, 0, 0.75);
    RGB_print_RGB_from_hsv(0, 0, 0.5);
    RGB_print_RGB_from_hsv(0, 1, 0.5);
    RGB_print_RGB_from_hsv(60, 1, 0.5);
    RGB_print_RGB_from_hsv(120, 1, 0.5);
    RGB_print_RGB_from_hsv(300, 1, 0.5);
    RGB_print_RGB_from_hsv(180, 1, 0.5);
    RGB_print_RGB_from_hsv(240, 1, 0.5);
}

int
main1()
{
    // test_RGB_from_hsv();
    return 0;
}
void
RGB_generate_colors(size_t howmany, struct RGB *colors)
{
}

int
cmp_string(const void *str1, const void *str2)
{
    return strcmp((char *) str1, (char *) str2);
}

int
cmp_double(const void *d1, const void *d2)
{
    return *((double *) d1) > *((double *) d2) ? 1 : -1;
}

int
cmp_idpair(const void *id1, const void *id2)
{
    size_t key1 = ((struct IDPAIR *) id1)->key;
    size_t key2 = ((struct IDPAIR *) id2)->key;
    if (key1 == key2)
        return 0;
    else
        return key1 > key2 ? 1 : -1;
}

/*
 * This function parses the meta table
 * 1) read in a meta table,
 * 2) parse the meta talbe into integer ids and
 * 3) an array of unique labels
 * 4) then sort the sample_id - lable-id by sample_id.
 *
 * scol, lcol are column numbers of sample names and lables, from 1 to ...
 * */
void
IBDTOTAL_meta_read(struct IBDTOTAL *ibdtot, const char *fn_sample_map, int scol,
    int lcol, struct METAINFO *meta)
{
    size_t nsam = ibdtot->num_sample;
    meta->npairs = nsam;

    FILE *fp = fopen(fn_sample_map, "r");
    assert(fp != NULL);

    /* first round read the table to get the size of max len of label names */
    char *linePtr = NULL, *token = NULL;
    size_t lineSize = 0;
    meta->lsize = 0;
    while (getline(&linePtr, &lineSize, fp) > 0 && *linePtr != '\n') {
        token = strtok(linePtr, "\t \n");
        for (int i = 2; i <= lcol; i++) // find the correct column
        {
            token = strtok(NULL, "\t \n");
            if (token == NULL)
                fprintf(stderr, "FATAL: you have empty fields in the meta file %s",
                    fn_sample_map);
        }
        int str_len = strlen(token);
        if (meta->lsize < str_len)
            meta->lsize = str_len;
        fprintf(stderr, "run %d, token: %s \n", __LINE__, token);
    }
    meta->lsize += 1; // considering the additional '\0'

    /* alloc memory */
    char *label_names = calloc(nsam, meta->lsize);
    char *label_names2 = calloc(nsam, meta->lsize);
    meta->pairs = calloc(nsam, sizeof(*meta->pairs));

    /* second round read the actual content of the sample name and labels name
     * parse the sample names to ids based on the sample names in ibdtot */
    fseek(fp, 0, SEEK_SET);
    for (size_t i = 0; i < nsam; i++) {
        assert(getline(&linePtr, &lineSize, fp) > 0);
        int icol = 1;
        for (token = strtok(linePtr, "\t \n"); icol <= (lcol > scol ? lcol : scol);
             icol++) { // bug fixed
            if (icol == lcol) {
                strcpy(label_names + i * meta->lsize, token);
            } else if (icol == scol) {
                // samples are already sorted before this.
                char *p = bsearch(
                    token, ibdtot->samples, nsam, ibdtot->sample_name_len, cmp_string);
                assert(p != NULL);
                meta->pairs[i].key = (p - ibdtot->samples) / ibdtot->sample_name_len;
            } else {
            };
            token = strtok(NULL, "\t \n"); // bug fixed
        };
    }
    free(linePtr);
    linePtr = NULL;

    /* make a copy of labels and sort the label, find a unique list of labels */
    memcpy(label_names2, label_names, nsam * meta->lsize);
    qsort(label_names2, nsam, meta->lsize, cmp_string);

    meta->n_uniq = 1;
    size_t last_uniq_id = 0;

    for (size_t i = 1; i < nsam; i++) {
        if (strcmp(label_names2 + meta->lsize * i,
                label_names2 + meta->lsize * last_uniq_id)
            == 0)
            *(label_names2 + meta->lsize * i) = '\0';
        else {
            last_uniq_id++;
            meta->n_uniq++;
            strcpy(label_names2 + meta->lsize * last_uniq_id,
                label_names2 + meta->lsize * i);
        }
    }

    label_names2 = realloc(label_names2, meta->n_uniq * meta->lsize);
    meta->uniq_labels = label_names2;
    label_names2 = NULL;

    /* based on the uniq l abels array, assign label id to each sample */
    fprintf(stderr, "%ld\n", nsam);
    for (size_t i = 0; i < nsam; i++) {
        // printf("%ld\t%s\n", i, label_names + i * meta->lsize);
        char *p = bsearch(label_names + i * meta->lsize, meta->uniq_labels, meta->n_uniq,
            meta->lsize, cmp_string);
        meta->pairs[i].id = (p - meta->uniq_labels) / meta->lsize;
    }

    free(label_names);
    fclose(fp);
    fp = NULL;

    /* sort id pair by sample id */
    qsort(meta->pairs, nsam, sizeof(*meta->pairs), cmp_idpair);
}

void
IBDTOTAL_make_toy_data(struct IBDTOTAL *ibdtot, const char *fn_toy_meta_to_write)
{
    ibdtot->num_sample = 100;
    ibdtot->sample_name_len = 20;
    size_t num_pairs = ibdtot->num_sample * (ibdtot->num_sample - 1) / 2;
    ibdtot->samples = calloc(ibdtot->num_sample, ibdtot->sample_name_len);
    ibdtot->cM = calloc(num_pairs, sizeof(double));

    for (size_t i = 0; i < ibdtot->num_sample; i++)
        sprintf(ibdtot->samples + i * ibdtot->sample_name_len, "s%ld", i);
    for (size_t i = 0; i < num_pairs; i++)
        ibdtot->cM[i] = rand() % 1000 * 0.1;

    qsort(ibdtot->samples, ibdtot->num_sample, ibdtot->sample_name_len, cmp_string);

    FILE *fp = fopen(fn_toy_meta_to_write, "w");
    assert(fp != NULL);
    for (size_t i = 0; i < ibdtot->num_sample; i++)
        fprintf(fp, "%s\tlabel_%d\n", ibdtot->samples + i * ibdtot->sample_name_len,
            rand() % 8);

    fclose(fp);

    fp = NULL;
}

int
test_IBDTOTAL_meta_read()
{
    struct IBDTOTAL ibdtot;
    struct METAINFO meta;

    IBDTOTAL_make_toy_data(&ibdtot, "toy_meta.txt");

    IBDTOTAL_meta_read(&ibdtot, "toy_meta.txt", 1, 2, &meta);
    printf("stderr: %d\n", ibdtot.num_sample);

    for (size_t i = 0; i < ibdtot.num_sample; i++) {
        fprintf(stderr, "%s\t%ld\t%s\t%ld\n",
            ibdtot.samples + ibdtot.sample_name_len * meta.pairs[i].key,
            meta.pairs[i].key, meta.uniq_labels + meta.lsize * meta.pairs[i].id,
            meta.pairs[i].id);
        fflush(stderr);
    }

    return 0;
}

void
IBDTOTAL_writefile(struct IBDTOTAL *ibdtot, const char *filename)
{
    FILE *fp = fopen(filename, "wb");
    assert(fp != NULL);
    int num_pairs = ibdtot->num_sample * (ibdtot->num_sample - 1) / 2;
    int num_chars = ibdtot->num_sample * ibdtot->sample_name_len;
    fwrite(ibdtot, sizeof(struct IBDTOTAL), 1, fp);
    fwrite(ibdtot->cM, sizeof(double), num_pairs, fp);
    fwrite(ibdtot->samples, sizeof(char), num_chars, fp);
    fclose(fp);
    fp = NULL;
}

// out: ibdtot , in: filename
void
IBDTOTAL_readfile(struct IBDTOTAL *ibdtot, const char *filename)
{
    // read ibdtotal structure for size information
    FILE *fp = fopen(filename, "rb");
    assert(fp != NULL);
    fread(ibdtot, sizeof(struct IBDTOTAL), 1, fp);
    ibdtot->cM = NULL;
    ibdtot->samples = NULL;

    // alloc according to size information
    int num_pairs = ibdtot->num_sample * (ibdtot->num_sample - 1) / 2;
    int num_chars = ibdtot->num_sample * ibdtot->sample_name_len;
    ibdtot->cM = (double *) malloc(sizeof(double) * num_pairs);
    assert(ibdtot->cM != NULL);
    ibdtot->samples = (char *) malloc(sizeof(char) * num_chars);
    assert(ibdtot->samples != NULL);

    // read in arrays according to size info
    fread(ibdtot->cM, sizeof(double), num_pairs, fp);
    fread(ibdtot->samples, sizeof(char), num_chars, fp);

    fclose(fp);
    fp = NULL;
}

void
IBDTOTAL_init(struct IBDTOTAL *ibdtot, const char *samples_file_name)
{
    /* read in sample names */
    FILE *fp = fopen(samples_file_name, "r");
    assert(fp != NULL);
    char buff[100];
    char *p = buff;
    size_t buff_size = 99;

    // determine max sample name len and how many samples;
    ibdtot->sample_name_len = 0;
    ibdtot->num_sample = 0;
    ibdtot->sample_name_len = 0;
    while (getline(&p, &buff_size, fp) > 0) {
        ibdtot->num_sample += 1;
        if (strlen(p) + 1 > ibdtot->sample_name_len)
            ibdtot->sample_name_len = strlen(p) + 1;
        buff_size = 99;
    }
    ibdtot->samples
        = (char *) calloc(ibdtot->num_sample, sizeof(char) * ibdtot->sample_name_len);
    // printf("size(p): %ld\n", ibdtot->sample_name_len);
    // printf("size(p): %d\n", ibdtot->num_smaple);

    // read in the sample sames
    fseek(fp, 0, SEEK_SET);
    for (int i = 0; i < ibdtot->num_sample; i++) {
        p = ibdtot->samples + i * ibdtot->sample_name_len;
        buff_size = ibdtot->sample_name_len;
        assert(getline(&p, &buff_size, fp) > 0);
        // remove '\n'
        for (; *p != '\n'; p++) {
        };
        *p = '\0';
    }
    fclose(fp);
    fp = NULL;

    // sorting sample names
    qsort(ibdtot->samples, ibdtot->num_sample, ibdtot->sample_name_len, cmp_string);

    // allocate mem for cM
    size_t num_pairs = ibdtot->num_sample * (ibdtot->num_sample - 1) / 2;
    ibdtot->cM = (double *) calloc(num_pairs, sizeof(double));
};

void
IBDTOTAL_free(struct IBDTOTAL *ibdtot, const char *samples_file_name)
{
    if (ibdtot->cM != NULL) {
        free(ibdtot->cM);
        ibdtot->cM = NULL;
    }
    if (ibdtot->samples != NULL) {
        free(ibdtot->samples);
        ibdtot->samples = NULL;
    }
}
int
IBDTOTAL_find_smaple_intID(struct IBDTOTAL *ibdtot, const char *str)
{
    int i;
    char *p = (char *) bsearch(
        str, ibdtot->samples, ibdtot->num_sample, ibdtot->sample_name_len, cmp_string);
    i = (p - ibdtot->samples) / ibdtot->sample_name_len;
    return i;
}

int
IBDTOTAL_sampleIds_to_index(struct IBDTOTAL *ibdtot, int id1, int id2)
{

    /* Map two sample ids (in the imaginary upper matrix) to the index of an
     * element for this pair in an array This can save much RAM. _ _ _ _ _ _ | 0 1
     * 1 1 1 1 | | 0 0 1 1 1 1 |  i th row, j th col:  count all 1 s from 0th to
     * i-1 th row | 0 0 0 1 1 1 |                       and then add j minus # of
     * 0 s in the j th current line | 0 0 0 0 1 1 | | 0 0 0 0 0 1 | | 0 0 0 0 0 0
     * |
     *             - - - - - -
     */
    assert(id1 != id2);
    int temp;
    if (id1 > id2) {
        temp = id1;
        id1 = id2;
        id2 = temp;
    }

    return (2 * ibdtot->num_sample - id1 - 1) * id1 / 2 + id2 - id1 - 1;
}

void
IBDTOTAL_add_cM(struct IBDTOTAL *ibdtot, const char *str1, const char *str2, double cM)
{
    int id1 = IBDTOTAL_find_smaple_intID(ibdtot, str1);
    int id2 = IBDTOTAL_find_smaple_intID(ibdtot, str2);
    int idtemp;
    if (id1 > id2) {
        idtemp = id2;
        id2 = id1;
        id1 = idtemp;
    }

    // printf("%d\t%s\t%d\t%s\t%lf\n", id1, str1, id2, str2, cM);
    int index = IBDTOTAL_sampleIds_to_index(ibdtot, id1, id2);
    ibdtot->cM[index] += cM;
}

double
IBDTOTAL_get_cM(struct IBDTOTAL *ibdtot, int id1, int id2)
{
    int idtemp;
    if (id1 > id2) {
        idtemp = id2;
        id2 = id1;
        id1 = idtemp;
    }

    int index = IBDTOTAL_sampleIds_to_index(ibdtot, id1, id2);
    return ibdtot->cM[index];
}

void
IBDTOTAL_print(struct IBDTOTAL *ibdtot)
{
    /*
    for(int i =0; i <ibdtot->num_smaple; i++)
    {
            printf("%s\n", ibdtot->samples + i * ibdtot->sample_name_len);
    }
    */

    // printf("-----------\n");

    double cM;
    double genome_size = 3545.8;
    for (int i = 0; i < ibdtot->num_sample - 1; i++)
        for (int j = i + 1; j < ibdtot->num_sample; j++) {
            cM = IBDTOTAL_get_cM(ibdtot, i, j);
            if (cM < 1.9)
                continue;
            // printf("%d\t%s\t%d\t%s\t%lf\n",
            printf("%s\t%s\t%g\t%g\n",
                // i,
                ibdtot->samples + i * ibdtot->sample_name_len,
                // j,
                ibdtot->samples + j * ibdtot->sample_name_len, cM, cM / genome_size);
        }
}

/*
 *   https://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries
 */
int
WriteIMAGE(const char *filename, int w, int h, unsigned char *img)
{
    FILE *f;
    int filesize = 54 + 3 * w * h; // w is your image width, h is image height, both int

    unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };
    unsigned char bmpinfoheader[40]
        = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0 };
    unsigned char bmppad[3] = { 0, 0, 0 };

    bmpfileheader[2] = (unsigned char) (filesize);
    bmpfileheader[3] = (unsigned char) (filesize >> 8);
    bmpfileheader[4] = (unsigned char) (filesize >> 16);
    bmpfileheader[5] = (unsigned char) (filesize >> 24);

    bmpinfoheader[4] = (unsigned char) (w);
    bmpinfoheader[5] = (unsigned char) (w >> 8);
    bmpinfoheader[6] = (unsigned char) (w >> 16);
    bmpinfoheader[7] = (unsigned char) (w >> 24);
    bmpinfoheader[8] = (unsigned char) (h);
    bmpinfoheader[9] = (unsigned char) (h >> 8);
    bmpinfoheader[10] = (unsigned char) (h >> 16);
    bmpinfoheader[11] = (unsigned char) (h >> 24);

    f = fopen(filename, "wb");
    fwrite(bmpfileheader, 1, 14, f);
    fwrite(bmpinfoheader, 1, 40, f);
    for (int i = 0; i < h; i++) {
        fwrite(img + (w * (h - i - 1) * 3), 3, w, f);
        fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
    }
    fclose(f);

    return 0;
}

// find the percentile in ibdtot->cM
// 1. first copy to cM_temp array;
// 2. sort the temp array
// 3. first the first non zero ( >1)
// 3. use percentage to calculate index
// 4. get percentile by index
// Note: para: percentage: 0.0 - 100.0
double
IBDTOTAL_get_cM_percentile(struct IBDTOTAL *ibdtot, double percentage)
{
    double percentile = 0, *cM_temp = NULL;
    size_t index = 0;
    size_t num_zeros = 0, i, num_nonzeros;
    size_t num_pairs = ibdtot->num_sample * (ibdtot->num_sample - 1) / 2;

    assert(percentage >= 0);
    assert(percentage <= 100);

    cM_temp = (double *) malloc(sizeof(double) * num_pairs);
    assert(cM_temp != NULL);
    memcpy(cM_temp, ibdtot->cM, sizeof(double) * num_pairs);

    qsort(cM_temp, num_pairs, sizeof(double), cmp_double);

    assert(cM_temp[0] >= 0);

    // calculate num_zeros
    for (size_t i = 0; i < num_pairs && cM_temp[i] < 0.1; i++, num_zeros++) {
    };

    num_nonzeros = num_pairs - num_zeros;

    index = (int) (percentage * num_nonzeros / 100);
    if (index >= num_nonzeros)
        index = num_nonzeros - 1;

    assert(index < num_nonzeros);

    index += num_zeros;
    percentile = cM_temp[index];

    // output_percentile
    int num_perc = 1000;

    FILE *fp = fopen("res_ibdtotal_quantile.txt", "w");
    assert(fp != NULL);
    for (int i = 0; i <= num_perc; i++) {
        double perc = 1.0 * i / num_perc;
        size_t index = (size_t)(num_pairs * perc);
        if (index >= num_pairs)
            index = num_pairs - 1;
        fprintf(fp, "%ld\t%g\t%g\n", index, perc, cM_temp[index]);
    }
    fclose(fp);
    fp = NULL;

    // output_histogram

    fp = fopen("res_ibdtotal_hist.txt", "w");
    assert(fp != NULL);
    double bin_upper = 1;
    size_t counter = 0;
    for (size_t i = 0; i < num_pairs; i++) {
        if (cM_temp[i] <= bin_upper)
            counter++;
        else {
            fprintf(fp, "%g\t%g\t%ld\n", bin_upper - 1, bin_upper, counter);
            counter = 0;
            bin_upper += 1;
            i--;
        }
    }
    fclose(fp);
    fp = NULL;
    free(cM_temp);
    cM_temp = NULL;

    fprintf(stderr, "%lf\n", percentile);
    return percentile;
}

int
IBDTOTAL_generate_map(
    struct IBDTOTAL *ibdtot, char *filename, double cutoff_cM, int *sample_orders)
{
    int w = ibdtot->num_sample / 3;
    int h = w;
    int index, i1, i2;
    double value, temp;
    double values9[9];
    int is_saturated;
    unsigned char *img = (unsigned char *) calloc(w * h, sizeof(*img) * 3);

    for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++) {
            int touched_diag = 0;
            // 1. use the top left cell
            // index = IBDTOTAL_sampleIds_to_index(ibdtot, 3*i, 3*j);
            // value = ibdtot->cM[index];

            // 2. use the median of the nine cells
            // save values of 9 cells to an array, and then sort
            for (i1 = 0; i1 < 3; i1++)
                for (i2 = 0; i2 < 3; i2++) {
                    int row = 3 * i + i1;
                    int col = 3 * j + i2;
                    if (sample_orders != NULL) {
                        row = sample_orders[row];
                        col = sample_orders[col];
                    }

                    if (row == col) {
                        touched_diag = 1;
                        break;
                    }
                    index = IBDTOTAL_sampleIds_to_index(ibdtot, row, col);
                    values9[3 * i1 + i2] = ibdtot->cM[index];
                }

            if (touched_diag)
                value = 0;
            else
                // qsort(values9,9,sizeof(double),cmp_double);
                value = values9[4]; // median

            // tranform the value ; cutoff value is a percentile number
            value = sqrt(value / cutoff_cM) * 255;
            if (value > 255) {
                value = 255;
                is_saturated = 1;
            } else {
                is_saturated = 0;
            }
            value = 255 - value; // high will be darker

            // make pixel
            // if (i > j) // upper triangle
            if (0) {
                img[3 * (i * w + j) + 2] = (unsigned char) 255; // r
                img[3 * (i * w + j) + 1] = (unsigned char) 255; // g
                img[3 * (i * w + j) + 0] = (unsigned char) 255; // b
            } else if (is_saturated == 0) {
                img[3 * (i * w + j) + 2] = (unsigned char) value; // r
                img[3 * (i * w + j) + 1] = (unsigned char) value; // g
                img[3 * (i * w + j) + 0] = (unsigned char) value; // b
            } else {
                img[3 * (i * w + j) + 2] = (unsigned char) 255; // r
                img[3 * (i * w + j) + 1] = (unsigned char) 0;   // g
                img[3 * (i * w + j) + 0] = (unsigned char) 0;   // b
            }
        }

    WriteIMAGE(filename, w, h, img);

    free(img);
    return 0;
}

// TODO: work on color pallete
void
IBDTOTAL_generate_map_legend(
    struct IBDTOTAL *ibdtot, const char *image_fn, int *sample_orders)
{
    int num_samples = ibdtot->num_sample;
    struct RGB {
        unsigned char b;
        unsigned char g;
        unsigned char r;
    };

    struct RGB color_8v1_A_NAD = { 176, 121, 44 };   // 1148
    struct RGB color_8v1_B_NAD = { 54, 127, 251 };   // 880
    struct RGB color_C1_hisp = { 59, 158, 57 };      // 463
    struct RGB color_Emerge_ids = { 53, 45, 210 };   // 4852
    struct RGB color_LARGEPD = { 185, 109, 147 };    // 1501
    struct RGB color_NWD = { 78, 86, 138 };          // 17560
    struct RGB color_PAGEII = { 192, 125, 223 };     // 5213
    struct RGB color_SIGMA_OMNI = { 127, 127, 127 }; // 1147

    int width = 5000, height = num_samples;
    unsigned char *img = malloc(3 * width * height * sizeof(*img));

    for (int row = 0; row < num_samples; row++) {
        int index;
        if (sample_orders != NULL)
            index = sample_orders[row];
        else
            index = row;
        struct RGB color = { 255, 255, 255 };
        if (index < 1148)
            color = color_8v1_A_NAD;
        else if (index < 1148 + 880)
            color = color_8v1_B_NAD;
        else if (index < 1148 + 880 + 463)
            color = color_C1_hisp;
        else if (index < 1148 + 880 + 463 + 4852)
            color = color_Emerge_ids;
        else if (index < 1148 + 880 + 463 + 4852 + 1501)
            color = color_LARGEPD;
        else if (index < 1148 + 880 + 463 + 4852 + 1501 + 17560)
            color = color_NWD;
        else if (index < 1148 + 880 + 463 + 4852 + 1501 + 17560 + 5213)
            color = color_PAGEII;
        else if (index < 1148 + 880 + 463 + 4852 + 1501 + 17560 + 5213 + 1147)
            color = color_SIGMA_OMNI;
        else {
        }
        for (int col = 0; col < width; col++)
            memcpy(img + 3 * (row * width + col), &color, sizeof(struct RGB));
        //  printf("index=%d, r=%uz, g=%uz, b=%uz\n",index,  color.r, color.g,
        //  color.b);
    }
    WriteIMAGE(image_fn, width, height, img);
}
void
IBDTOTAL_generate_label_matrix(struct IBDTOTAL *ibdtot, const char *image_fn,
    struct METAINFO *meta, int *sample_orders)

{
    int num_samples = ibdtot->num_sample;
    int num_groups = meta->n_uniq;
    size_t height = num_samples;
    size_t gr_width = 200, color_width = 150;
    size_t width = gr_width * num_groups;

    unsigned char *img = malloc(3 * width * height * sizeof(*img));
    unsigned char *line = malloc(3 * color_width * 1 * sizeof(*img));

    memset(
        img, (unsigned char) 255, 3 * width * height * sizeof(*img)); // white background
    memset(line, 0, 3 * color_width * 1 * sizeof(*img));              // black line

    for (int row = 0; row < num_samples; row++) {
        int index;
        if (sample_orders != NULL)
            index = sample_orders[row];
        else
            index = row;
        size_t label_id = meta->pairs[index].id;

        for (int g = 0; g < num_groups; g++)
            if (label_id == g)
                memcpy(img + 3 * (row * width + g * gr_width), line,
                    3 * color_width * 1 * sizeof(*img));
    }
    WriteIMAGE(image_fn, width, height, img);
    free(img);
    free(line);
    img = NULL;
    line = NULL;

    size_t *label_counts = calloc(meta->n_uniq, sizeof(*label_counts));
    for (size_t i = 0; i < meta->npairs; i++)
        label_counts[meta->pairs[i].id] += 1;

    char newname[100];
    memset(newname, 100, sizeof(char));
    strcpy(newname, image_fn);
    strtok(newname, ".");
    strcat(newname, "_block_label.txt");
    FILE *fp = fopen(newname, "w");
    assert(fp != NULL);
    for (size_t block = 0; block < meta->n_uniq; block++)
        fprintf(fp, "%s\t%ld\n", meta->uniq_labels + block * meta->lsize,
            label_counts[block]);
    fclose(fp);
    fp = NULL;

    free(label_counts);
}

void
IBDTOTAL_generate_label_bar(struct IBDTOTAL *ibdtot, const char *image_fn,
    struct METAINFO *meta, int *sample_orders)
{

    int num_samples = ibdtot->num_sample;

    int width = 1000, height = num_samples;
    unsigned char *img = malloc(3 * width * height * sizeof(*img));

    for (int row = 0; row < num_samples; row++) {
        int index;
        if (sample_orders != NULL)
            index = sample_orders[row];
        else
            index = row;
        size_t label_id = meta->pairs[index].id;
        struct RGB color
            = RGB_from_hsv(360.0 * label_id / meta->n_uniq, 1, 0.5 * (label_id % 2));

        for (int col = 0; col < width; col++)
            memcpy(img + 3 * (row * width + col), &color, sizeof(struct RGB));
        //  printf("index=%d, r=%uz, g=%uz, b=%uz\n",index,  color.r, color.g,
        //  color.b);
    }
    WriteIMAGE(image_fn, width, height, img);
    free(img);
    img = NULL;

    // write label block image
    unsigned char *img_sorted_by_label = malloc(3 * width * height * sizeof(*img));
    size_t *label_counts = calloc(meta->n_uniq, sizeof(*label_counts));
    for (size_t i = 0; i < meta->npairs; i++)
        label_counts[meta->pairs[i].id] += 1;

    size_t nrow = 0;
    for (size_t block = 0; block < meta->n_uniq; block++) {
        struct RGB color
            = RGB_from_hsv(360.0 * block / meta->n_uniq, 1, 0.5 * (block % 2));
        for (size_t pixel = 0; pixel < width * label_counts[block]; pixel++)
            memcpy(img_sorted_by_label + 3 * (nrow * width + pixel), &color,
                sizeof(struct RGB));
        nrow += label_counts[block];
    }

    char newname[100];
    strcpy(newname, image_fn);
    strtok(newname, ".");
    strcat(newname, "_block.bmp");
    WriteIMAGE(newname, width, height, img_sorted_by_label);

    // write unique lable list in the same order of the color block
    memset(newname, 100, sizeof(char));
    strcpy(newname, image_fn);
    strtok(newname, ".");
    strcat(newname, "_block_label.txt");
    FILE *fp = fopen(newname, "w");
    assert(fp != NULL);
    for (size_t block = 0; block < meta->n_uniq; block++)
        fprintf(fp, "%s\t%ld\n", meta->uniq_labels + block * meta->lsize,
            label_counts[block]);
    fclose(fp);
    fp = NULL;

    free(img_sorted_by_label);
    img_sorted_by_label = NULL;
    free(label_counts);
    label_counts = NULL;
}

void
test_IBDTOTAL_generate_label_bar()
{
    struct IBDTOTAL ibdtot;
    IBDTOTAL_make_toy_data(&ibdtot, "toy_meta.txt");
    struct METAINFO meta;
    IBDTOTAL_meta_read(&ibdtot, "toy_meta.txt", 1, 2, &meta);

    FILE *fp_order = fopen("./sample_orders.txt", "r"); // argc[3] is sample order files
    assert(fp_order != NULL);
    int *sample_orders = calloc(ibdtot.num_sample, sizeof(*sample_orders));
    assert(sample_orders != NULL);
    for (int i = 0; i < ibdtot.num_sample; i++) {
        fscanf(fp_order, "%d", sample_orders + i);
        // fprintf(stderr, "%d " , sample_orders[i]);
    }
    fclose(fp_order);
    fp_order = NULL;
    IBDTOTAL_generate_label_bar(&ibdtot, "1.bmp", &meta, sample_orders);
}

int
main2()
{
    test_IBDTOTAL_generate_label_bar();
    return 0;
}

int
main(int argc, char *argv[])
{
    struct IBDTOTAL ibdtotal;

    char buff[1000], *p = buff, *token = NULL;
    size_t size = 1000;
    size_t lineCount = 0;
    int colCount = 0;
    char *s1, *s2;
    double cM, cutoff_cM;
    char *beg = NULL, *sep = NULL;
    int chr;

    if (argc != 2 && argc != 7) {
        fprintf(stderr, 
			"\nUsage 1: read ibd from stdin, stdout = ibdtotal," 
			" also save binary matrix to tmp_ibd.dat\n\n"
			" 	cat xxx.ibd | ./ibdtotal <samples_list file>  >ibd_total.txt\n\n"
			"\nUsage 2: read from binary matrix file\n\n"
			"	./ibdtotal <samples_list file> [<matrix_file> "
                        "<samples_orders_file> <meta_table> <scol> <lcol>]\n\n");
        exit(1);
    }
    /* read merged ibd from stdin*/
    else if (argc == 2) {
        IBDTOTAL_init(&ibdtotal, argv[1]); // argv[1] = samples_list file

        // parse lines and save it to cM, sample pairs into arrays
        while (getline(&p, &size, stdin) > 0) {
            lineCount++;
            // if(lineCount % 10000 == 0) fprintf(stderr, "line: %ld\n", lineCount);
            token = strtok(p, "\t");
            colCount = 0;
            while (colCount <= 9) {
                colCount++;
                if (colCount == 1)
                    s1 = token;
                else if (colCount == 3)
                    s2 = token;
                else if (colCount == 5)
                    chr = atoi(token);
                else if (colCount == 9)  // update cM to 9th col, score to be 8th col
                    cM = strtod(token, NULL);
                else {
                };
                token = strtok(NULL, "\t");
            }
            IBDTOTAL_add_cM(&ibdtotal, s1, s2, cM);
        }
        // write to file
        IBDTOTAL_writefile(&ibdtotal, "tmp_ibd.dat");
        // call print functions
        IBDTOTAL_print(&ibdtotal);
    }
    /* read ibdtotal from matrix*/
    else if (argc == 7) {
        IBDTOTAL_init(&ibdtotal, argv[1]); // argv[1] = samples_list file
        fprintf(stderr, "reading file\n");
        IBDTOTAL_readfile(&ibdtotal, argv[2]);
        fprintf(stderr, "done reading file\n");

        int *sample_orders = NULL;
        // reading order files
        FILE *fp_order = fopen(argv[3], "r"); // argc[3] is sample order files
        assert(fp_order != NULL);
        sample_orders = calloc(ibdtotal.num_sample, sizeof(*sample_orders));
        assert(sample_orders != NULL);
        for (int i = 0; i < ibdtotal.num_sample; i++) {
            fscanf(fp_order, "%d", sample_orders + i);
            // fprintf(stderr, "%d " , sample_orders[i]);
        }
        fclose(fp_order);
        fp_order = NULL;

        fprintf(stderr, "start\n");

        // calculate percentile
        // --- cutoff_cM = IBDTOTAL_get_cM_percentile(&ibdtotal, 90);
        // --- fprintf(stderr, "percentage= %lf, percentile %lf\n", 90.0, cutoff_cM);
        // --- // cutoff_cM = 50;
        // --- // generate bmp
        // --- IBDTOTAL_generate_map(&ibdtotal, "img_raw.bmp", cutoff_cM, NULL);
        // --- if (sample_orders != NULL)
        // ---   IBDTOTAL_generate_map(&ibdtotal, "img_clust.bmp", cutoff_cM,
        // ---                         sample_orders);

        // generate sample color
        struct METAINFO meta;
        fprintf(stderr, "meta_read_start\n");
        IBDTOTAL_meta_read(&ibdtotal, argv[4], atoi(argv[5]), atoi(argv[6]), &meta);
        fprintf(stderr, "meta_read_end\n");
        fprintf(stderr, "generate_label_bar start\n");
        IBDTOTAL_generate_label_bar(
            &ibdtotal, "img_clust_label.bmp", &meta, sample_orders);
        fprintf(stderr, "generate_label_bar end\n");

        // generate label matrix
        IBDTOTAL_generate_label_matrix(
            &ibdtotal, "img_clust_label_matrix.bmp", &meta, sample_orders);

        free(sample_orders);
        sample_orders = NULL;
    }

    return 0;
}
