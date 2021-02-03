/* ibdmap:
 * 	calculate genetic positions from physical positions based
 * 	on plink format map file.
 *
 * Input:
 * 	1. a unique list of physicial positions from stdin
 * 		- One non-negative integer number per line
 * 	2. the first command arguments indicating the plink format map file.
 * 		- Assume 4 non-empty columns (tab delimited) per row.
 * 		- 3rd column is the position in cM.
 * 		- 4th column is the position in bp.
 *
 * Output:
 * 	print two-column tab-delimited table to stdout:
 * 	1. 1st column contains positions in bp (from stdin input).
 * 	2. 2nd column contains the corresponding calculated position
 * 	in cM.
 * */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    size_t *bp;
    double *cM;
    size_t nmemb;
    size_t ncap;
    size_t inc;
} maps_t;

int
size_cmp(const void *size1, const void *size2)
{
    if (*((size_t *) size1) > *((size_t *) size2))
        return 1;
    else if (*((size_t *) size1) < *((size_t *) size2))
        return -1;
    else
        return 0;
}

void
maps_free(maps_t *self)
{
    if (self->bp != NULL) {
        free(self->bp);
        self->bp = NULL;
    }
    if (self->cM != NULL) {
        free(self->cM);
        self->cM = NULL;
    }
}

void
maps_read_uniq_positions(maps_t *self, const char *filename)
{
    size_t *bp = NULL;
    double *cM = NULL;
    size_t ncap = 0;
    size_t nmemb = 0;
    size_t inc = 10000;

    char buff[1000], *p = buff;
    size_t buff_size = 1000;
    FILE *fp = stdin;
    int use_stdin = 1;
    if (strcmp(filename, "stdin") != 0) {
        use_stdin = 0;
        fp = fopen(filename, "r");
    }

    while (getline(&p, &buff_size, fp) > 0 && *p != '\n') {
        if (nmemb + 1 >= ncap) {
            // alloc memeory
            ncap += inc;
            bp = realloc(bp, ncap * sizeof(*bp));
            cM = realloc(cM, ncap * sizeof(*cM));
            assert(bp != NULL && "cannot alloc memory");
            assert(cM != NULL && "cannot alloc memory");
        }
        char *tok = strtok(p, "\n");
        bp[nmemb] = strtol(tok, NULL, 10);
        nmemb++;
    }
    ncap = nmemb;
    bp = realloc(bp, ncap * sizeof(*bp));
    cM = realloc(cM, ncap * sizeof(*cM));

    if (use_stdin == 0) {
        fclose(fp);
    }
    fp = NULL;

    // qsort
    qsort(bp, ncap, sizeof(*bp), size_cmp);

    // update self
    self->bp = bp;
    bp = NULL;
    self->cM = cM;
    cM = NULL;
    self->ncap = ncap;
    self->nmemb = nmemb;
    self->inc = 0;
}

int
test_maps_read_uniq_positions()
{
    maps_t maps;
    maps_read_uniq_positions(&maps, "./positions_sorted.txt");
    for (size_t i = 0; i < maps.nmemb; i++)
        printf("%ld\n", maps.bp[i]);
    return 0;
}

void
maps_read_map_file(maps_t *self, const char *filename)
{
    size_t *bp = NULL;
    double *cM = NULL;
    size_t ncap = 0;
    size_t nmemb = 0;
    char *tok = NULL;

    char buff[1000], *p = buff;
    size_t buff_size = 1000;
    FILE *fp = fopen(filename, "r");

    // get size info
    while (getline(&p, &buff_size, fp) > 0 && *p != '\n')
        ncap++;

    // alloc memeory
    bp = calloc(ncap, sizeof(*bp));
    cM = calloc(ncap, sizeof(*cM));
    assert(bp != NULL && "cannot alloc memory");
    assert(cM != NULL && "cannot alloc memory");

    // get actual positions
    fseek(fp, 0, SEEK_SET);
    while (getline(&p, &buff_size, fp) > 0 && *p != '\n') {
        tok = strtok(p, " \t\n"); // chr
        assert(tok != NULL && "map file format error, should use plink format!");
        tok = strtok(NULL, " \t\n"); // snpid
        assert(tok != NULL && "map file format error, should use plink format!");
        tok = strtok(NULL, " \t\n"); // cM, double
        assert(tok != NULL && "map file format error, should use plink format!");
        cM[nmemb] = strtod(tok, NULL);
        tok = strtok(NULL, " \t\n"); // bp, long
        assert(tok != NULL && "map file format error, should use plink format!");
        bp[nmemb] = strtol(tok, NULL, 10);
        nmemb++;
    }

    fclose(fp);
    fp = NULL;

    // update self
    self->bp = bp;
    bp = NULL;
    self->cM = cM;
    cM = NULL;
    self->ncap = ncap;
    self->nmemb = nmemb;
    self->inc = 0;
}

void
maps_update_cM_from_map(maps_t *self, maps_t *ref_map)
{
    size_t ind = 0;
    size_t nmemb = self->nmemb;
    double *cM = self->cM;
    size_t *bp = self->bp;

    size_t ind_ref = 0;
    size_t nmemb_ref = ref_map->nmemb;
    double *cM_ref = ref_map->cM;
    size_t *bp_ref = ref_map->bp;

    double left_ext_slope;
    if (bp[0] == 0) {
        left_ext_slope = (cM_ref[1] - cM_ref[0]) / (bp_ref[1] - bp_ref[0]);
    } else {
        left_ext_slope = cM_ref[0] / bp_ref[0];
    }
    double right_ext_slope = (cM_ref[nmemb_ref - 1] - cM_ref[nmemb_ref - 2])
                             / (bp_ref[nmemb_ref - 1] - bp_ref[nmemb_ref - 2]);

    for (ind = 0; ind < nmemb; ind++) {
        if (bp[ind] <= bp_ref[0]) {
            /* smaller than first position in ref map */
            double y0 = cM_ref[0];
            size_t x0 = bp_ref[0];
            size_t x = bp[ind];
            double y = y0 - left_ext_slope * (x0 - x);
            cM[ind] = y;
        } else if (bp[ind] >= bp_ref[nmemb_ref - 1]) {
            /* greater than last position in ref map */
            double y_max = cM_ref[nmemb_ref - 1]; // fixed a bug here
            size_t x_max = bp_ref[nmemb_ref - 1];
            size_t x = bp[ind];
            double y = y_max + right_ext_slope * (x - x_max);
            cM[ind] = y;

        } else {
            /* between two position in a map */
            // find the range
            while (bp_ref[ind_ref] <= bp[ind])
                ind_ref++;
            double y2 = cM_ref[ind_ref];
            double y1 = cM_ref[ind_ref - 1];
            size_t x2 = bp_ref[ind_ref];
            size_t x1 = bp_ref[ind_ref - 1];
            size_t x = bp[ind];
            assert(
                x >= x1 && x <= x2 && "lhs position should smaller than rhs position");
            double slope = (y2 - y1) / (x2 - x1);
            double y = y1 + (x - x1) * slope;
            cM[ind] = y;
        }
    }
}

int
test_maps_read_map_file()
{
    maps_t maps;
    maps_read_map_file(&maps, "./plink.chr18.GRCh38.map");
    for (size_t i = 0; i < maps.nmemb; i++)
        printf("%ld\t%lf\n", maps.bp[i], maps.cM[i]);
    return 0;
}

int
test_maps_update_cM_from_map()
{
    maps_t maps, refmaps;
    maps_read_uniq_positions(&maps, "./positions_sorted.txt");
    maps_read_map_file(&refmaps, "./plink.chr18.GRCh38.map");
    maps_update_cM_from_map(&maps, &refmaps);
    for (size_t i = 0; i < maps.nmemb; i++)
        printf("%ld\t%lf\n", maps.bp[i], maps.cM[i]);
    return 0;
}

int
main(int argc, char *argv[])
{
    char *usage = "Usage: cat uniq_position_list | ./ibdmap plink_map.map\n";
    if (argc < 2) {
        fprintf(stderr, "%s", usage);
        exit(1);
    }
    maps_t maps, refmaps;
    maps_read_uniq_positions(&maps, "stdin");
    maps_read_map_file(&refmaps, argv[1]);
    maps_update_cM_from_map(&maps, &refmaps);
    for (size_t i = 0; i < maps.nmemb; i++)
        printf("%ld\t%.4f\n", maps.bp[i], maps.cM[i]);
    return 0;
}
