#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

typedef struct {
    long bp_per_cM;
    double time_scale;
    long genome_size_cM;
} msmeta_t;

typedef struct {
    int nsam;
    int nrep;
    double rho;
    double theta;
    long nsites;

} mscmd_t;

typedef struct {
    long nsegsites;
    long *positions;
    char *matrix;
    int nsam;
    long chr_length;
} mschr_t;

int
msmeta_set_bp_per_cM(msmeta_t *self, long bp_per_cM)
{
    self->bp_per_cM = bp_per_cM;
    return 0;
}
int
mschr_alloc(mschr_t *self, int nsam, long nseg, long nsites)
{
    self->chr_length = nsites;
    self->nsam = nsam;
    self->nsegsites = nseg;
    self->matrix = calloc(self->nsam * self->nsegsites, sizeof(char));
    self->positions = calloc(self->nsegsites, sizeof(long));
    return 0;
}
int
mschr_free(mschr_t *self)
{
    if (self->matrix)
        free(self->matrix);
    if (self->positions)
        free(self->positions);
    self->matrix = NULL;
    self->positions = NULL;
    return 0;
}
typedef struct {
    msmeta_t *meta;
    mscmd_t *cmd;
    mschr_t *chrs;
    int num_chr;
} msfile_t;

int
mscmd_parse(mscmd_t *self, char *mscmd)
{
    int ret = 0;
    char *tok = NULL;
    tok = strtok(mscmd, " \n");
    while (tok != NULL) {
        if (strcmp(tok, "ms") == 0) {
            tok = strtok(NULL, " \n");
            assert(tok != NULL);
            if (strcmp(tok, "-threads") == 0) {
                /* do this to skip the threads */
                tok = strtok(NULL, " \n");
                assert(tok != NULL);
            	tok = strtok(NULL, " \n");
            	assert(tok != NULL);
            }
            self->nsam = atoi(tok);

            tok = strtok(NULL, " \n");
            assert(tok != NULL);
            self->nrep = atoi(tok);
        } else if (strcmp(tok, "-threads") == 0) {
	    /* do this to skip the threads */
            tok = strtok(NULL, " \n");
            assert(tok != NULL);
        } else if(strcmp(tok, "-t") ==0 ) {
            tok = strtok(NULL, " \n");
            assert(tok != NULL);
            self->theta = strtod(tok, NULL);
            assert(self->theta > 0);

	} else if (strcmp(tok, "-r") == 0) {
            tok = strtok(NULL, " \n");
            assert(tok != NULL);
            self->rho = strtod(tok, NULL);
            assert(self->rho > 0);

            tok = strtok(NULL, " \n");
            assert(tok != NULL);
            self->nsites = strtol(tok, NULL, 10);
            assert(self->nsites > 0);
        } else {
        }
        tok = strtok(NULL, " \n");
    }
    return ret;
}

void
test_mscmd_parse()
{
    char mscmd[] = "ms 300 14 -t 300.0 -r 19999.973333333335 750000 -seed 212943  "
                   "[3.2rc Build:162]\n";
    mscmd_t cmd;
    mscmd_parse(&cmd, mscmd);
    printf("nsam = %d, nrep = %d, theta = %g, rho = %g, nsites = %ld \n", cmd.nsam,
        cmd.nrep, cmd.theta, cmd.rho, cmd.nsites);
}

int
msfile_alloc(msfile_t *self)
{
    self->cmd = calloc(1, sizeof(mscmd_t));
    self->meta = calloc(1, sizeof(msmeta_t));
    self->chrs = NULL;
    self->num_chr = 0;
    return 0;
};

int
msfile_free(msfile_t *self)
{
    if (self->cmd)
        free(self->cmd);
    if (self->meta)
        free(self->meta);
    if (self->chrs) {
        for (int i = 0; i < self->num_chr; i++) {
            mschr_free(self->chrs + i);
        }
        free(self->chrs);
    }
    self->cmd = NULL;
    self->meta = NULL;
    self->chrs = NULL;
    return 0;
};

int
msfile_parser(msfile_t *self, FILE *msfile)
{
    char *linePtr = NULL;
    size_t buf_size = 0;
    mschr_t *cur_chr = NULL;
    char *cur_allel = NULL;
    char *tok = NULL;
    long nseg;

    while (getline(&linePtr, &buf_size, msfile) > 0) {
        if (strncmp(linePtr, "ms", 2) == 0) {
            mscmd_parse(self->cmd, linePtr);
            self->num_chr = self->cmd->nrep;
            self->chrs = calloc(self->num_chr, sizeof(mschr_t));
        } else if (strncmp(linePtr, "//", 2) == 0) {
            if (cur_chr == NULL) {
                cur_chr = self->chrs;
            } else {
                cur_chr++;
            }
            assert(cur_chr <= self->chrs + self->num_chr - 1);
        } else if (strncmp(linePtr, "[", 1) == 0) {
            // ignore for now
            ;
        } else if (strncmp(linePtr, "segsites", 8) == 0) {
            tok = strtok(linePtr, " ");
            tok = strtok(NULL, " ");
            nseg = strtol(tok, NULL, 10);
            assert(nseg > 0);
            mschr_alloc(cur_chr, self->cmd->nsam, nseg, self->cmd->nsites);
            cur_allel = cur_chr->matrix;
        } else if (strncmp(linePtr, "positions", 9) == 0) {
            tok = strtok(linePtr, " ");
            tok = strtok(NULL, " ");
            for (long i = 0; i < cur_chr->nsegsites; i++) {
                assert(tok != NULL);
                cur_chr->positions[i] = (long) (strtod(tok, NULL) * cur_chr->chr_length);
                // to avoid collision
                if (i >= 1) {
                    if (cur_chr->positions[i] <= cur_chr->positions[i - 1]) {
                        cur_chr->positions[i] = cur_chr->positions[i - 1] + 1;
                    }
                    assert(cur_chr->positions[i] > 0);
                    tok = strtok(NULL, " ");
                }
	    }
        }
        // need to distinguish allele from the random number line
        else if ((*linePtr == '0' || *linePtr == '1') && cur_chr != NULL) {
            assert(strlen(linePtr) >= cur_chr->nsegsites);
            memcpy(cur_allel, linePtr, cur_chr->nsegsites);

            cur_allel += cur_chr->nsegsites;
        } else {
            continue;
        }
    }
    return 0;
}

int
msfile_print_vcf(msfile_t *self, FILE *f_vcf, int het, FILE *f_map)
{
    mschr_t *cur_chr = NULL;
    int chrno = 0;
    char ATGC[] = "ATGC";
    int ref, alt;
    long position;
    /* header */
    fprintf(f_vcf, "##fileformat=VCFv4.2\n"
                   "##source=macs_tree_parser by B.G.\n");
    for (cur_chr = self->chrs; cur_chr < self->chrs + self->num_chr; cur_chr++) {
        fprintf(f_vcf, "##contig=<ID=%ld,length=%ld>\n", cur_chr - self->chrs + 1,
            cur_chr->chr_length);
    }
    fprintf(f_vcf, "##INFO=<ID=PR,Number=1,Type=Flag,Description=\"coverted from ms "
                   "format data\">\n"
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    if (het != 0) { // het
        assert(self->cmd->nsam % 2 == 0);
        for (long i = 0; i < self->cmd->nsam / 2; i++)
            fprintf(f_vcf, "\t%ld", i);
    } else { // hom
        for (long i = 0; i < self->cmd->nsam; i++)
            fprintf(f_vcf, "\t%ld", i);
    }
    fprintf(f_vcf, "\n");

    /* body */
    for (cur_chr = self->chrs; cur_chr <= self->chrs + self->num_chr - 1; cur_chr++) {
        chrno = cur_chr - self->chrs + 1;
        for (long site = 0; site < cur_chr->nsegsites; site++) {
            position = cur_chr->positions[site];
            ref = rand() % 4;
            alt = ref;
            while (alt == ref)
                alt = rand() % 4;
            //	printf intial columns
            fprintf(f_vcf, "%d\t%ld\tsnp_%d_%ld\t%c\t%c\t.\t.\tPR\tGT", chrno, position,
                chrno, position, ATGC[ref], ATGC[alt]);
            // ** for map file
            fprintf(f_map, "%d\tsnp_%d_%ld\t%g\t%ld\n", chrno, chrno, position,
                1.0 * position / self->meta->bp_per_cM, position);

            //	printf allele columns
            if (het != 0) // het
                for (int hapid = 0; hapid < cur_chr->nsam; hapid += 2)
                    fprintf(f_vcf, "\t%c|%c",
                        cur_chr->matrix[hapid * cur_chr->nsegsites + site],
                        cur_chr->matrix[(hapid + 1) * cur_chr->nsegsites + site]);
            else // homo
                for (int hapid = 0; hapid < cur_chr->nsam; ++hapid)
                    fprintf(f_vcf, "\t%c|%c",
                        cur_chr->matrix[hapid * cur_chr->nsegsites + site],
                        cur_chr->matrix[hapid * cur_chr->nsegsites + site]);
            fprintf(f_vcf, "\n");
        }
    }
    return 0;
}

void
test_msfile_parser()
{
    FILE *fp = fopen("./test.ms", "r");
    msfile_t msfile;
    msfile_alloc(&msfile);
    msfile_parser(&msfile, fp);
    msfile_print_vcf(&msfile, stdout, 0, stderr);
    msfile_free(&msfile);
    fclose(fp);
    return;
}
int
main(int argc, char *argv[])
{
    msfile_t msfile;
    long bp_per_cM;
    int is_het;
    if (argc != 3) {
        fprintf(stderr, " Usage:\n\n  msms .... | ms2vcf <bp_per_cM> <het=1, hom=0>");
        exit(1);
    }
    bp_per_cM = strtol(argv[1], NULL, 10);
    assert(bp_per_cM > 0);
    is_het = atoi(argv[2]);
    msfile_alloc(&msfile);
    msmeta_set_bp_per_cM(msfile.meta, bp_per_cM);
    //msfile_parser(&msfile, fp);
    msfile_parser(&msfile, stdin);
    msfile_print_vcf(&msfile, stdout, is_het, stderr);
    msfile_free(&msfile);
    return 0;
}
