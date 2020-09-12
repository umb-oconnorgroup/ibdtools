#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int str_cmp(const void *str1, const void *str2) {
  return strcmp((char *)str1, (char *)str2);
}

int size_cmp(const void *size1, const void *size2) {
  if (*((size_t *)size1) > *((size_t *)size2))
    return 1;
  else if (*((size_t *)size1) < *((size_t *)size2))
    return -1;
  else
    return 0;
}

typedef struct {
  char *names;
  int name_size;
  size_t nmemb;
  size_t ncap;
  size_t inc;
} samples_t;

void samples_read(samples_t *self, const char *filename) {
  size_t name_size = 0;
  size_t nmemb;
  size_t ncap = 0;
  char *names = NULL;

  char buff[1000], *p = buff;
  size_t buff_size = 1000;
  FILE *fp = fopen(filename, "r");
  assert(fp != NULL);

  // get size info
  while (getline(&p, &buff_size, fp) > 0 && *p != '\n') {
    int chars_used = strlen(p) + 1;
    if (chars_used > name_size)
      name_size = chars_used;
    ncap++;
  }

  // alloc memory
  names = calloc(name_size * ncap, sizeof(*names));

  // read in samples_names
  fseek(fp, 0, SEEK_SET);
  for (size_t i = 0; i < ncap; i++) {
    int ret;
    p = names + i * name_size;
    ret = getline(&p, &name_size, fp);
    assert(ret > 0);

    // remove '\n'
    for (char *p2 = p; p2 < p + name_size; p2++)
      if (*p2 == '\n') {
        *p2 = '\0';
        break;
      }
  }

  fclose(fp);
  fp = NULL;

  // sort: linux sort until genrate different order from c strcmp
  qsort(names, ncap, name_size, str_cmp);

  // update self
  self->names = names;
  names = NULL;
  self->name_size = name_size;
  self->ncap = ncap;
  self->nmemb = ncap;
  self->inc = 1;
}

int test_sample_reads() {
  char *filename = "./samples_sorted.txt";
  samples_t samples;
  samples_read(&samples, filename);
  for (size_t i = 0; i < samples.nmemb; i++)
    printf("%s\n", samples.names + i * samples.name_size);
  return 0;
}

typedef struct {
  size_t *bp;
  double *cM;
  size_t nmemb;
  size_t ncap;
  size_t inc;
} maps_t;

void maps_free(maps_t *self) {
  if (self->bp != NULL) {
    free(self->bp);
    self->bp = NULL;
  }
  if (self->cM != NULL) {
    free(self->cM);
    self->cM = NULL;
  }
}

void maps_read_uniq_positions(maps_t *self, const char *filename) {
  size_t *bp = NULL;
  double *cM = NULL;
  size_t ncap = 0;

  char buff[1000], *p = buff;
  size_t buff_size = 1000;
  FILE *fp = fopen(filename, "r");

  // get size info
  while (getline(&p, &buff_size, fp) > 0 && *p != '\n')
    ncap++;

  // alloc memeory
  bp = calloc(ncap, sizeof(*bp));
  cM = calloc(ncap, sizeof(*cM));
  assert(bp != NULL);
  assert(cM != NULL);

  // get actual positions
  fseek(fp, 0, SEEK_SET);
  for (size_t i = 0; i < ncap; i++)
    fscanf(fp, "%ld", bp + i);

  fclose(fp);
  fp = NULL;

  // qsort
  qsort(bp, ncap, sizeof(*bp), size_cmp);

  // some vcf have repeat position, which will fail the assertion. comment it out
  // there could be some duplication but will affection the lookup of poistions
  // for (size_t i = 1; i < ncap; i++)
  // {
  //   fprintf(stderr, "%ld\t%ld\n", bp[i-1], bp[i]);
  //   assert(bp[i - 1] < bp[i]);
  // }

  // update self
  self->bp = bp;
  bp = NULL;
  self->cM = cM;
  cM = NULL;
  self->ncap = ncap;
  self->nmemb = ncap;
  self->inc = 0;
}

int test_maps_read_uniq_positions() {
  maps_t maps;
  maps_read_uniq_positions(&maps, "./positions_sorted.txt");
  for (size_t i = 0; i < maps.nmemb; i++)
    printf("%ld\n", maps.bp[i]);
  return 0;
}

void maps_read_map_file(maps_t *self, const char *filename) {
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
  assert(bp != NULL);
  assert(cM != NULL);

  // get actual positions
  fseek(fp, 0, SEEK_SET);
  while (getline(&p, &buff_size, fp) > 0 && *p != '\n') {
    tok = strtok(p, " \t\n");    // chr
    tok = strtok(NULL, " \t\n"); // snpid
    tok = strtok(NULL, " \t\n"); // cM, double
    assert(tok != NULL);
    cM[nmemb] = strtod(tok, NULL);
    tok = strtok(NULL, " \t\n"); // bp, long
    assert(tok != NULL);
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

void maps_update_cM_from_map(maps_t *self, maps_t *ref_map) {
  size_t ind = 0;
  size_t nmemb = self->nmemb;
  double *cM = self->cM;
  size_t *bp = self->bp;

  size_t ind_ref = 0;
  size_t nmemb_ref = ref_map->nmemb;
  double *cM_ref = ref_map->cM;
  size_t *bp_ref = ref_map->bp;

  double left_ext_slope = (cM_ref[1] - cM_ref[0]) / (bp_ref[1] - bp_ref[0]);
  double right_ext_slope = (cM_ref[nmemb_ref - 1] - cM_ref[nmemb_ref - 2]) /
                           (bp_ref[nmemb_ref - 1] - bp_ref[nmemb_ref - 2]);

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
      assert(x >= x1 && x <= x2);
      double slope = (y2 - y1) / (x2 - x1);
      double y = y1 + (x - x1) * slope;
      cM[ind] = y;
    }
  }
}

int test_maps_read_map_file() {
  maps_t maps;
  maps_read_map_file(&maps, "./plink.chr18.GRCh38.map");
  for (size_t i = 0; i < maps.nmemb; i++)
    printf("%ld\t%lf\n", maps.bp[i], maps.cM[i]);
  return 0;
}

int test_maps_update_cM_from_map() {
  maps_t maps, refmaps;
  maps_read_uniq_positions(&maps, "./positions_sorted.txt");
  maps_read_map_file(&refmaps, "./plink.chr18.GRCh38.map");
  maps_update_cM_from_map(&maps, &refmaps);
  for (size_t i = 0; i < maps.nmemb; i++)
    printf("%ld\t%lf\n", maps.bp[i], maps.cM[i]);
  return 0;
}

typedef struct {
  int id1;
  int id2;
  size_t start;
  size_t end;
} ibdrec_t;

typedef struct {
  ibdrec_t *lines;
  size_t nmemb;
  size_t ncap;
  size_t inc;

  int chr;
  samples_t samples;
  maps_t maps;
} ibd_t;

int ibdrec_cmp(const void *r1, const void *r2) {
  ibdrec_t *p1 = (ibdrec_t *)r1, *p2 = (ibdrec_t *)r2;
  long long result;
  if ((result = p1->id1 - p2->id1) != 0)
    return result > 0 ? 1 : -1;
  else if ((result = p1->id2 - p2->id2) != 0)
    return result > 0 ? 1 : -1;
  else if ((result = p1->start - p2->start) != 0)
    return result > 0 ? 1 : -1;
  else if ((result = p1->end - p2->end) != 0)
    return result > 0 ? 1 : -1;
  else
    return 0;
}

void ibd_read(ibd_t *self, const char *fn_samples, const char *fn_positions,
              const char *fn_refmap, const char *fn_ibds)

{
  samples_t samples;
  maps_t maps, refmaps;

  ibdrec_t *lines = calloc(10000, sizeof(*lines));
  size_t nmemb = 0;
  size_t ncap = 10000;
  size_t inc = 10000;
  int chr = -1;

  char buff[1000], *p = buff, *tok = NULL, *pName1 = NULL, *pName2 = NULL;
  size_t buff_size = 1000;

  // get sample map information
  samples_read(&samples, fn_samples);
  maps_read_uniq_positions(&maps, fn_positions);
  maps_read_map_file(&refmaps, fn_refmap);
  maps_update_cM_from_map(&maps, &refmaps);

  // for(size_t i=1; i< maps.nmemb ; i++)
	//	  assert(maps.cM[i-1]<=maps.cM[i]);

  FILE *fp = NULL;
  if (strcmp(fn_ibds, "stdin") == 0)
    fp = stdin;
  else
    fp = fopen(fn_ibds, "r");

  int id1, id2, temp;
  while (getline(&p, &buff_size, fp) > 0 && *p != '\n') {
    tok = strtok(p, "\t "); // id1
    // convert name to index
    pName1 = (char *)bsearch(tok, samples.names, samples.nmemb,
                             samples.name_size, str_cmp);
    if (pName1 == NULL)
      fprintf(stderr, "|%s|\n", tok);
    assert(pName1 != NULL);
    id1 = (pName1 - samples.names) / samples.name_size;

    tok = strtok(NULL, "\t "); // hap1, ignore

    tok = strtok(NULL, "\t "); // id2
    pName2 = (char *)bsearch(tok, samples.names, samples.nmemb,
                             samples.name_size, str_cmp);
    if (pName2 == NULL)
      fprintf(stderr, "|%s|\n", tok);
    assert(pName2 != NULL);
    id2 = (pName2 - samples.names) / samples.name_size;
    // printf("id1=%d, id2=%d\n", id1, id2);

    tok = strtok(NULL, "\t "); // hap2, ignore

    tok = strtok(NULL, "\t "); // chr
    if (nmemb == 0)
      chr = atoi(tok);
    else
      assert(chr == atoi(tok));

    tok = strtok(NULL, "\t "); // start
    lines[nmemb].start = strtol(tok, NULL, 10);

    tok = strtok(NULL, "\t "); // end
    lines[nmemb].end = strtol(tok, NULL, 10);

    // swap id1 and id2
    assert(id1 != id2);
    if (id1 == id2) {
      fprintf(stderr, "%d, %d\n", id1, id2);
    }
    if (id1 > id2) {
      temp = id1; // stupid bug here fixed
      id1 = id2;
      id2 = temp;
    }
    lines[nmemb].id1 = id1;
    lines[nmemb].id2 = id2;

    /*
    fprintf(stderr, "%s, %d, %ld, \n %s, %d, %ld, \n%ld, %ld\nsample name %d\n",
                    pName1, id1, pName1 - samples.names,pName2, id2,
    pName2-samples.names, lines[nmemb].start, lines[nmemb].end,
            samples.name_size);
            */

    nmemb++;

    if (nmemb + 10 > ncap) {
      ncap += inc;
      lines = realloc(lines, sizeof(*lines) * ncap);
      inc = ncap;
    }
  }
  if (stdin != fp) {
    fclose(fp);
    fp = NULL;
  }

  // update self
  self->lines = lines;
  lines = NULL;
  self->nmemb = nmemb;
  self->ncap = ncap;
  self->inc = inc;
  self->chr = chr;
  self->samples = samples;
  samples.names = NULL;
  self->maps = maps;
  maps.cM = NULL;
}

void ibd_merge(ibd_t *self) {
  ibdrec_t *lines = self->lines;
  size_t nmemb = self->nmemb;
  maps_t *maps = &self->maps;
  int need_merge = 0;

  // sort
  qsort(lines, nmemb, sizeof(*lines), ibdrec_cmp);

  // merge
  for (size_t i = 1; i < nmemb; i++) {
    size_t last_end = lines[i - 1].end;
    size_t this_beg = lines[i].start;
    // if ids are not equal, no merge
    if (lines[i - 1].id1 != lines[i].id1 || lines[i - 1].id2 != lines[i].id2)
      continue;

    // overlapping
    if (last_end >= this_beg)
      need_merge = 1;
    // close enough
    //  There could be a problem for true ibd. True ibd coordinate is from tree interval length
    //  However, the vcf position is calculated from the scaled pos (0,1) so the position of
    //  true ibd does not match to that found in the vcf most likely due to precision issue.
    else {
      size_t *p1_bp = (size_t *)bsearch(&last_end, maps->bp, maps->nmemb,
                                        sizeof(*maps->bp), size_cmp);
      assert(p1_bp != NULL);
      size_t *p2_bp = (size_t *)bsearch(&this_beg, maps->bp, maps->nmemb,
                                        sizeof(*maps->bp), size_cmp);
      assert(p2_bp != NULL);
      /*NOTE: these are the mering criteria */
      if (p2_bp - p1_bp <= 1 &&
          maps->cM[p2_bp - maps->bp] - maps->cM[p1_bp - maps->bp] <= 0.6)
        need_merge = 1;
      else
        need_merge = 0;
    }

    // merge
    if (need_merge == 1) {
      /*
    fprintf(stderr, "prev\t%d\t%d\t%ld\t%ld\n", lines[i-1].id1, lines[i-1].id2,
    lines[i-1].start, lines[i-1].end); fprintf(stderr,
    "this\t%d\t%d\t%ld\t%ld\n", lines[i].id1, lines[i].id2, lines[i].start,
    lines[i].end);
    */
      size_t this_end = lines[i].end;
      // label the previous id1 = -1
      lines[i - 1].id1 = -1;
      // extend the current record;
      lines[i].start = lines[i - 1].start;
      if (last_end > this_end)
        lines[i].end = last_end;
      /*
      fprintf(stderr, "MERGE\t%d\t%d\t%ld\t%ld\n", lines[i].id1, lines[i].id2,
      lines[i].start, lines[i].end);
      */
    }
  }
}

void ibd_print_merged_ibd(ibd_t *self) {
  char *names = self->samples.names;
  size_t name_size = self->samples.name_size;

  ibdrec_t *lines = self->lines;
  size_t nmemb = self->nmemb;

  maps_t *maps = &self->maps;

  for (size_t i = 0; i < nmemb; i++) {

    if(lines[i].id1 < 0) continue; // for those marked garbage after merging

    char *pName1, *pName2;
    size_t *pbp = NULL;
    double cM_start, cM_end;

    pName1 = names + name_size * lines[i].id1;

    pName2 = names + name_size * lines[i].id2;
    pbp = (size_t *)bsearch(&(lines[i].start), maps->bp, maps->nmemb,
                            sizeof(size_t), size_cmp);
    assert(pbp!=NULL);
    cM_start = maps->cM[pbp - maps->bp];

    pbp = (size_t *)bsearch(&(lines[i].end), maps->bp, maps->nmemb,
                            sizeof(size_t), size_cmp);
    assert(pbp!=NULL);
    cM_end = maps->cM[pbp - maps->bp];

    if(cM_end < cM_start)
    {
    	fprintf(stderr, "%s\t1\t%s\t1\t%d\t%ld\t%ld\t%lf\t%lf\n", pName1, pName2,
            self->chr, lines[i].start, lines[i].end, cM_start, cM_end );
    }

    fprintf(stdout, "%s\t1\t%s\t1\t%d\t%ld\t%ld\t%d\t%lf\n", pName1, pName2,
            self->chr, lines[i].start, lines[i].end, 3, cM_end - cM_start);
  }
}

int test_ibd_read_from_stdin() {
  ibd_t ibd;
  ibd_read(&ibd, "./samples_sorted.txt", "./positions_sorted.txt",
           "./plink.chr18.GRCh38.map", "stdin");
  ibd_merge(&ibd);
  ibd_print_merged_ibd(&ibd);
  return 0;
}

int test_ibd_read_from_file() {
  ibd_t ibd;
  ibd_read(&ibd, "./samples_sorted.txt", "./positions_sorted.txt",
           "./plink.chr18.GRCh38.map", "./hapibd_18.ibd");
  ibd_merge(&ibd);
  ibd_print_merged_ibd(&ibd);
  return 0;
}
int main(int argc, char *argv[]) {
  // test_maps_read_uniq_positions();
  // test_maps_read_map_file();
  // test_maps_update_cM_from_map();
  // test_ibd_read_from_file();

  char *usage = " zcat xxx.ibd.gz | ./ibdmerg <sample_list.txt> "
                "<bp_positions_list> <plink.map>";
  if (argc < 4) {
    fprintf(stderr, "\n%s\n", usage);
    exit(-1);
  }

  
  ibd_t ibd;
  ibd_read(&ibd, argv[1], argv[2], argv[3], "stdin");
  ibd_merge(&ibd);
  ibd_print_merged_ibd(&ibd);

  return 0;
}
