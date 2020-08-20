/*
NWD112649	1	NWD278543	1	21	16453940
16634657	2.06	2.059
NWD112649	2	NWD784564	1	21	16453940
16636429	2.06	2.06
NWD112649	1	NWD777183	2	21	16453940
16667976	2.08	2.077
NWD112649	1	NWD715103	2	21	16453940
16676880	2.08	2.083
NWD112649	1	NWD843745	1	21	16453940
16634657	2.06	2.059
NWD112649	0	NWD741447	0	21	16453940
16640949	2.07	2.068
NWD112649	1	NWD778585	1	21	16453940
16634657	2.06	2.059
NWD112649	1	NWD261612	2	21	16453940
16636429	2.06	2.06
NWD112649	1	NWD555559	1	21	16453940
16638138	2.06	2.063
NWD112649	1	NWD463256	1	21	16453940
16667976	2.08	2.077
*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// those segments are non-overlapping
struct IBDTOTAL {
  double *cM;    // array of size equal: num_pairs = n*(n-1)/2
  char *samples; // size of n
  int num_smaple;
  size_t sample_name_len;
};

void IBDTOTAL_writefile(struct IBDTOTAL *ibdtot, const char *filename) {
  FILE *fp = fopen(filename, "wb");
  assert(fp != NULL);
  int num_pairs = ibdtot->num_smaple * (ibdtot->num_smaple - 1) / 2;
  int num_chars = ibdtot->num_smaple * ibdtot->sample_name_len;
  fwrite(ibdtot, sizeof(struct IBDTOTAL), 1, fp);
  fwrite(ibdtot->cM, sizeof(double), num_pairs, fp);
  fwrite(ibdtot->samples, sizeof(char), num_chars, fp);
  fclose(fp);
  fp = NULL;
}

// out: ibdtot , in: filename
void IBDTOTAL_readfile(struct IBDTOTAL *ibdtot, const char *filename) {
  // read ibdtotal structure for size information
  FILE *fp = fopen(filename, "rb");
  assert(fp != NULL);
  fread(ibdtot, sizeof(struct IBDTOTAL), 1, fp);
  ibdtot->cM = NULL;
  ibdtot->samples = NULL;

  // alloc according to size information
  int num_pairs = ibdtot->num_smaple * (ibdtot->num_smaple - 1) / 2;
  int num_chars = ibdtot->num_smaple * ibdtot->sample_name_len;
  ibdtot->cM = (double *)malloc(sizeof(double) * num_pairs);
  assert(ibdtot->cM != NULL);
  ibdtot->samples = (char *)malloc(sizeof(char) * num_chars);
  assert(ibdtot->samples != NULL);

  // read in arrays according to size info
  fread(ibdtot->cM, sizeof(double), num_pairs, fp);
  fread(ibdtot->samples, sizeof(char), num_chars, fp);
  fclose(fp);
  fp = NULL;
}

int cmp_string(const void *str1, const void *str2) {
  return strcmp((char *)str1, (char *)str2);
}

int cmp_double(const void *d1, const void *d2) {
  return *((double *)d1) > *((double *)d2) ? 1 : -1;
}

void IBDTOTAL_init(struct IBDTOTAL *ibdtot, const char *samples_file_name) {
  /* read in sample names */
  FILE *fp = fopen(samples_file_name, "r");
  assert(fp != NULL);
  char buff[100];
  char *p = buff;
  size_t buff_size = 99;

  // determine max sample name len and how many samples;
  ibdtot->sample_name_len = 0;
  ibdtot->num_smaple = 0;
  ibdtot->sample_name_len = 0;
  while (getline(&p, &buff_size, fp) > 0) {
    ibdtot->num_smaple += 1;
    if (strlen(p) + 1 > ibdtot->sample_name_len)
      ibdtot->sample_name_len = strlen(p) + 1;
    buff_size = 99;
  }
  ibdtot->samples = (char *)calloc(ibdtot->num_smaple,
                                   sizeof(char) * ibdtot->sample_name_len);
  // printf("size(p): %ld\n", ibdtot->sample_name_len);
  // printf("size(p): %d\n", ibdtot->num_smaple);

  // read in the sample sames
  fseek(fp, 0, SEEK_SET);
  for (int i = 0; i < ibdtot->num_smaple; i++) {
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
  qsort(ibdtot->samples, ibdtot->num_smaple, ibdtot->sample_name_len,
        cmp_string);

  // allocate mem for cM
  size_t num_pairs = ibdtot->num_smaple * (ibdtot->num_smaple - 1) / 2;
  ibdtot->cM = (double *)calloc(num_pairs, sizeof(double));
};

void IBDTOTAL_free(struct IBDTOTAL *ibdtot, const char *samples_file_name) {
  if (ibdtot->cM != NULL) {
    free(ibdtot->cM);
    ibdtot->cM = NULL;
  }
  if (ibdtot->samples != NULL) {
    free(ibdtot->samples);
    ibdtot->samples = NULL;
  }
}
int IBDTOTAL_find_smaple_intID(struct IBDTOTAL *ibdtot, const char *str) {
  int i;
  char *p = (char *)bsearch(str, ibdtot->samples, ibdtot->num_smaple,
                            ibdtot->sample_name_len, cmp_string);
  i = (p - ibdtot->samples) / ibdtot->sample_name_len;
  return i;
}

int IBDTOTAL_sampleIds_to_index(struct IBDTOTAL *ibdtot, int id1, int id2) {

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

  return (2 * ibdtot->num_smaple - id1 - 1) * id1 / 2 + id2 - id1 - 1;
}

void IBDTOTAL_add_cM(struct IBDTOTAL *ibdtot, const char *str1,
                     const char *str2, double cM) {
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

double IBDTOTAL_get_cM(struct IBDTOTAL *ibdtot, int id1, int id2) {
  int idtemp;
  if (id1 > id2) {
    idtemp = id2;
    id2 = id1;
    id1 = idtemp;
  }

  int index = IBDTOTAL_sampleIds_to_index(ibdtot, id1, id2);
  return ibdtot->cM[index];
}

void IBDTOTAL_print(struct IBDTOTAL *ibdtot) {
  /*
  for(int i =0; i <ibdtot->num_smaple; i++)
  {
          printf("%s\n", ibdtot->samples + i * ibdtot->sample_name_len);
  }
  */

  // printf("-----------\n");

  double cM;
  double genome_size = 3545.8;
  for (int i = 0; i < ibdtot->num_smaple - 1; i++)
    for (int j = i + 1; j < ibdtot->num_smaple; j++) {
      cM = IBDTOTAL_get_cM(ibdtot, i, j);
      if (cM < 1.9)
        continue;
      // printf("%d\t%s\t%d\t%s\t%lf\n",
      printf("%s\t%s\t%g\t%g\n",
             // i,
             ibdtot->samples + i * ibdtot->sample_name_len,
             // j,
             ibdtot->samples + j * ibdtot->sample_name_len, 
	     cM, cM / genome_size);
    }
}

/*
 *   https://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries
 */
int WriteIMAGE(const char *filename, int w, int h, unsigned char *img) {
  FILE *f;
  int filesize =
      54 + 3 * w * h; // w is your image width, h is image height, both int

  unsigned char bmpfileheader[14] = {'B', 'M', 0, 0,  0, 0, 0,
                                     0,   0,   0, 54, 0, 0, 0};
  unsigned char bmpinfoheader[40] = {40, 0, 0, 0, 0, 0, 0,  0,
                                     0,  0, 0, 0, 1, 0, 24, 0};
  unsigned char bmppad[3] = {0, 0, 0};

  bmpfileheader[2] = (unsigned char)(filesize);
  bmpfileheader[3] = (unsigned char)(filesize >> 8);
  bmpfileheader[4] = (unsigned char)(filesize >> 16);
  bmpfileheader[5] = (unsigned char)(filesize >> 24);

  bmpinfoheader[4] = (unsigned char)(w);
  bmpinfoheader[5] = (unsigned char)(w >> 8);
  bmpinfoheader[6] = (unsigned char)(w >> 16);
  bmpinfoheader[7] = (unsigned char)(w >> 24);
  bmpinfoheader[8] = (unsigned char)(h);
  bmpinfoheader[9] = (unsigned char)(h >> 8);
  bmpinfoheader[10] = (unsigned char)(h >> 16);
  bmpinfoheader[11] = (unsigned char)(h >> 24);

  f = fopen(filename, "wb");
  fwrite(bmpfileheader, 1, 14, f);
  fwrite(bmpinfoheader, 1, 40, f);
  for (int i = 0; i < h; i++) {
    fwrite(img + (w * (h - i - 1) * 3), 3, w, f);
    fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
  }

  free(img);
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
double IBDTOTAL_get_cM_percentile(struct IBDTOTAL *ibdtot, double percentage) {
  double percentile = 0, *cM_temp = NULL;
  size_t index = 0;
  size_t num_zeros = 0, i, num_nonzeros;
  size_t num_pairs = ibdtot->num_smaple * (ibdtot->num_smaple - 1) / 2;

  assert(percentage >= 0);
  assert(percentage <= 100);

  cM_temp = (double *)malloc(sizeof(double) * num_pairs);
  assert(cM_temp != NULL);
  memcpy(cM_temp, ibdtot->cM, sizeof(double) * num_pairs);

  qsort(cM_temp, num_pairs, sizeof(double), cmp_double);

  assert(cM_temp[0]>=0);

  // calculate num_zeros
  for (size_t i = 0; i < num_pairs && cM_temp[i] < 0.1; i++, num_zeros++) {
  };

  num_nonzeros = num_pairs - num_zeros;

  index = (int)(percentage * num_nonzeros / 100);
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
    long index = (long)(num_pairs * perc);
    if (index >= num_pairs)
      index = num_pairs - 1;
    fprintf(fp, "%g\t%g\n", perc, cM_temp[index]);
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

int IBDTOTAL_generate_map(struct IBDTOTAL *ibdtot, char *filename,
                          double cutoff_cM, int *sample_orders) {
  int w = ibdtot->num_smaple / 3;
  int h = w;
  int index, i1, i2;
  double value, temp;
  double values9[9];
  int is_saturated;
  unsigned char *img = (unsigned char *)calloc(w * h, sizeof(*img) * 3);

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
        img[3 * (i * w + j) + 2] = (unsigned char)255; // r
        img[3 * (i * w + j) + 1] = (unsigned char)255; // g
        img[3 * (i * w + j) + 0] = (unsigned char)255; // b
      } else if (is_saturated == 0) {
        img[3 * (i * w + j) + 2] = (unsigned char)value; // r
        img[3 * (i * w + j) + 1] = (unsigned char)value; // g
        img[3 * (i * w + j) + 0] = (unsigned char)value; // b
      } else {
        img[3 * (i * w + j) + 2] = (unsigned char)255; // r
        img[3 * (i * w + j) + 1] = (unsigned char)0;   // g
        img[3 * (i * w + j) + 0] = (unsigned char)0;   // b
      }
    }

  WriteIMAGE(filename, w, h, img);

  return 0;
}

// TODO: work on color pallete
void IBDTOTAL_generate_map_legend(struct IBDTOTAL *ibdtot, const char *image_fn,
                                  int *sample_orders) {
  int num_samples = ibdtot->num_smaple;
  struct RGB {
    unsigned char b;
    unsigned char g;
    unsigned char r;
  };

  struct RGB color_8v1_A_NAD = {176, 121, 44};   // 1148
  struct RGB color_8v1_B_NAD = {54, 127, 251};   // 880
  struct RGB color_C1_hisp = {59, 158, 57};      // 463
  struct RGB color_Emerge_ids = {53, 45, 210};   // 4852
  struct RGB color_LARGEPD = {185, 109, 147};    // 1501
  struct RGB color_NWD = {78, 86, 138};          // 17560
  struct RGB color_PAGEII = {192, 125, 223};     // 5213
  struct RGB color_SIGMA_OMNI = {127, 127, 127}; // 1147

  int width = 5000, height = num_samples;
  unsigned char *img = malloc(3 * width * height * sizeof(*img));

  for (int row = 0; row < num_samples; row++) {
    int index;
    if (sample_orders != NULL)
      index = sample_orders[row];
    else
      index = row;
    struct RGB color = {255, 255, 255};
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

int main(int argc, char *argv[]) {
  struct IBDTOTAL ibdtotal;

  char buff[1000], *p = buff, *token = NULL;
  size_t size = 1000;
  size_t lineCount = 0;
  int colCount = 0;
  char *s1, *s2;
  double cM, cutoff_cM;
  char *beg = NULL, *sep = NULL;

  if (argc == 1 || argc > 4) {
    fprintf(stderr, "\nUsage: ./ibdtotal <samples_list file> [<matrix_file> "
                    "<samples_orders_file>]\n\n");
    exit(1);
  }
  /* read merged ibd from stdin*/
  else if (argc == 2) {
    IBDTOTAL_init(&ibdtotal, argv[1]); // argv[1] = samples_list file

    // parse lines and save it cM, sample pairs into arrays
    while (getline(&p, &size, stdin) > 0) {
      lineCount++;
      // if(lineCount % 10000 == 0) fprintf(stderr, "line: %ld\n", lineCount);
      token = strtok(p, "\t");
      colCount = 0;
      while (colCount <= 8) {
        colCount++;
        if (colCount == 1)
          s1 = token;
        else if (colCount == 3)
          s2 = token;
        else if (colCount == 8)
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
  else if (argc > 3) {
    IBDTOTAL_init(&ibdtotal, argv[1]); // argv[1] = samples_list file
    fprintf(stderr, "reading file\n");
    IBDTOTAL_readfile(&ibdtotal, argv[2]);
    fprintf(stderr, "done reading file\n");

    int *sample_orders = NULL;
    if (argc == 4) {
      // reading order files
      FILE *fp_order = fopen(argv[3], "r"); // argc[3] is sample order files
      assert(fp_order != NULL);
      sample_orders = calloc(ibdtotal.num_smaple, sizeof(*sample_orders));
      assert(sample_orders != NULL);
      for (int i = 0; i < ibdtotal.num_smaple; i++) {
        fscanf(fp_order, "%d", sample_orders + i);
        // fprintf(stderr, "%d " , sample_orders[i]);
      }
      fclose(fp_order);
      fp_order = NULL;
    }

    fprintf(stderr, "start");

    // calculate percentile
    cutoff_cM = IBDTOTAL_get_cM_percentile(&ibdtotal, 90);
    fprintf(stderr, "percentage= %lf, percentile %lf\n", 90.0, cutoff_cM);
    // cutoff_cM = 50;
    // generate bmp
    IBDTOTAL_generate_map(&ibdtotal, "img_raw.bmp", cutoff_cM, NULL);
    if (sample_orders != NULL)
      IBDTOTAL_generate_map(&ibdtotal, "img_clust.bmp", cutoff_cM,
                            sample_orders);
    // generate sample color
    IBDTOTAL_generate_map_legend(&ibdtotal, "img_raw_legend.bmp", NULL);
    if (sample_orders != NULL)
      IBDTOTAL_generate_map_legend(&ibdtotal, "img_clust_legend.bmp",
                                   sample_orders);
  }

  return 0;
}
