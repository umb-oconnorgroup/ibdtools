#include <algorithm>
#include <argp.h>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

using namespace std;

// -----------------------------------------
// For argument parsing

const char *argp_program_version = "ibdtotal 0.1";
const char *argp_program_bug_address = "gbinux@gmail.com";
static char doc[] =
    "Summarize sample-pair total ibd and save it to a binary matrix.\n";
static char args_doc[] = "SAMPLE_FILE OUT_FILE";
static struct argp_option options[] = {
    {0, 0, 0, 0, "Required arguments:"},
    {"SAMPLE_FILE", 1, 0, OPTION_DOC | OPTION_NO_USAGE,
     "sample file that list all samples in IBD files, one line per sample"},
    {"OUT_FILE", 1, 0, OPTION_DOC | OPTION_NO_USAGE,
     "binary file to save the lower trigular total IBD matrix"},

    {0, 0, 0, 0,
     "Optional arguments: at least one of -i or -m should be specified"},

    {"IBD_file", 'i', "FILE", 0,
     "IBD file with a format of  5 columns delimited by tabs. Columns are sn1, "
     "sn2, start_bp, end_bp, length_in_cM"},
    {"Matrix_files", 'm', "FILE [FILE...]", 0,
     "Existing Ibd matrix file to combine; can be used to merge with partial "
     "total IBD matrix calculated from other chromosomes. ATTENTION: -m option "
     "should be specified at the END!!"},
    //{"verbose", 'v', 0, 0, "Show processing updates, default false"},
    {0}};
struct arguments {
  const char *ibd_fn;
  const char *sample_fn;
  const char *out_fn;
  vector<const char *> binary_in_files;
  bool verbose;
};
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct arguments *arguments = (struct arguments *)state->input;
  switch (key) {
  case 'm':
    arguments->binary_in_files.push_back(arg);
    while (state->next < state->argc &&
           strncmp(state->argv[state->next], "-", 1) != 0) {
      arguments->binary_in_files.push_back(state->argv[state->next]);
      state->next++;
    }
    break;
  case 'i':
    arguments->ibd_fn = arg;
    break;
  case 'v':
    arguments->verbose = true;
    break;
  case ARGP_KEY_ARG:
    switch (state->arg_num) {
    case 0:
      arguments->sample_fn = arg;
      break;
    case 1:
      arguments->out_fn = arg;
      break;
    default:
      argp_usage(state);
    }
    break;
  case ARGP_KEY_END:
    if (state->arg_num < 2)
      argp_usage(state);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

class IbdTotal {
  // data
  vector<uint16_t> cM; // matrix of cM = (uint16_t)(cM * 10)
  unordered_map<string, uint32_t> map_sample_name_2_id;
  vector<string> vec_sample_name;

  // Files
  struct arguments args;

public:
  IbdTotal(struct arguments args) : args(args) {}

  void read_sample_names() {
    // 1. Read sample names to a vector; sort them and
    // 2. Add them to an unordered map to map sample name to sample id;
    // 3. Determine the size of the cM matrix and intialize all elements to 0;
    if (args.sample_fn == NULL)
      return;
    char buff[1000];
    const int size = 1000;
    size_t nsam;
    gzFile gzf_sample = gzopen(args.sample_fn, "r");
    assert(gzf_sample != Z_NULL);

    while (gzgets(gzf_sample, buff, size) != NULL) {
      // remove trailing '\n'
      char *p = buff;
      p = next_field(p);
      vec_sample_name.push_back(buff);
    }
    gzclose(gzf_sample);
    sort(vec_sample_name.begin(), vec_sample_name.end());
    for (size_t i = 0; i < vec_sample_name.size(); i++) {
      map_sample_name_2_id[vec_sample_name[i]] = i;
    }
    assert(vec_sample_name.size() > 0);
    nsam = vec_sample_name.size();
    cM.resize(nsam * (nsam - 1) / 2, 0);
    cerr << "num_sample: " << nsam << '\n';
  }
  void print_map_sample2id() {
    for (auto &e : map_sample_name_2_id)
      cerr << "MAP: " << e.first << ": " << e.second << "\n";
  }

  void read_partial_sum_from_binary_file() {
    // Load the paritial sum from previous runs. Allows merging the total from
    // this chromosome to sum of previous chromosomes.
    if (args.binary_in_files.size() == 0)
      return;
    vector<uint16_t> vec_temp(cM.size());
    for (auto &fn : args.binary_in_files) {
      FILE *fp = fopen(fn, "rb");
      assert(fp != NULL && "binary_in_file can't be opened!");
      size_t num_element_read =
          fread(&vec_temp[0], sizeof(vec_temp[0]), vec_temp.size(), fp);
      assert(num_element_read == cM.size());
      cout << "Before append: cM[0] = " << cM[0] << '\n';
      transform(cM.begin(), cM.end(), vec_temp.begin(), cM.begin(),
                plus<uint16_t>());
      cout << "After append: cM[0] = " << cM[0] << '\n';
      fclose(fp);
    }
  }

  void read_ibd_records() {
    // read ibd_records and update the matrix by adding the length to elements
    // of cM matrix

    if (args.ibd_fn == NULL)
      return;

    size_t rec_counter = 0;
    struct ibd_rect_t {
      char buff[1000];
      size_t buff_size;
      const char *sample1;
      const char *sample2;
      uint16_t cM_x10_int;
    } rec[2] = {0};
    rec[0].buff_size = 1000;
    rec[1].buff_size = 1000;

    ibd_rect_t *p_last_rec = rec + 0, *p_curr_rec = rec + 1, *p_tmp_rec = NULL;
    uint32_t sample_id1;
    uint32_t sample_id2;
    gzFile gzf_ibd = gzopen(args.ibd_fn, "r");
    assert(gzf_ibd != Z_NULL);
    while (gzgets(gzf_ibd, p_curr_rec->buff, p_curr_rec->buff_size) != NULL) {
      // cerr << p_curr_rec->buff;
      char *p = p_curr_rec->buff;
      p_curr_rec->sample1 = p; // col1

      p = next_field(p);
      p_curr_rec->sample2 = p; // col2

      // skip two fields       // col3, col4
      p = next_field(p);
      p = next_field(p);

      p = next_field(p);
      p_curr_rec->cM_x10_int = (uint16_t)(strtod(p, NULL) * 10); // col5
      if (p_last_rec->sample1 == NULL ||
          strcmp(p_curr_rec->sample1, p_last_rec->sample1) != 0 ||
          strcmp(p_curr_rec->sample2, p_last_rec->sample2) != 0) {
        sample_id1 = map_sample_name_2_id.at(p_curr_rec->sample1);
        sample_id2 = map_sample_name_2_id.at(p_curr_rec->sample2);
      }
      // add length to correct matrix element
      cM_at(sample_id1, sample_id2) += p_curr_rec->cM_x10_int;

      // swap p_curr_rec with p_last_rec
      p_tmp_rec = p_curr_rec;
      p_curr_rec = p_last_rec;
      p_last_rec = p_tmp_rec;

      // counter
      rec_counter++;
      if (rec_counter % 1000000 == 0)
        cerr << "\r Num of IBD records read: " << rec_counter << '\n';
    }
    gzclose(gzf_ibd);
  }

  void write_result() {
    assert(args.out_fn != NULL);
    FILE *fp_binary_out = fopen(args.out_fn, "wb");
    assert(fp_binary_out != NULL);
    size_t num_element_written =
        fwrite(&cM[0], sizeof(cM[0]), cM.size(), fp_binary_out);
    assert(num_element_written == cM.size());
    fclose(fp_binary_out);

    // print_cM();
    cerr << "Number of Item written: " << num_element_written << '\n';
  }

private:
  // Helper functions ------------

  uint16_t &cM_at(uint32_t sample_id1, uint32_t sample_id2) {
    // return the element reference from sample ids
    // ...
    // Matrix layout:
    // * - - - - -
    // 0 * - - - -
    // 0 0 * - - -
    // 0 0 0 * - -
    // 0 0 0 0 * -
    // 0 0 0 0 0 *

    size_t row = (sample_id1 > sample_id2) ? sample_id1 : sample_id2;
    size_t col = (sample_id1 < sample_id2) ? sample_id1 : sample_id2;
    size_t index = (row - 1) * row / 2 + col;
    // printf(",%lu", index);
    return cM[index];
  }

  char *next_field(char *p) {
    // Find next tab/newline and change it to '\0' and move to the next char.
    // If find '\0'; return NULL;
    assert(p != NULL && "Cant call next_field with NULL pointer");
    while (*p != '\t' && *p != '\n' && *p != '\0')
      p++;
    if (*p == '\0')
      return NULL;
    *p = '\0';
    p++;
    return p;
  }

  void print_cM() {
    for (size_t i = 0; i < vec_sample_name.size(); i++) {
      for (size_t j = 0; j < vec_sample_name.size(); j++) {
        if (j >= i)
          printf("%d\t", 0);
        else
          printf("%hu\t", cM_at(i, j));
      }
      printf("\n");
    }
  }
};

int main(int argc, char *argv[]) {
  static arguments args;
  args.ibd_fn = NULL;
  args.sample_fn = NULL;
  args.out_fn = NULL;
  args.binary_in_files.resize(0);
  args.verbose = false;

  argp_parse(&argp, argc, argv, 0, 0, &args);

  if (args.binary_in_files.size() == 0 && args.ibd_fn == NULL) {
    cerr << "ERROR: Neigher IBD file or Matrix Files are specified! \n";
    abort();
  }

  IbdTotal tot(args);
  tot.read_sample_names();
  tot.read_partial_sum_from_binary_file();
  tot.read_ibd_records();
  tot.write_result();
  return 0;
}
