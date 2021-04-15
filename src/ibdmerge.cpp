#include <algorithm>
#include <argp.h>
#include <cassert>
#include <cstdio>
#include <execution>
#include <fstream>
#include <htslib/vcf.h>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

#define UINT32T_EMPTY ~(uint32_t)(0)

using namespace std;

// -----------------------------------------
// For argument parsing

const char *argp_program_version = "ibdmerge 0.1";
const char *argp_program_bug_address = "gbinux@gmail.com";
static char doc[] =
    "ibdmerge merges IBD segments following Browning's tool "
    "merge-ibd-segments.jar. This is written to handle large datasets.\n";
static char args_doc[] = "IBD_FILE VCF_FILE MAP_FILE";
static struct argp_option options[] = {
    {0, 0, 0, 0, "Required arguments:"},
    {"IBD_FILE", 1, 0, OPTION_DOC | OPTION_NO_USAGE,
     "IBD files with format `sn1:sn2\tstart\tend...`, sorted by sort -k1,1 "
     "-k2,2n -k3,3n"},
    {"VCF_FILE", 1, 0, OPTION_DOC | OPTION_NO_USAGE,
     "VCF file (tab delimited) used to call IBD. The VCF file is expected to only have "
     "biallelic sites"},
    {"MAP_FILE", 1, 0, OPTION_DOC | OPTION_NO_USAGE,
     "MAP file (space delimited) used to call IBD"},

    {0, 0, 0, 0, "Optional arguments:"},
    {"max_cM", 'm', "float", 0,
     "max length of gap between IBD segments (cM), default 0.6"},
    {"discord", 'd', "int", 0,
     "max number of genotypes in IBD gap that are inconsistent with IBD, "
     "default 1"},
    {"out", 'o', "filename", 0,
     "IBD output file after merging, default stdout"},
    {"verbose", 'v', 0, 0, "Show processing updates, default false"},
    {0}};
struct arguments {
  const char *ibd_fn;
  const char *vcf_fn;
  const char *map_fn;
  const char *out_fn;
  float max_cM;
  int discord;
  bool verbose;
};
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct arguments *arguments = (struct arguments *)state->input;
  switch (key) {
  case 'm':
    arguments->max_cM = atof(arg);
    break;
  case 'd':
    arguments->discord = atoi(arg);
    break;
  case 'o':
    arguments->out_fn = arg;
    break;
  case 'v':
    arguments->verbose = true;
    break;
  case ARGP_KEY_ARG:
    switch (state->arg_num) {
    case 0:
      arguments->ibd_fn = arg;
      break;
    case 1:
      arguments->vcf_fn = arg;
      break;
    case 2:
      arguments->map_fn = arg;
      break;
    default:
      argp_usage(state);
    }
    break;
  case ARGP_KEY_END:
    if (state->arg_num < 3)
      argp_usage(state);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

// -----------------------------------------
// For input, merging and output

class GeneticMap {

public:
  vector<uint32_t> bp;
  vector<float> cm;
  float last_slope;

  void parse_genetic_map(const char *genetic_map_file) {
    string line, field;
    ifstream ifs(genetic_map_file, ios_base::in);
    if (!ifs.is_open()) {
      cerr << "Cannot open map file \n";
      exit(1);
    }
    int line_count = 0;
    bp.push_back(0);
    cm.push_back(0);
    cerr << "Parsing genetic map ...";
    for (; getline(ifs, line); line_count++) {
      istringstream ss(line);
      int field_count = 0;
      for (; getline(ss, field, ' '); // plink map delimiters are space
           field_count++) {
        if (field_count == 2) {
          cm.push_back(stof(field));
        } else if (field_count == 3) {
          bp.push_back(stoul(field));
        } else { // blank
        }
      }
      assert(field_count == 4 && "genetic map does not have 4 columns");
    };
    ifs.close();
    assert(bp.size() >= 2 && "bp vector should be not less than 2");
    last_slope = (cm[cm.size() - 1] - cm[cm.size() - 2]) /
                 (bp[bp.size() - 1] - bp[bp.size() - 2]);

    cerr << "  done!\n";
  }

  float get_cm(uint32_t bp_pos) {
    assert(bp_pos >= 0 && "bp_pos should not be less than 0");
    float cm_pos = 0;
    if (bp_pos >= bp[bp.size() - 1]) {
      cm_pos = (bp_pos - bp[bp.size() - 1]) * last_slope + cm[cm.size() - 1];
    } else {
      vector<uint32_t>::iterator it = upper_bound(bp.begin(), bp.end(), bp_pos);
      size_t idx = distance(bp.begin(), it);
      float slope = (cm[idx] - cm[idx - 1]) / (bp[idx] - bp[idx - 1]);
      cm_pos = (bp_pos - bp[idx - 1]) * slope + cm[idx - 1];
    }
    return cm_pos;
  }
};

class Variants {

public:
  unordered_map<uint32_t, uint32_t>
      pos2id_map; // position mapped to position id
  vector<uint32_t> bp_vec;
  vector<float> cm_vec;

  void set_size(uint32_t npos) {
    pos2id_map.reserve(npos);
    cm_vec.resize(npos, 0);
    bp_vec.resize(npos, 0);
  }

  void add_pos(uint32_t pos) {
    this->bp_vec.push_back(pos);
    uint32_t pos_id = bp_vec.size() - 1;
    this->pos2id_map[pos] = pos_id;
  }

  uint32_t get_pos_id(uint32_t pos) { return pos2id_map.at(pos); }

  float get_cM(uint32_t pos_id) { return cm_vec[pos_id]; }

  void calc_genetic_pos(GeneticMap &gmap) {
    if (cm_vec.size() <= bp_vec.size())
      cm_vec.resize(bp_vec.size());
    transform(bp_vec.begin(), bp_vec.end(), cm_vec.begin(),
              [&](uint32_t pos) -> float { return gmap.get_cm(pos); });
  }

  pair<uint32_t, uint32_t> get_pos_id_between(uint32_t left_bp_pos,
                                              uint32_t right_bp_pos) {
    assert(left_bp_pos <= right_bp_pos);
    auto left = upper_bound(bp_vec.begin(), bp_vec.end(), left_bp_pos);
    auto right = lower_bound(bp_vec.begin(), bp_vec.end(), right_bp_pos);
    return pair<uint32_t, uint32_t>(distance(bp_vec.begin(), left),
                                    distance(bp_vec.begin(), right));
  }
};

class Samples {

public:
  unordered_map<string, uint32_t>
      name2id_map; // sample name mapped to sample id.
  vector<string> name_vec;

  void add_sample(const char *key) {
    name_vec.push_back(key);
    uint32_t id = name_vec.size() - 1;
    name2id_map[name_vec[id]] = id;
  }

  uint32_t get_sample_id(const char *sample_name) {
    return name2id_map[sample_name];
  }
};

class Genotypes {

public:
  vector<uint64_t> genotype_vec;
  size_t size_per_row;
  uint32_t npos, nsam;

  size_t shrink_to_fit() {
    genotype_vec.shrink_to_fit();
    return genotype_vec.size() * sizeof(uint64_t);
  }

  size_t set_size(uint32_t num_variants, uint32_t num_samples, int num_ploidy) {
    size_t total_bytes = 0;
    size_t nsam_per_64bit = 64 / 2 / num_ploidy;

    // Currently, genotype info is expected to take all bits (no gaps)
    assert((64 / 2) % num_ploidy == 0);
    size_per_row = num_samples / nsam_per_64bit;
    if (num_samples % nsam_per_64bit != 0)
      size_per_row++;
    genotype_vec.resize((size_t)size_per_row * (size_t)num_variants, 0);
    nsam = num_samples;
    npos = num_variants;
    total_bytes =
        (size_t)size_per_row * (size_t)num_variants * sizeof(uint64_t);
    return total_bytes;
  }

  int get_allele(uint32_t pos_id, uint32_t sample_id, int allele_id) {
    uint32_t bit_offset = sample_id * 4 + allele_id * 2;
    uint32_t element_idx = size_per_row * pos_id + bit_offset / 64;
    return (genotype_vec[element_idx] >> (bit_offset % 64)) & 3;
  }

  void set_allele(uint32_t pos_id, uint32_t sample_id, int allele_id,
                  uint64_t allele) {
    uint32_t bit_offset = sample_id * 4 + allele_id * 2;
    uint32_t element_idx = size_per_row * pos_id + bit_offset / 64;
    genotype_vec.at(element_idx) |= ((allele & 3) << (bit_offset % 64));
  }

  void print_genotypes() {
    for (uint32_t i = 0; i < npos; i++) {
      for (uint32_t j = 0; j < nsam; j++) {
        cout << get_allele(i, j, 0);
        cout << get_allele(i, j, 1);
        cout << ' ';
      }
      cout << '\n';
    }
  }
};

class Vcf {

public:
  Samples sample;
  Variants variant;
  Genotypes genotype;
  size_t mem_use_estimate;
  bool verbose;

  void parse_vcf(string vcf_file_name) {
    // intialization
    htsFile *fp = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = NULL;
    int32_t res = 0, count = 0, rec_count = 0;
    uint32_t *p_allele = NULL, nsam = 0, ploidy = 0;

    cerr << "Parsing VCF file ...\n";

    // read the file the second time to save genotype data
    fp = bcf_open(vcf_file_name.c_str(), "r");
    assert(fp && "Can't open vcf file");

    if (!(hdr = bcf_hdr_read(fp)))
      assert(hdr && "Can't read header of the vcf file");
    nsam = hdr->n[BCF_DT_SAMPLE];

    rec = bcf_init();
    assert(rec && "Can't init bcf record");

    cerr << "Parsing VCF file ...";

    for (int i = 0; i < hdr->n[BCF_DT_SAMPLE]; i++)
      this->sample.add_sample(hdr->id[BCF_DT_SAMPLE][i].key);

    while ((res = bcf_read(fp, hdr, rec)) == 0) {
      assert(res == 0 && "Can't read record");
      if (verbose) {
        cerr << "\r    Read variants and genotypes: " << rec_count;
        if (rec->errcode)
          cerr << " Errorcode" << rec->errcode;
      }

      res = bcf_get_genotypes(hdr, rec, &p_allele, &count);
      if (ploidy == 0)
        ploidy = count / nsam;
      assert(res != 0 && "Can't get genotypes");
      // set alleles
      mem_use_estimate = genotype.set_size(rec_count + 1, count / 2, ploidy);

      if (verbose) {
        cerr << "    Estimated memory usage: "
             << 1.0 * mem_use_estimate / 1024 / 1024 / 1024 << " Gb";
      }
      for (int32_t i = 0; i < count; i++) {
        this->genotype.set_allele(rec_count, i / 2, i % 2,
                                  p_allele[i] >> 1 & 3);
      }
      // set positions
      this->variant.add_pos(rec->pos +
                            1); // pos is 0-based, while vcf is 1-based. wield
      rec_count++;
    };
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    assert(bcf_close(fp) == 0);
    fp = NULL;
    cerr << "\n    Done parsing VCF file\n";

    // release memory used memory
    mem_use_estimate = genotype.shrink_to_fit();
  }
};

class IbdMerger {
  gzFile gz_ibd_in, gz_ibd_out;
  bool gz_ibd_in_eof;
  size_t max_lines;
  vector<uint32_t> vSid1, vSid2;
  vector<uint32_t> vIndex;
  vector<uint32_t> vStarts, vEnds, vRange;
  Vcf &vcf_ref;
  GeneticMap &gmap_ref;

  // criteria
  int max_disc_sites;
  float max_cM;

  // buffers
  char buffer[1000];
  string field;
  istringstream iss;

  // counters
  uint64_t merged_count, printed_count;

  // print level
  bool verbose;

public:
  IbdMerger(gzFile gzFile_ibd_in, gzFile gzFile_ibd_out, size_t max_lines,
            int max_disc_sites, float max_cM, Vcf &vcf, GeneticMap &gmap,
            bool verbose)
      : gz_ibd_in(gzFile_ibd_in), gz_ibd_out(gzFile_ibd_out),
        max_lines(max_lines), vcf_ref(vcf), gmap_ref(gmap),
        max_disc_sites(max_disc_sites), max_cM(max_cM), verbose(verbose) {
    vSid1.reserve(max_lines);
    vSid2.reserve(max_lines);
    vStarts.reserve(max_lines);
    vEnds.reserve(max_lines);
    vRange.reserve(max_lines);
    field.reserve(1000);
    gz_ibd_in_eof = false;
    merged_count = 0;
    printed_count = 0;
  }

  ~IbdMerger() {}

  bool is_ibd_in_eof() { return gz_ibd_in_eof; }

  void read_sorted_ibd_into_buffer() {
    char *gzret = NULL;
  uint32_t sample_id1, sample_id2;
    while (vStarts.size() < max_lines &&
           (gzret = gzgets(gz_ibd_in, buffer, 1000)) != Z_NULL) {
      iss.clear();
      iss.str(buffer);
      getline(iss, field, ':');
      sample_id1 = vcf_ref.sample.name2id_map.at(field);
      getline(iss, field, '\t');
      sample_id2 = vcf_ref.sample.name2id_map.at(field);
      // sample names are only saved when there is a new pair of samples,
      // and the vIndex saves the starting index for each sample pair
      if (vSid1.size() == 0 || *(vSid1.end() - 1) != sample_id1 ||
          *(vSid2.end() - 1) != sample_id2) {
        vSid1.push_back(sample_id1);
        vSid2.push_back(sample_id2);
        vIndex.push_back(vStarts.size());
      }
      // bp positions are saved for each record.
      getline(iss, field, '\t');
      vStarts.push_back(stoul(field));
      getline(iss, field, '\t');
      vEnds.push_back(stoul(field));
    };

    if (gzeof(gz_ibd_in))
      gz_ibd_in_eof = true;
    else if (gzret == Z_NULL) {
      cerr << "Error on gzgets!\n";
      exit(1);
    }
  }

  void write_merged_ibd_to_file() {
    // last merged sample pair_id
    uint32_t num_sam_pair, sp_id, num_rec, end_id, end_rec_id, rec_id;
    float cM;
    num_sam_pair = vIndex.size();
    num_rec = vStarts.size();
    end_id = gz_ibd_in_eof ? num_sam_pair : (num_sam_pair - 1);
    for (sp_id = 0; sp_id < end_id; sp_id++) {
      string &sn1 = vcf_ref.sample.name_vec.at(vSid1[sp_id]);
      string &sn2 = vcf_ref.sample.name_vec.at(vSid2[sp_id]);
      end_rec_id = (sp_id == num_sam_pair - 1) ? num_rec : vIndex[sp_id + 1];
      for (rec_id = vIndex[sp_id]; rec_id < end_rec_id; rec_id++) {
        uint32_t &s_pos = vStarts[rec_id];
        uint32_t &e_pos = vEnds[rec_id];
        if (s_pos != UINT32T_EMPTY) {
          cM = gmap_ref.get_cm(e_pos) - gmap_ref.get_cm(s_pos);
          if (gz_ibd_out == Z_NULL)
            cout << sn1 << '\t' << sn2 << '\t' << s_pos << '\t' << e_pos << '\t'
                 << cM << '\n';
          else
            gzprintf(gz_ibd_out, "%s\t%s\t%lu\t%lu\t%0.5g\n", sn1.c_str(),
                     sn2.c_str(), s_pos, e_pos, cM);
          printed_count++;
        } else {
          merged_count++;
        }
      }
    }
    // verbose output to cerr
    if (verbose)
      cerr << "\r    Merged segment: " << merged_count
           << "; printed segment: " << printed_count;
    // recycle un-processed records.
    if (!gz_ibd_in_eof) {
      // only keep the last elements
      vSid1.erase(vSid1.begin(), vSid1.end() - 1);
      vSid2.erase(vSid2.begin(), vSid2.end() - 1);
      // only keep records for last sample pair
      vStarts.erase(vStarts.begin(),
                    vStarts.begin() + vIndex[num_sam_pair - 1]);
      vEnds.erase(vEnds.begin(), vEnds.begin() + vIndex[num_sam_pair - 1]);
      // only one index point to 0
      vIndex.resize(1, 0);
    }
  }

  void print_counters() {
    cerr << "\n    Merged segment: " << merged_count
         << "; printed segment: " << printed_count;
  }

  bool is_discordant_excessive(uint32_t sample_id1, uint32_t sample_id2,
                               uint32_t left_bp_pos, uint32_t right_bp_pos) {
    vector<uint32_t> &bp_vec = vcf_ref.variant.bp_vec;
    Genotypes &genotype = vcf_ref.genotype;

    // cal position ids between (left_bp_pos, right_bp_pos)
    assert(left_bp_pos <= right_bp_pos);
    auto left = upper_bound(bp_vec.begin(), bp_vec.end(), left_bp_pos);
    auto right = lower_bound(bp_vec.begin(), bp_vec.end(), right_bp_pos);
    auto first_pos_id = distance(bp_vec.begin(), left);
    auto last_pos_id = distance(bp_vec.begin(), right);

    int n_discordant = 0;
    for (uint32_t pos_id = first_pos_id; pos_id < last_pos_id; pos_id++) {
      int a0, a1, b0, b1;
      a0 = genotype.get_allele(pos_id, sample_id1, 0);
      a1 = genotype.get_allele(pos_id, sample_id1, 1);
      b0 = genotype.get_allele(pos_id, sample_id2, 0);
      b1 = genotype.get_allele(pos_id, sample_id2, 1);
      if (a0 == a1 && b0 == b1 && a1 != b1 && a1 != 0 && b1 != 0) {
        n_discordant++;
        if (n_discordant > max_disc_sites)
          return true;
      }
    }
    return false;
  }

  void merge_for_sample_pair(uint32_t i) {
    uint32_t first, last; // [inclusive, non-includsive)
    uint32_t bp_pos1, bp_pos2;
  uint32_t sample_id1, sample_id2;
    first = vIndex[i];
    last = (i == vIndex.size() - 1) ? vStarts.size() : vIndex[i + 1];
    sample_id1 = vSid1[i];
    sample_id2 = vSid2[i];

    // j point to prev record; k point to curr record;
    for (uint32_t j = first, k = first + 1; k < last; k++) {
      bp_pos1 = vEnds[j];
      bp_pos2 = vStarts[k];
      // if overlapping or non-overlapping but close enough and have discord
      // sites no greater than threshold, then merge
      if (bp_pos1 >= bp_pos2 ||
          ((gmap_ref.get_cm(bp_pos2) - gmap_ref.get_cm(bp_pos1) < max_cM) &&
           !is_discordant_excessive(sample_id1, sample_id2, bp_pos1,
                                    bp_pos2))) {
        // merge k record into j record, mark k's start pos to UINT32T_EMPTY to
        // represent empty record
        vEnds[j] = vEnds[j] > vEnds[k] ? vEnds[j] : vEnds[k];
        vStarts[k] = UINT32T_EMPTY;
      } else {
        // no merge, than use k as the prev record
        j = k;
      }
    }
  }

  void merge_ibd_in_buffer() {
    // Not using the the last element because the groups of ibd belonging to
    // this pair may be completed readin if not at the end of file
    vRange.resize(vIndex.size());
    iota(vRange.begin(), vRange.end(), 0);

    if (!gz_ibd_in_eof)
      for_each(execution::par, vRange.begin(), vRange.end() - 1,
               [this](uint32_t i) { this->merge_for_sample_pair(i); });
    else
      for_each(execution::par, vRange.begin(), vRange.end(),
               [this](uint32_t i) { this->merge_for_sample_pair(i); });
  }
};

int main(int argc, char *argv[]) {
  struct arguments arguments;
  arguments.max_cM = 0.6;
  arguments.discord = 1;
  arguments.out_fn = "stdout";
  arguments.ibd_fn = NULL;
  arguments.map_fn = NULL;
  arguments.vcf_fn = NULL;
  arguments.verbose = false;

  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  ios_base::sync_with_stdio(false);
  GeneticMap gmap;
  gmap.parse_genetic_map(arguments.map_fn);
  size_t mem_use_estimate = 0;

  Vcf vcf;
  vcf.verbose = arguments.verbose;
  vcf.parse_vcf(arguments.vcf_fn);
  vcf.variant.calc_genetic_pos(gmap);

  mem_use_estimate += vcf.mem_use_estimate;

  gzFile fp_in = gzopen(arguments.ibd_fn, "r");
  gzFile fp_out = Z_NULL;
  if (strcmp(arguments.out_fn, "stdout") != 0) {
    fp_out = gzopen(arguments.out_fn, "w");
    gzbuffer(fp_out, 1024 * 1024 * 100);
  }

  FILE *fp_tmp = fopen(arguments.ibd_fn, "r");
  fseek(fp_tmp, 0, SEEK_END);
  size_t file_size = ftell(fp_tmp);
  fclose(fp_tmp);
  size_t max_lines = 100000000;
  if (file_size /20 < max_lines) max_lines = file_size/20;

  IbdMerger ibdmerger(fp_in, fp_out, max_lines, arguments.discord,
                      arguments.max_cM, vcf, gmap, arguments.verbose);
  mem_use_estimate += max_lines * (15 + 15 + 4 + 4) + 1024 * 1024 + 100;
  if (arguments.verbose) {
    cerr << "\n    Estimated memory usage: "
         << 1.0 * mem_use_estimate / 1024 / 1024 / 1024 << "  Gb\n";
  }

  cerr << "Merging IBD segments... \n";
  while (!ibdmerger.is_ibd_in_eof()) {
    ibdmerger.read_sorted_ibd_into_buffer();
    ibdmerger.merge_ibd_in_buffer();
    ibdmerger.write_merged_ibd_to_file();
  }
  ibdmerger.print_counters();
  cerr << "\n    Done merging IBD segments! \n";

  if (fp_in != Z_NULL)
    gzclose(fp_in);
  if (fp_out != Z_NULL)
    gzclose(fp_out);
  return 0;
}
