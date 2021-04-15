#include <algorithm>
#include <argp.h>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <htslib/vcf.h>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

using namespace std;

// -----------------------------------------
// For argument parsing

const char *argp_program_version = "ibdmerge 0.1";
const char *argp_program_bug_address = "gbinux@gmail.com";
static char doc[] =
    "ibdmerge merges IBD segments following Browning's tool "
    "merge-ibd-segments.jar. This is rewritten to handle large datasets.\n";
static char args_doc[] = "IBD_FILE VCF_FILE MAP_FILE";
static struct argp_option options[] = {
    {0, 0, 0, 0, "Required arguments:"},
    {"IBD_FILE", 1, 0, OPTION_DOC | OPTION_NO_USAGE,
     "IBD files with format `sn1:sn2\tstart\tend...`, sorted by sort -k1,1 "
     "-k2,2n -k3,3n"},
    {"VCF_FILE", 1, 0, OPTION_DOC | OPTION_NO_USAGE,
     "VCF file used to call IBD. The VCF file is expected to only have "
     "biallelic sites"},
    {"MAP_FILE", 1, 0, OPTION_DOC | OPTION_NO_USAGE,
     "MAP file used to call IBD"},

    {0, 0, 0, 0, "Optional arguments:"},
    {"max_cM", 'm', "float", 0,
     "max length of gap between IBD segments (cM), default 0.6"},
    {"discord", 'd', "int", 0,
     "max number of genotypes in IBD gap that are inconsistent with IBD, "
     "default 1"},
    {"out", 'o', "filename", 0, "IBD output file after merging, default stdout"},
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

  void add_pos(uint32_t pos, uint32_t pos_id) {
    this->pos2id_map[pos] = pos_id;
    this->bp_vec[pos_id] = pos;
  }

  uint32_t get_pos_id(uint32_t pos) { return pos2id_map.at(pos); }

  float get_cM(uint32_t pos_id) { return cm_vec[pos_id]; }

  void calc_genetic_pos(GeneticMap &gmap) {
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

  void set_size(uint32_t nsam) {
    name2id_map.reserve(nsam);
    name_vec.resize(nsam);
  }

  void add_sample(const char *key, uint32_t id) {
    name_vec[id] = key;
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
  bool verbose;

  void parse_vcf(string vcf_file_name) {
    // intialization
    htsFile *fp = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = NULL;
    int32_t res = 0, count = 0, rec_count = 0;
    uint32_t *p_allele = NULL, npos = 0, nsam = 0;
    size_t estimated_mem = 0;

    // read the file first time to get no. of lines, no of samples.
    cerr << "Parsing VCF file ...\n";
    fp = bcf_open(vcf_file_name.c_str(), "r");
    assert(fp && "Can't open vcf file");

    if (!(hdr = bcf_hdr_read(fp)))
      assert(hdr && "Can't read header of the vcf file");
    nsam = hdr->n[BCF_DT_SAMPLE];

    rec = bcf_init1();
    assert(rec && "Can't init bcf record");

    while ((res = bcf_read(fp, hdr, rec)) == 0) {
      if (verbose)
        cerr << "\r    Count variants: " << npos;
      npos++;
    }
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    assert(bcf_close(fp) == 0);
    fp = NULL;

    cerr << "\n    Detected nsam = " << nsam << ", npos = " << npos;

    // reserve space for samples, positions, genotypes.
    this->sample.set_size(nsam);
    this->variant.set_size(npos);
    estimated_mem = this->genotype.set_size(npos, nsam, 2);
    cerr << "\n    Estimated memory usage: "
         << estimated_mem * 1.0 / 1024 / 1024 / 1024 << " Gb \n";

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
      this->sample.add_sample(hdr->id[BCF_DT_SAMPLE][i].key, i);

    while ((res = bcf_read(fp, hdr, rec)) == 0) {
      assert(res == 0 && "Can't read record");
      if (verbose) {
        cerr << "\r    Read variants and genotypes: " << rec_count;
        if (rec->errcode)
          cerr << " Errorcode" << rec->errcode;
      }

      res = bcf_get_genotypes(hdr, rec, &p_allele, &count);
      assert(res != 0 && "Can't get genotypes");
      // set alleles
      for (int32_t i = 0; i < count; i++) {
        this->genotype.set_allele(rec_count, i / 2, i % 2,
                                  p_allele[i] >> 1 & 3);
      }
      // set positions
      this->variant.add_pos(
          rec->pos + 1,
          rec_count); // pos is 0-based, while vcf is 1-based. wield
      rec_count++;
      bcf_clear(rec);
    };
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    assert(bcf_close(fp) == 0);
    fp = NULL;
    cerr << "\n    Done parsing VCF file\n";
  }
};

class IbdRec {

public:
  string line, field, sn1, sn2;
  istringstream iss;
  uint32_t sample_id1, sample_id2;
  uint32_t s_pos, e_pos;
  string chr_name;
  Vcf &vcf_ref;
  GeneticMap &gmap_ref;
  gzFile out_file;

  IbdRec(Vcf &vcf, GeneticMap &gmap, const char *out_fn)
      : vcf_ref(vcf), gmap_ref(gmap) {
    line.reserve(10000);
    field.reserve(10000);
    if (strcmp("stdout", out_fn) == 0)
      out_file = Z_NULL;
    else {
      out_file = gzopen(out_fn, "w");
      if (out_file == Z_NULL) {
        cerr << "Cannot open output file to write\n";
        exit(1);
      }
    }
  }

  ~IbdRec() {
    int ret;
    if (out_file != Z_NULL) {
      ret = gzclose(out_file);
      if (ret != Z_OK) {
        cerr << "Close file failured with code: " << ret << '\n';
        exit(1);
      }
    }
  }

  int read_rec(gzFile &fp) {
    char buf[1000];
    if (Z_NULL == gzgets(fp, buf, 1000))
      return -1;
    else {
      line = buf;
      iss.clear();
      iss.str(line);
      getline(iss, sn1, ':');    // sample 1
      getline(iss, sn2, '\t');   // sample 1
      getline(iss, field, '\t'); // start_pos
      s_pos = stoul(field);
      getline(iss, field, '\t'); // end_pos
      e_pos = stoul(field);

      sample_id1 = -1;
      sample_id2 = -1;

      return 0;
    }
  }

  void calc_sample_id() {
    sample_id1 = vcf_ref.sample.get_sample_id(sn1.c_str());
    sample_id2 = vcf_ref.sample.get_sample_id(sn2.c_str());
  }

  void calc_sample_id(IbdRec &prev_rec) {
    if (prev_rec.sn1 == sn1 && prev_rec.sn2 == sn2) {
      sample_id1 = prev_rec.sample_id1;
      sample_id2 = prev_rec.sample_id2;

    } else {
      sample_id1 = vcf_ref.sample.get_sample_id(sn1.c_str());
      sample_id2 = vcf_ref.sample.get_sample_id(sn2.c_str());
    }
  }

  bool is_discordant_site_greater_than(IbdRec &prev_rec, int max_disc_sites) {
    int n_discordant = 0;
    auto [left, right] =
        vcf_ref.variant.get_pos_id_between(prev_rec.e_pos, s_pos);
    for (uint32_t pos_id = left; pos_id < right; pos_id++) {
      int a0, a1, b0, b1;
      a0 = vcf_ref.genotype.get_allele(pos_id, sample_id1, 0);
      a1 = vcf_ref.genotype.get_allele(pos_id, sample_id1, 1);
      b0 = vcf_ref.genotype.get_allele(pos_id, sample_id2, 0);
      b1 = vcf_ref.genotype.get_allele(pos_id, sample_id2, 1);
      if (a0 == a1 && b0 == b1 && a1 != b1 && a1 != 0 && b1 != 0) {
        n_discordant++;
        if (n_discordant > max_disc_sites)
          return true;
      }
    }
    return false;
  }

  float get_cm_dist_from_prev_rec(IbdRec &prev_rec) {
    return gmap_ref.get_cm(s_pos) - gmap_ref.get_cm(prev_rec.e_pos);
  }

  void print_rec() {
    float cM = gmap_ref.get_cm(e_pos) - gmap_ref.get_cm(s_pos);
    if (out_file == Z_NULL)
      cout << sn1 << '\t' << sn2 << '\t' << s_pos << '\t' << e_pos << '\t' << cM
           << '\n';
    else
      gzprintf(out_file, "%s\t%s\t%lu\t%lu\t%0.5g\n", sn1.c_str(), sn2.c_str(),
               s_pos, e_pos, cM);
  }

  bool check_and_merge(IbdRec &prev_rec, float max_cM, int max_disc_sites) {
    bool is_merged = false;
    if (sn1 == prev_rec.sn1 &&
        sn2 == prev_rec.sn2 &&      // if not the same sample pair, won't merge
        (prev_rec.e_pos >= s_pos || // overlapping, merge
         (get_cm_dist_from_prev_rec(prev_rec) <
              max_cM &&                     // no overlapping but close and
          !is_discordant_site_greater_than( // has less than 1 discordance sites
              prev_rec, max_disc_sites)))) {

      if (prev_rec.e_pos < e_pos)
        prev_rec.e_pos = e_pos; // if merge, keep the larger end

      is_merged = true;
    } else {
      prev_rec.print_rec();
      is_merged = false;
    }
    return is_merged;
  }
};

int main(int argc, char *argv[]) {
  int ret;
  uint64_t merge_count = 0, print_count = 0;

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

  Vcf vcf;
  vcf.verbose = arguments.verbose;
  vcf.parse_vcf(arguments.vcf_fn);
  vcf.variant.calc_genetic_pos(gmap);

  IbdRec rec1(vcf, gmap, arguments.out_fn), rec2(vcf, gmap, arguments.out_fn);
  IbdRec *prev = &rec1, *curr = &rec2, *tmp = NULL;

  gzFile fp = gzopen(arguments.ibd_fn, "r");

  cerr << "Merging IBD segments... \n";
  prev->read_rec(fp);
  prev->calc_sample_id();
  while (curr->read_rec(fp) == 0) {
    curr->calc_sample_id(*prev);
    ret = curr->check_and_merge(*prev, arguments.max_cM, arguments.discord);
    if (ret == true) {
      merge_count++;
    } else {
      tmp = prev;
      prev = curr;
      curr = tmp;
      print_count++;
    }
    if ((merge_count + print_count) % 1000000 == 0)
      cerr << "\r    Merged segment: " << merge_count
           << "; printed segment: " << print_count;
  }
  prev->print_rec();
  cerr << "\n    Done merging. Merged: " << merge_count
       << ", printed: " << print_count << '\n';

  gzclose(fp);
  return 0;
}
