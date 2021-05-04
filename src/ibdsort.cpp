#include <algorithm>
#include <argp.h>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <execution>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <memory_resource>
#include <numeric>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <zlib.h>

using namespace std;

// -----------------------------------------
// For argument parsing

const char *argp_program_version = "ibdsort 0.1";
const char *argp_program_bug_address = "gbinux@gmail.com";
static char doc[] =
    "ibdsort prepares formatted and sorted IBD data for ibdmerge. It uses "
    "external merge sort algorithm and is written to handle large datasets.\n";
static char args_doc[] = "IBD_FILE";
static struct argp_option options[] = {
    {0, 0, 0, 0, "Required arguments:"},
    {"IBD_FILE", 1, 0, OPTION_DOC | OPTION_NO_USAGE,
     "Input IBD file from hap_ibd for a single chromosome"},

    {0, 0, 0, 0, "Optional arguments:"},
    {"sample", 's', "filename", 0,
     "Sample file. If provided, IBD records will be sorted according to the "
     "alphabetical order of sample name; if not provodied, IBD samples will be "
     "sorted according to the order the program first encounters the sample "
     "name in the IBD records "
     "Default: NULL"},
    {"tmp_file_prefix", 'p', "string", 0,
     "The tmp file prefix. The file name will be {tmp_file_prefix}_{xx}. "
     "Default is using the IBD_FILE."},
    {"k_ways", 'k', "integer", 0, "Number of ways for merging. Default 8"},
    {"memory", 'm', "float", 0, "Maximum memory in Gb to use. Default 1"},
    {"out", 'o', "filename", 0,
     "IBD output file after formatting and sorting, default stdout"},
    {"keep_tmp_file", 'K', 0, 0, "Keep all tmp file for debuging"},
    {"verbose", 'v', 0, 0, "Show processing updates, default false"},
    {0}};
struct arguments {
  const char *ibd_in_fn;
  const char *sample_fn;
  const char *ibd_out_fn;
  const char *tmp_file_prefix;
  float mem_in_gb;
  int kways;
  bool to_keep_tmp_files;
  bool verbose;
};
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct arguments *arguments = (struct arguments *)state->input;
  switch (key) {
  case 'k':
    arguments->kways = atoi(arg);
    break;
  case 'm':
    arguments->mem_in_gb = atof(arg);
    break;
  case 'o':
    arguments->ibd_out_fn = arg;
    break;
  case 'p':
    arguments->tmp_file_prefix = arg;
    break;
  case 's':
    arguments->sample_fn = arg;
    break;
  case 'v':
    arguments->verbose = true;
    break;
  case 'K':
    arguments->to_keep_tmp_files = true;
    break;
  case ARGP_KEY_ARG:
    switch (state->arg_num) {
    case 0:
      arguments->ibd_in_fn = arg;
      break;
    default:
      argp_usage(state);
    }
    break;
  case ARGP_KEY_END:
    if (state->arg_num != 1)
      argp_usage(state);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

// -----------------------------------------

struct ibdrec_t {
  uint32_t id1;
  uint32_t id2;
  uint32_t start;
  uint32_t end;

  bool is_zero() { return id1 == 0 && id2 == 0 && start == 0 && end == 0; }
  friend bool operator<(const ibdrec_t &rec1, const ibdrec_t &rec2) {
    if (rec1.id1 != rec2.id1)
      return rec1.id1 < rec2.id1;
    if (rec1.id2 != rec2.id2)
      return rec1.id2 < rec2.id2;
    if (rec1.start != rec2.start)
      return rec1.start < rec2.start;
    return rec1.end < rec2.end;
  }

  friend ostream &operator<<(ostream &out, const ibdrec_t &rec) {
    out << rec.id1 << '\t' << rec.id2 << '\t' << rec.start << '\t' << rec.end
        << ' ';
    return out;
  }
};

//////////////////////////////////////////////////////////////////////////////////////////
// K way merging block info, including the block id and the iterator to a rec in
// the block
struct kblock_info_t {
  vector<ibdrec_t>::iterator it;
  int id;
  friend bool operator<(const kblock_info_t &pair1,
                        const kblock_info_t &pair2) {
    return *(pair1.it) < *(pair2.it);
  }
};

//////////////////////////////////////////////////////////////////////////////////////////
// Tmpfile containing a pointer to a file, corrsponding file id
struct TmpFile {
  FILE *fp;
  uint32_t fid;
  uint64_t num_rec;
  int pass;
  string fname;
  TmpFile() {
    fp = NULL;
    num_rec = 0;
    fid = ~0;
  };
  TmpFile(string &prefix, int tmp_file_id) {
    create_tmp_file(prefix, tmp_file_id);
  }
  void create_tmp_file(string &prefix, int tmp_file_id) {
    fname += prefix;
    fname += "_sort_tmp_";
    fname += to_string(tmp_file_id);
    fp = fopen(fname.c_str(), "wb+");
    fid = tmp_file_id;
    num_rec = 0;
  }
};

// custom memory resource
class ListMemResource : public std::pmr::memory_resource {

  std::pmr::memory_resource *res_; // Default resource };
  size_t size;
  struct block_t {
    void *N;
    void *P;
  };
  char *buffer;

public:
  ListMemResource(void *buffer, size_t size)
      : res_{std::pmr::get_default_resource()} {
    this->size = size;
    this->buffer = static_cast<char *>(buffer);
    memset(buffer, 0, size);
  }

private:
  void *do_allocate(std::size_t bytes, std::size_t alignment) override {
    block_t *pblock;
    size_t i;
    size_t block_size =
        (bytes / alignment + (bytes % alignment != 0)) * alignment;

    for (i = 0; i < size; i += block_size) {
      pblock = reinterpret_cast<block_t *>(buffer + i);
      if (pblock->N == NULL)
        break;
    }
    if (i >= size)
      return res_->allocate(bytes, alignment);
    else {
      return buffer + i;
    }
  }
  void do_deallocate(void *p, std::size_t bytes,
                     std::size_t alignment) override {
    if (p < buffer && p > buffer + size)
      return res_->deallocate(p, bytes, alignment);
    else {
      auto ptr = static_cast<block_t *>(p);
      ptr->N = NULL;
      ptr->P = NULL;
    }
  }
  bool
  do_is_equal(const std::pmr::memory_resource &other) const noexcept override {
    return (this == &other);
  }
};

// Use external merge sorting algorithm
// Ref: https://en.wikipedia.org/wiki/External_sorting

class IbdSorter {

  // memory
  size_t max_mem_in_num_rec;

  // For sorting chunks
  vector<ibdrec_t> vec_rec;

  // For K-way merging
  vector<vector<ibdrec_t>> vec_block;
  vector<ibdrec_t> vec_out;
  int K; // K-way merging

  // input and output files
  gzFile &gzf_ibd_in;
  gzFile &gzf_ibd_out;
  gzFile &gzf_samples;

  // Tmp files
  vector<TmpFile> vec_tmp_file;
  string prefix;
  size_t tmp_file_id_counter;
  bool to_keep_tmp_files;

  // sample name /sample id conversion
  unordered_map<string, uint32_t> sample_name2id_map;
  vector<string> vec_sample_names;

  // counter
  size_t num_out_rec;

  // verbose
  bool verbose;

public:
  IbdSorter(size_t max_mem_in_num_rec, int K, gzFile &in, gzFile &out,
            gzFile &sample_file, const char *prefix, bool to_keep_tmp_files,
            bool verbose)
      : max_mem_in_num_rec(max_mem_in_num_rec), K(K), gzf_ibd_in(in),
        gzf_ibd_out(out), gzf_samples(sample_file), prefix(prefix),
        to_keep_tmp_files(to_keep_tmp_files), verbose(verbose) {

    // reserve space for vectors
    // Only use half for the vector, the other half is needed for sorting
    // algorithm.
    vec_rec.reserve(max_mem_in_num_rec / 2);

    // intialize variables
    assert(prefix != NULL);
    tmp_file_id_counter = 0;
    num_out_rec = 0;
  }

  // in sort chunk stage, use half memory as a single block and other half for
  // sorting; in the merge stage, split all memory into K + 1 blocks.
  void mem_rearrange() {
    // release the memory of vec_rec
    vec_rec.resize(0);
    vec_rec.shrink_to_fit();

    // allocate memory of K blocks
    vector<size_t> vec_bin;
    size_t out_buffer_size = max_mem_in_num_rec * 2 / (K + 2);
    vec_out.reserve(out_buffer_size); // out rec buffer
    max_mem_in_num_rec -= out_buffer_size;

    size_t quotient, remainder;
    quotient = max_mem_in_num_rec / K;
    remainder = max_mem_in_num_rec % K;
    for (int i = 0; i < K; i++) {
      vec_block.push_back({}); // in rec buffer
      vec_block[i].reserve(quotient + (i < remainder ? 1 : 0));
      cerr << "chuck[" << i << "] size: " << quotient + (i < remainder ? 1 : 0)
           << '\n';
    }
  }

  // get sample id, add to map if not already existed.
  uint32_t get_sample_id(const string &sample_name) {
    auto it = sample_name2id_map.find(sample_name);
    uint32_t id;
    if (it == sample_name2id_map.end()) {
      if (verbose)
        cerr << sample_name << ": not Found\n";
      id = sample_name2id_map.size();
      sample_name2id_map[sample_name] = id;
      vec_sample_names.push_back(sample_name);
    } else {
      id = it->second;
      if (verbose)
        cerr << sample_name << " found! id = " << it->second
             << " key: " << it->first << " \n";
    }
    return id;
  };

  // add sample name to a vector, sort the vector and fill the sample name to id
  // map.
  void load_sample_names() {
    istringstream iss;
    string field;
    char buffer[1000];
    uint32_t id;

    if (gzf_samples != Z_NULL) {
      while (gzgets(gzf_samples, buffer, 1000) != NULL) {
        iss.clear();
        iss.str(buffer);
        getline(iss, field, '\n');
        if (verbose)
          cerr << "Add sample: " << field << '\n';
        vec_sample_names.push_back(field);
      }

      sort(vec_sample_names.begin(), vec_sample_names.end());

      for (auto &sn : vec_sample_names) {
        id = sample_name2id_map.size();
        sample_name2id_map[sn] = id;
      }
    }
    cerr << "Loaded " << vec_sample_names.size() << " sample name\n";
  }

  // load text data until the vector reach the pre-specified maximum capacity
  void load_text_data() {
    istringstream iss;
    string field;
    int field_no;
    char buffer[1000];
    ibdrec_t rec;
    char *ret;
    size_t num_rec_read = 0;

    while (vec_rec.size() < vec_rec.capacity() &&
           (ret = gzgets(gzf_ibd_in, buffer, 1000)) != NULL) {
      iss.clear();
      iss.str(buffer);
      field_no = 0;

      // read a record and fill the vec_dist vector
      for (int field_no = 0; getline(iss, field, '\t'); field_no++) {
        if (field_no == 0)
          rec.id1 = get_sample_id(field);
        else if (field_no == 2)
          rec.id2 = get_sample_id(field);
        else if (field_no == 5)
          rec.start = stoul(field);
        else if (field_no == 6) {
          rec.end = stoul(field);
          break;
        } else
          continue;
      }
      // to make sure id1 and id2 always occur in the same order
      if (rec.id1 > rec.id2)
        swap(rec.id1, rec.id2);

      // add to vector
      vec_rec.push_back(rec);
      num_rec_read++;
      if (num_rec_read % 1000000 == 0)
        cerr << '\r' << num_rec_read;
    }
    cerr << '\n';
    if (!ret)
      assert(gzeof(gzf_ibd_in));
  }
  // sort all data in the vector
  void sort_data() {
    sort(execution::par_unseq, vec_rec.begin(), vec_rec.end());
  }

  // dump the sorted data from vector and erase the element that is written to a
  // file
  void dump_vector(vector<ibdrec_t> &vec,
                   const vector<ibdrec_t>::iterator &first,
                   const vector<ibdrec_t>::iterator &last, TmpFile &tmp_file) {
    ibdrec_t *p = &(*first);
    fwrite(p, sizeof(ibdrec_t), distance(first, last), tmp_file.fp);
    tmp_file.num_rec += distance(first, last);
    cerr << "DUMPED VECTOR to file # " << tmp_file.fid
         << "\t# of record: " << distance(first, last)
         << "\tfirst record: " << *first << "\tlast record: " << *(last - 1)
         << '\n';
    vec.resize(0);
  }

  // call dump vector, and rewind the file
  void dump_sorted_chunk() {
    auto tmp_file = TmpFile(prefix, tmp_file_id_counter);
    tmp_file_id_counter++;
    dump_vector(vec_rec, vec_rec.begin(), vec_rec.end(), tmp_file);
    fflush(tmp_file.fp);
    rewind(tmp_file.fp);
    vec_tmp_file.push_back(tmp_file);
  }

  void sort_input_into_chunks() {
    // loading sample name should be outside the while loop;
    load_sample_names();
    while (!gzeof(gzf_ibd_in)) {
      load_text_data();
      sort_data();
      dump_sorted_chunk();
    }
  }

  void load_sorted_chunk(int id, TmpFile &tmp_file) {
    size_t num_rec_read;
    size_t num_remaining;
    vector<ibdrec_t> &vec = vec_block[id];
    num_remaining = vec.size();
    vec.resize(vec.capacity());
    num_rec_read = fread(&vec[num_remaining], sizeof(ibdrec_t),
                         vec.capacity() - num_remaining, tmp_file.fp);
    vec.resize(num_remaining + num_rec_read);
    if (num_rec_read > 0)
      cerr << "Load from file # " << tmp_file.fid
           << "\tnum_rec_read: " << num_rec_read
           << "\tfirst record: " << vec[num_remaining] << '\n';
  }

  void send_to_out_buffer(ibdrec_t &rec, TmpFile &tmp_file) {
    bool is_buffer_full = vec_out.size() >= vec_out.capacity();
    if (!is_buffer_full)
      vec_out.push_back(rec);
    else {
      dump_vector(vec_out, vec_out.begin(), vec_out.end(), tmp_file);
      // Be sure to push_back after dump_vector
      vec_out.push_back(rec);
    }
  }

  void flush_out_buffer(TmpFile &tmp_file) {
    dump_vector(vec_out, vec_out.begin(), vec_out.end(), tmp_file);
    // important to flush
    fflush(tmp_file.fp);
    rewind(tmp_file.fp);
  }

  TmpFile merge_K_chunks(vector<TmpFile>::iterator first,
                         vector<TmpFile>::iterator last, bool is_final_pass) {
    bool use_stdout = (gzf_ibd_out == Z_NULL);

    TmpFile tmp_file;
    if (!is_final_pass)
      tmp_file.create_tmp_file(prefix, tmp_file_id_counter++);

    // copy Tmp file to working vector
    size_t num_files = distance(first, last);
    vector<TmpFile> vec_working_files;
    vec_working_files.reserve(num_files);
    copy(first, last, back_inserter(vec_working_files));

    // NOTE: this is critical to improve the efficiency of the list. The idea is
    // preallocate memory block and make memory resource based on the block.
    // Then utilize the the pmr list to use this allocated memory resource. When
    // deleting a node, we just mark the sub block to a certain value; when
    // "allocate" memory for a new, we just find a sub block with mark value.
    // This way we avoid frequeny allocate and deallocate memory from the raw
    // memory resouce.
    vector<tuple<void *, void *, kblock_info_t>> vec_mem(num_files * 5 / 4);
    ListMemResource lst_rsc(&vec_mem[0],
                            sizeof(tuple<void *, void *, kblock_info_t>) *
                                vec_mem.size());
    pmr::list<kblock_info_t> lst_block_info(&lst_rsc);

    // load data, initialize block info list
    for (int i = 0; i < num_files; i++) {
      load_sorted_chunk(i, vec_working_files[i]);
      if (vec_block[i].size() == 0)
        continue;
      kblock_info_t kblock_info = {vec_block[i].begin(), i};
      // insert and keep in order
      auto pos = upper_bound(lst_block_info.begin(), lst_block_info.end(),
                             kblock_info);
      lst_block_info.insert(pos, kblock_info);
    }

    // while loop to find the minial record in K blocks and write it output file
    while (lst_block_info.size() > 0) {
      if (verbose)
        cerr << "lst_block_info size(): " << lst_block_info.size()
             << ", vec_tmp_file size: " << vec_tmp_file.size() << '\n';
      auto kblock_info =
          *lst_block_info.begin(); // for the block with minist record
      // assert(lst_block_info.begin() == min_element(lst_block_info.begin(),
      // lst_block_info.end()));
      ibdrec_t min_rec = *kblock_info.it;
      lst_block_info.erase(lst_block_info.begin()); // remove it from the list

      bool to_update_lst = true;
      if (!is_final_pass)
        send_to_out_buffer(*(kblock_info.it), tmp_file);
      else {
        if (use_stdout)
          cout << vec_sample_names[min_rec.id1] << '\t'
               << vec_sample_names[min_rec.id2] << '\t' << min_rec.start << '\t'
               << min_rec.end << '\n';
        // cout << "COUT: " << rec << '\n';
        else {
          gzprintf(gzf_ibd_out, "%s\t%s\t%lu\t%lu\n",
                   vec_sample_names[kblock_info.it->id1].c_str(),
                   vec_sample_names[kblock_info.it->id2].c_str(),
                   kblock_info.it->start, kblock_info.it->end);
        }
        num_out_rec++;
      }

      // move iterator to next record
      kblock_info.it++;

      // 1. if reaching to buffer end: a) also reaching to file end, the close
      // file;  b) not reaching file end, load new records, and iterator point
      // to the begin of the buffer. 2. If not reaching buffer end,  iterator
      // points to next record in the buffer.
      if (kblock_info.it >= vec_block[kblock_info.id].end()) {
        int id = kblock_info.id;
        FILE *fp = vec_working_files[id].fp;
        vec_block[id].resize(0);
        load_sorted_chunk(id, vec_working_files[id]);

        if (vec_block[id].size() > 0)
          kblock_info.it = vec_block[id].begin();
        else {
          to_update_lst = false;
          fclose(fp);
          vec_working_files[id].fp = NULL;
          cerr << "To KEEP TMP FILES: " << to_keep_tmp_files << '\n';
          if (!to_keep_tmp_files)
            remove(vec_working_files[id].fname.c_str());
        }
      }
      if (to_update_lst) {
        auto pos = upper_bound(lst_block_info.begin(), lst_block_info.end(),
                               kblock_info);
        // assert(lst_block_info.begin() == min_element(lst_block_info.begin(),
        // lst_block_info.end()));
        lst_block_info.insert(pos, kblock_info);
      }
    }

    // flush out buffer to file when the above loop is done
    if (!is_final_pass && !use_stdout)
      flush_out_buffer(tmp_file);

    return tmp_file;
  }

  void merge_all_chunks() {
    size_t vec_size;
    vector<TmpFile> vec_tmp_file_new;

    while ((vec_size = vec_tmp_file.size()) > 3) {
      vector<int> vec_index;
      for (int i = 0; i < vec_size; i += K)
        vec_index.push_back(i);

      transform(vec_index.begin(), vec_index.end(),
                back_inserter(vec_tmp_file_new), [&](int id) -> TmpFile {
                  auto first = vec_tmp_file.begin() + id;
                  auto last = vec_tmp_file.begin() + id + K;
                  if (last > vec_tmp_file.end())
                    last = vec_tmp_file.end();
                  return merge_K_chunks(first, last, false);
                });

      swap(vec_tmp_file, vec_tmp_file_new);
      vec_tmp_file_new.resize(0);
    }

    // last pass.
    // If not use vec_out buffer, release it and set the gzbuffer
    if (gzf_ibd_out == Z_NULL) {
      vec_out.resize(0);
      vec_out.shrink_to_fit();
    }
    merge_K_chunks(vec_tmp_file.begin(), vec_tmp_file.end(), true);

    cerr << "Total # output record: " << num_out_rec << '\n';
  }
};

int main(int argc, char *argv[]) {

  gzFile fp_ibd_in, fp_ibd_out, fp_sample;
  struct arguments arguments;
  size_t max_lines;

  arguments.kways = 8;
  arguments.mem_in_gb = 1;
  arguments.sample_fn = NULL;
  arguments.ibd_out_fn = NULL;
  arguments.tmp_file_prefix = NULL;
  arguments.to_keep_tmp_files = false;
  arguments.verbose = false;
  argp_parse(&argp, argc, argv, 0, 0, &arguments);
  if (arguments.tmp_file_prefix == NULL)
    arguments.tmp_file_prefix = arguments.ibd_in_fn;

  fp_ibd_in = gzopen(arguments.ibd_in_fn, "r");
  assert(fp_ibd_in != Z_NULL);

  max_lines = arguments.mem_in_gb * 1024 * 1024 * 1024 / sizeof(ibdrec_t);

  if (arguments.ibd_out_fn) {
    fp_ibd_out = gzopen(arguments.ibd_out_fn, "w");
    assert(fp_ibd_out != Z_NULL);
  } else
    fp_ibd_out = Z_NULL;

  if (arguments.sample_fn) {
    fp_sample = gzopen(arguments.sample_fn, "r");
    cerr << "Sample file name" << arguments.sample_fn << '\n';
    assert(fp_sample != Z_NULL);
  } else
    fp_sample = Z_NULL;

  IbdSorter ibdsorter(max_lines, arguments.kways, fp_ibd_in, fp_ibd_out,
                      fp_sample, arguments.tmp_file_prefix,
                      arguments.to_keep_tmp_files, arguments.verbose);
  ibdsorter.sort_input_into_chunks();
  ibdsorter.mem_rearrange();
  ibdsorter.merge_all_chunks();

  gzclose(fp_ibd_in);
  gzclose(fp_ibd_out);
  gzclose(fp_sample);

  cerr << sizeof(kblock_info_t) << '\n';

  return 0;
}
