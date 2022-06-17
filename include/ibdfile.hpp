#ifndef __ibdfile_hpp__
#define __ibdfile_hpp__

#include "common.hpp"

class MetaFile;

// IbdFile object represent an packed IBD file
class IbdFile
{
    std::string filename;
    BGZF *fp;

    // buffer
    std::vector<ibd_rec1_t> ibd_vec;

    // pointer
    MetaFile *meta;

  public:
    IbdFile(){};
    IbdFile(
        const char *ibd_fn, MetaFile *meta_ = NULL, size_t max_rec = 1 * 1024 * 1024);

    std::vector<ibd_rec1_t> &get_vec();

    const std::string &get_filename();

    BGZF *get_fp();

    MetaFile *get_meta();

    void open(const char *mode = "r");

    void close();

    void from_raw_ibd(const char *raw_ibd_in, int col_sample1 = 0, int col_sample2 = 2,
        int col_start = 5, int col_end = 6, int col_hap1 = 1, int col_hap2 = 3);
    void to_raw_ibd(const char *raw_ibd_fn,
        size_t line_buffer_capcity = 10 * 1024 * 1024, const char *subpop_fn = NULL);

    void write_to_file(size_t max = 0);

    void write_to_file(IbdFile &other, size_t max = 0);

    // if new_capcity >= current capcity then enlarge the capacity
    // if append is false
    // @return  true if vector is full; false if vector not full, which
    // indicates end of file
    bool read_from_file(bool append = false, size_t new_capacity = 0);

    void delete_file_from_disk();

    bool has_equal_vector(const IbdFile &other);

    void head(size_t n);

    void tail(size_t n);

    void summary();

    void print_to_file(std::ofstream &ofs);
};

#endif
