#ifndef __ibdmerger_hpp__
#define __ibdmerger_hpp__

#include "./common.hpp"
#include "./ibdfile.hpp"
#include "./metafile.hpp"

class IbdMerger
{
    IbdFile in;
    IbdFile out;

    MetaFile meta;

    size_t max_rec_allowed_by_ram;
    std::string out_mode;

    struct range_t {
        size_t first, last;
    };

    std::vector<range_t> range_vec;

    // threshold
    float max_cm;
    size_t max_snp;

    std::vector<ibd_rec1_t> &vec{ in.get_vec() };

  public:
    IbdMerger(const char *in_fn, const char *out_fn, const char *out_mode,
        const char *meta_fn, size_t max_rec_allowed_by_ram, size_t max_snp = 1,
        float max_cm = 0.6);

    void merge();

  private:
    size_t find_group_start(size_t id);

    bool is_mergeable(size_t id1, size_t id2);
    void merge_range(size_t first, size_t last);
};

#endif
