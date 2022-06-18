#ifndef __ibdsplitter_hpp__
#define __ibdsplitter_hpp__

#include "./common.hpp"
#include "metafile.hpp"
#include "positions.hpp"
#include <algorithm>
#include <memory>
class IbdFile;

class IbdSplitter
{
    std::vector<region_label_t> &labels;
    std::unique_ptr<IbdFile> in;
    std::vector<IbdFile> out_files;
    std::unique_ptr<MetaFile> meta;

    float cm_threshold;

  public:
    IbdSplitter(const char *in_fn, const char *out_fn_prefix, const char *meta_fn,
        std::vector<region_label_t> &labels_, uint32_t max_label = 1,
        float cm_threshold = 2.0, size_t in_ibdrec_vec_capacity = 1 * 1024 * 1024,
        size_t out_ibdrec_vec_capacity = 1 * 1024 * 1024);

    void print_labels();

    // driving code
    void split();

    void split_ibd_rec(ibd_rec1_t rec);

    void transfer(ibd_rec1_t rec, uint32_t rid);

  private:
    bool
    is_longer_than_2cm(uint32_t pid_l, uint32_t pid_r)
    {
        auto &pos = meta->get_positions();
        return pos.get_cm(pid_r) - pos.get_cm(pid_l) >= cm_threshold;
    }
};

#endif
