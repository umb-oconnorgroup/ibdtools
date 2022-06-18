#ifndef __ibdsorter_hpp__
#define __ibdsorter_hpp__
#include <common.hpp>
#include <memory>
#include <string>
#include <vector>

class IbdFile;

// Note:
// Be sure to call constructor
// Explitly use auto & instead of auto for reference
// File id should not use vector size as it needs to erase already merged items.
class IbdSorter
{
    std::unique_ptr<IbdFile> in;
    std::unique_ptr<IbdFile> out;
    std::vector<IbdFile> chunks;

    size_t max_rec_allowed_by_ram;
    std::string chunk_prefix;
    std::string out_mode;
    size_t counter;

    ibd_rec_cmp_t Cmp;

  private:
    void merge_first_k_chunk(IbdFile &out_put, size_t kways);

  public:
    IbdSorter(const char *in_fn, const char *out_fn, const char *out_mode,
        const char *chunk_fn_prefix, size_t max_rec_allowed_by_ram,
        ibd_rec_cmp_t Cmp = std::less<ibd_rec1_t>());

    void sort_into_chunks();
    void merge_chunks(size_t kways = 10);
};

#endif
