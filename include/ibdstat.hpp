#ifndef __ibdstat_hpp__
#define __ibdstat_hpp__
#include "ibdfile.hpp"
#include "metafile.hpp"

class IbdStat
{
    IbdFile ibd;
    MetaFile meta;
    size_t genome_sz_cM;

  public:
    IbdStat(const char *ibd_fn, const char *meta_fn, float genome_sz_cM = 3545.83,
        size_t max_rec = 1 * 1024 * 1024);

    std::vector<size_t> get_stat();
};

#endif
