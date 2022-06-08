#ifndef __ibdstat_hpp__
#define __ibdstat_hpp__
#include "common.hpp"
#include "ibdfile.hpp"
#include "metafile.hpp"
#include "positions.hpp"
#include <algorithm>

class IbdStat
{
    IbdFile ibd;
    MetaFile meta;
    size_t genome_sz_cM;

  public:
    IbdStat(const char *ibd_fn, const char *meta_fn, float genome_sz_cM = 3545.83,
        size_t max_rec = 1 * 1024 * 1024)
    {

        BGZF *fp = bgzf_open(meta_fn, "r");
        verify(fp != NULL && "can't open meta file");
        meta.read_from_file(fp, false);
        bgzf_close(fp);
        fp = NULL;

        ibd = IbdFile(ibd_fn, &meta, max_rec);

        this->genome_sz_cM = ceil(genome_sz_cM);
    }

    std::vector<size_t>
    get_stat()
    {
        auto &pos = ibd.get_meta()->get_positions();
        auto &vec = ibd.get_vec();

        std::vector<size_t> counter; // For return
        counter.resize(genome_sz_cM * 1.5);

        ibd.open("r");

        bool read_full;
        do {
            read_full = ibd.read_from_file();
            for (auto rec : vec) {
                float cM = pos.get_cm(rec.get_pid2()) - pos.get_cm(rec.get_pid1());
                counter[cM]++;
            }

        } while (read_full);

        size_t i;
        for (i = counter.size() - 1; i >= 0; i--) {
            if (counter[i] != 0)
                break;
        }

        counter.resize(i + 1);

        ibd.close();

        return counter;
    }
};

#endif
