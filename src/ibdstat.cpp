
#include "ibdstat.hpp"
#include "chromosomes.hpp"
#include "common.hpp"
#include "genotypes.hpp"
#include "gmap.hpp"
#include "ibdfile.hpp"
#include "metafile.hpp"
#include "positions.hpp"
#include "samples.hpp"
#include <algorithm>

IbdStat::IbdStat(
    const char *ibd_fn, const char *meta_fn, float genome_sz_cM, size_t max_rec)
{

    BGZF *fp = bgzf_open(meta_fn, "r");
    my_assert(fp != NULL, "can't open meta file");
    meta.read_from_file(fp, false);
    bgzf_close(fp);
    fp = NULL;

    ibd = IbdFile(ibd_fn, &meta, max_rec);

    this->genome_sz_cM = ceil(genome_sz_cM);
}

std::vector<size_t>
IbdStat::get_stat()
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
