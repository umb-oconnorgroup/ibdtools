#ifndef __metafile_hpp__
#define __metafile_hpp__

#include <htslib/bgzf.h>
#include <htslib/vcf.h>

#include <iostream>

#include "chromosomes.hpp"
#include "common.hpp"
#include "genotypes.hpp"
#include "gmap.hpp"
#include "positions.hpp"
#include "samples.hpp"

class MetaFile
{
    Chromosomes chromosomes;
    GeneticMap gmap;
    Samples samples;
    Positions positions;
    Genotypes genotypes;

  public:
    MetaFile() {}
    void parse_files(const char *vcf_fn, const char *gmap_fn, bool parse_genotype = true,
        std::string chr_name = "0");
    void write_to_file(BGZF *fp);

    bool is_equal(MetaFile &other);

    void read_from_file(BGZF *fp, bool read_genotype = true);

    void print();

    // getters
    Genotypes &
    get_genotypes()
    {
        return genotypes;
    }
    Samples &
    get_samples()
    {
        return samples;
    }
    Positions &
    get_positions()
    {
        return positions;
    }
    GeneticMap &
    get_genetic_map()
    {
        return gmap;
    }
    Chromosomes &
    get_chromosomes()
    {
        return chromosomes;
    }
};

#endif
