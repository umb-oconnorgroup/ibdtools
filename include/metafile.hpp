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
    void
    parse_files(const char *vcf_fn, const char *gmap_fn, bool parse_genotype = true,
        std::string chr_name = "0")
    {
        // ScopedTimer timer("vcffile-parse_files");

        htsFile *htsfp = bcf_open(vcf_fn, "r");
        assert(htsfp != NULL && "Cannot open vcf file");

        // assert(hts_set_threads(htsfp, 3) == 0);

        bcf_hdr_t *header = bcf_hdr_read(htsfp);
        assert(header != NULL && "Cannot read header");
        auto nsam = header->n[BCF_DT_SAMPLE];

        genotypes = Genotypes(nsam);
        chromosomes.add(chr_name, 0, 0);
        positions = Positions(chromosomes.get_id(chr_name));
        GeneticMap gmap(-1, gmap_fn);

        // add names to the Samples object
        for (int32_t i = 0; i < nsam; i++)
            samples.add(header->id[BCF_DT_SAMPLE][i].key);

        bcf1_t *rec = bcf_init();
        assert(rec != NULL && "Cannot initialize bcf record");

        // loop over each record
        int32_t *dest = NULL;
        int32_t count = 0;
        while (bcf_read(htsfp, header, rec) == 0) {
            const int32_t max_ploidy = 2;

            // 1. get positions
            // This took a while to debug but bcftools internally use 0-based, while
            // ibd/vcf use 1-based;
            uint32_t bp_pos = rec->pos + 1;

            // some vcf file have duplicated variants for a given position; Only the
            // first variant at this position is considered for genotype information
            bool success = positions.add(bp_pos, gmap.get_cm(bp_pos));

            if (parse_genotype && success) {
                // 2. get an array of alleles
                bcf_get_genotypes(header, rec, &dest, &count);
                assert(max_ploidy * nsam == count && "Max ploidy is not 2");

                // 3. loop over allele array and store the alleles
                for (int32_t i = 0; i < count; i += 2) {
                    uint8_t allele[2];

                    allele[0] = bcf_gt_allele(dest[i]);
                    if (bcf_gt_is_missing(dest[i]))
                        allele[0] = 2;

                    allele[1] = bcf_gt_allele(dest[i + 1]);
                    if (bcf_gt_is_missing(dest[i + 1]))
                        allele[1] = 2;

                    genotypes.add(allele[0] | (allele[1] << 2));
                }
            }
        }
        if (dest != NULL)
            free(dest);

        bcf_hdr_destroy(header);
        bcf_destroy(rec);
        bcf_close(htsfp);
    }

    void
    write_to_file(BGZF *fp)
    {
        chromosomes.write_to_file(fp);
        gmap.write_to_file(fp);
        samples.write_to_file(fp);
        positions.write_to_file(fp);
        genotypes.write_to_file(fp);
    }

    bool
    is_equal(MetaFile &other)
    {
        return get_chromosomes().is_equal(other.get_chromosomes())
               && get_genetic_map().is_equal(other.get_genetic_map())
               && get_genotypes().is_equal(other.get_genotypes())
               && get_positions().is_equal(other.get_positions())
               && get_samples().is_equal(other.get_samples());
    }

    void
    read_from_file(BGZF *fp)
    {
        chromosomes.read_from_file(fp);
        gmap.read_from_file(fp);
        samples.read_from_file(fp);
        positions.read_from_file(fp);
        genotypes.read_from_file(fp);
    }

    void
    print()
    {
        samples.print();
        positions.print();
        genotypes.print();
    }

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
