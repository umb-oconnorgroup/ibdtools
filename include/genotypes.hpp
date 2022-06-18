#ifndef __genotypes_hpp__
#define __genotypes_hpp__

#include <htslib/hts.h>
#include <vector>

// genotype of a given sample at a given site is represented by 4 bits
// First 2 bits is for the first haploid; second 2 bits for the 2nd haploid
// For each 2 bits, 0: Ref, 1:Alt, 2: missing or unphased
//
class Genotypes
{
    std::vector<uint8_t> genotype_vec;
    // for quick access
    size_t ncol; // number of samples per position
    // next genotype use low 4 bits
    size_t use_low_bits;

  public:
    Genotypes(size_t ncol_) : ncol(ncol_) { use_low_bits = 1; }
    Genotypes() {}

    void
    add(uint8_t two_alleles)
    { // or 4 bits
        if (use_low_bits)
            // add a char and use low 4 bits for even number
            genotype_vec.push_back(two_alleles);
        else {
            // use high 4 bits for odd number
            genotype_vec.back() |= (two_alleles << 4);
        }
        use_low_bits = use_low_bits == 1 ? 0 : 1;
    }

    uint8_t
    at(uint32_t pid, uint32_t sid)
    {
        // find genotype id (unit is 4 bits)
        size_t genotype_id = pid;
        genotype_id *= ncol;
        genotype_id += sid;

        size_t byte_id = genotype_id >> 1;
        // determine which 4 bits to use
        bool is_high_bits = (genotype_id & 0x1);

        // get the bype containing the 4 bits
        uint8_t res = genotype_vec[byte_id];

        return is_high_bits ? (res >> 4) : (res & 0xf);
    }

    std::vector<std::string> get_haplotypes();

    void print();

    void write_to_file(BGZF *fp);

    void read_from_file(BGZF *fp);

    bool
    is_equal(Genotypes &other)
    {
        return genotype_vec == other.genotype_vec && use_low_bits == other.use_low_bits;
    }
};

#endif
