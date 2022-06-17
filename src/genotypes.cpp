#include "genotypes.hpp"
#include "common.hpp"
// genotype of a given sample at a given site is represented by 4 bits
// First 2 bits is for the first haploid; second 2 bits for the 2nd haploid
// For each 2 bits, 0: Ref, 1:Alt, 2: missing or unphased
//
std::vector<std::string>
Genotypes::get_haplotypes()
{
    std::vector<std::string> hap_vec;
    size_t num_alleles_per_row = ncol * 2;
    size_t num_alleles = genotype_vec.size() * 4 - (use_low_bits ? 0 : 2);
    size_t nrow = num_alleles / num_alleles_per_row;

    hap_vec.reserve(nrow);
    for (size_t row = 0; row < nrow; row++) {
        std::string line;
        line.reserve(num_alleles_per_row);

        for (size_t col = 0; col < ncol; col++) {
            uint8_t gt = at(row, col);
            line += std::to_string(gt & 0b11); // first allele
            line += std::to_string(gt >> 2);   // second allele
        }
        hap_vec.push_back(line);
    }

    return hap_vec;
}

void
Genotypes::print()
{
    size_t gt_counter = 0;
    for (auto byte : genotype_vec) {
        std::cout << (byte & 0b11) << (byte >> 2 & 0b11) << ' ';
        if (++gt_counter == ncol) {
            gt_counter = 0;
            std::cout << '\n';
        }
        std::cout << (byte >> 4 & 0b11) << (byte >> 6) << ' ';
        if (++gt_counter == ncol) {
            gt_counter = 0;
            std::cout << '\n';
        }
    }
}

void
Genotypes::write_to_file(BGZF *fp)
{
    write_vector_to_file(genotype_vec, fp);
    write_element_to_file(use_low_bits, fp);
    write_element_to_file(ncol, fp);
}

void
Genotypes::read_from_file(BGZF *fp)
{
    read_vector_from_file(genotype_vec, fp);
    read_element_from_file(use_low_bits, fp);
    read_element_from_file(ncol, fp); // needed to calculate byte-id from sid and pid
}
