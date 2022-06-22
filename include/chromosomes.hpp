#ifndef __chromosomes_hpp__
#define __chromosomes_hpp__

#include "common.hpp"
#include <htslib/bgzf.h>
#include <string>
#include <vector>

class Chromosomes
{
    std::vector<std::string> names_vec;
    std::vector<uint32_t> len_bp_vec;
    std::vector<float> len_cM_vec;

  public:
    void add(std::string name, uint32_t bp_length, float cM_length);

    std::string &
    get_name(const ssize_t id)
    {
        return names_vec.at(id);
    }

    size_t
    get_id(const std::string name)
    {
        auto it = std::find(names_vec.begin(), names_vec.end(), name);
        my_assert(it != names_vec.end(), "");
        return std::distance(names_vec.begin(), it);
    }

    void write_to_file(BGZF *fp);

    void read_from_file(BGZF *fp);

    bool
    is_equal(Chromosomes &other)
    {
        return names_vec == other.names_vec && len_bp_vec == other.len_bp_vec
               && len_cM_vec == other.len_cM_vec;
    }
};

#endif
