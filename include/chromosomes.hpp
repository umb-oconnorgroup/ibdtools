#ifndef __chromosomes_hpp__
#define __chromosomes_hpp__

#include "common.hpp"

class Chromosomes
{
    std::vector<std::string> names_vec;
    std::vector<uint32_t> len_bp_vec;
    std::vector<float> len_cM_vec;

  public:
    void
    add(std::string name, uint32_t bp_length, float cM_length)
    {
        auto it = std::find(names_vec.begin(), names_vec.end(), name);

        if (it == names_vec.end()) {
            // add new chromomes
            names_vec.push_back(name);
            len_bp_vec.push_back(bp_length);
            len_cM_vec.push_back(cM_length);
        } else {
            // update existing one
            size_t id = distance(names_vec.begin(), it);
            len_bp_vec[id] = bp_length;
            len_cM_vec[id] = cM_length;
        }
    }

    std::string &
    get_name(const ssize_t id)
    {
        return names_vec.at(id);
    }

    size_t
    get_id(const std::string name)
    {
        auto it = std::find(names_vec.begin(), names_vec.end(), name);
        verify(it != names_vec.end());
        return std::distance(names_vec.begin(), it);
    }

    void
    write_to_file(BGZF *fp)
    {
        write_vector_to_file(names_vec, fp);
        write_vector_to_file(len_bp_vec, fp);
        write_vector_to_file(len_cM_vec, fp);
    }

    void
    read_from_file(BGZF *fp)
    {
        read_vector_from_file(names_vec, fp);
        read_vector_from_file(len_bp_vec, fp);
        read_vector_from_file(len_cM_vec, fp);
    }

    bool
    is_equal(Chromosomes &other)
    {
        return names_vec == other.names_vec && len_bp_vec == other.len_bp_vec
               && len_cM_vec == other.len_cM_vec;
    }
};

#endif
