#include "chromosomes.hpp"
#include "common.hpp"

void
Chromosomes::add(std::string name, uint32_t bp_length, float cM_length)
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

void
Chromosomes::write_to_file(BGZF *fp)
{
    write_vector_to_file(names_vec, fp);
    write_vector_to_file(len_bp_vec, fp);
    write_vector_to_file(len_cM_vec, fp);
}

void
Chromosomes::read_from_file(BGZF *fp)
{
    read_vector_from_file(names_vec, fp);
    read_vector_from_file(len_bp_vec, fp);
    read_vector_from_file(len_cM_vec, fp);
}
