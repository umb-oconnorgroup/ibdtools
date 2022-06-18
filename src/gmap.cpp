#include "gmap.hpp"
#include "common.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

GeneticMap::GeneticMap(int chrom_id_, const char *gmap_fn) : chrom_id(chrom_id_)
{
    // add a initial position
    bp_pos_vec.push_back(0);
    cm_pos_vec.push_back(0.0);
    // slope_vec has 1 less element than pos vectors until the end

    std::string line, field;
    std::ifstream ifs(gmap_fn);
    size_t bp_pos;
    long double cm_pos;
    size_t line_counter = 0;
    while (std::getline(ifs, line, '\n')) {
        if (line_counter < 10) {
            exit_on_false((line.npos == line.find_first_of('\t')),
                "Error in parsing plink map. Found tab but should use space "
                "as column delimiter",
                __FILE__, __LINE__);
        }
        std::istringstream iss(line);
        std::getline(iss, field, ' ');
        std::getline(iss, field, ' ');
        std::getline(iss, field, ' ');
        cm_pos = std::stold(field);
        std::getline(iss, field, ' ');
        bp_pos = std::stoul(field);
        add_position(bp_pos, cm_pos);
        line_counter += 1;
    }
    add_final_slope();
}

void
GeneticMap::print_range_info(size_t lower_id)
{
    exit_on_false(lower_id >= 0, "", __FILE__, __LINE__);
    if (lower_id >= bp_pos_vec.size() - 1) {
        std::cout << std::setprecision(10) << "bp [" << bp_pos_vec.back() << ", Inf)"
                  << " cm [" << cm_pos_vec.back()
                  << ", Inf) with a slope: " << slope_vec.back() << '\n';
    } else {
        std::cout << "bp [" << bp_pos_vec[lower_id] << ", " << bp_pos_vec[lower_id + 1]
                  << ")"
                  << " cm [" << cm_pos_vec[lower_id] << ", " << cm_pos_vec[lower_id + 1]
                  << ") with a slope: " << slope_vec[lower_id] << '\n';
    }
}

void
GeneticMap::write_to_file(BGZF *fp)
{
    write_element_to_file(chrom_id, fp);
    write_vector_to_file(bp_pos_vec, fp);
    write_vector_to_file(cm_pos_vec, fp);
    write_vector_to_file(slope_vec, fp);
}

void
GeneticMap::read_from_file(BGZF *fp)
{
    read_element_from_file(chrom_id, fp);
    read_vector_from_file(bp_pos_vec, fp);
    read_vector_from_file(cm_pos_vec, fp);
    read_vector_from_file(slope_vec, fp);
}

void
GeneticMap::print()
{
    for (size_t i = 0; i < bp_pos_vec.size(); i++) {
        std::cout << bp_pos_vec[i] << '\t' << cm_pos_vec[i] << '\t' << slope_vec[i]
                  << '\n';
    }
}
