#ifndef __gmap_hpp__
#define __gmap_hpp__

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include "common.hpp"

class GeneticMap
{
    int chrom_id;
    std::vector<uint32_t> bp_pos_vec;
    std::vector<float> cm_pos_vec;
    std::vector<float> slope_vec; // cm per bp

  public:
    GeneticMap() {}
    GeneticMap(int chrom_id_, const char *gmap_fn) : chrom_id(chrom_id_)
    {
        // add a initial position
        bp_pos_vec.push_back(0);
        cm_pos_vec.push_back(0.0);
        // slope_vec has 1 less element than pos vectors until the end

        std::string line, field;
        std::ifstream ifs(gmap_fn);
        uint32_t bp_pos;
        float cm_pos;
        size_t line_counter = 0;
        while (std::getline(ifs, line, '\n')) {
            std::istringstream iss(line);
            std::getline(iss, field, ' ');
            std::getline(iss, field, ' ');
            std::getline(iss, field, ' ');
            cm_pos = std::stod(field);
            std::getline(iss, field, ' ');
            bp_pos = std::stoul(field);
            add_position(bp_pos, cm_pos);
        }
        add_final_slope();
    }

    uint32_t
    get_bp(float cm)
    {
        auto id = std::distance(cm_pos_vec.begin(),
            std::upper_bound(cm_pos_vec.begin(), cm_pos_vec.end(), cm));
        id -= 1; // move to map point before the position
        auto slope = slope_vec[id];
        return (cm - cm_pos_vec[id]) / slope + bp_pos_vec[id];
    }

    float
    get_cm(uint32_t bp)
    {
        auto id = std::distance(bp_pos_vec.begin(),
            std::upper_bound(bp_pos_vec.begin(), bp_pos_vec.end(), bp));
        id -= 1; // move to map point before the position
        auto slope = slope_vec[id];
        return (bp - bp_pos_vec[id]) * slope + cm_pos_vec[id];
    }

    void
    write_to_file(BGZF *fp)
    {
        write_element_to_file(chrom_id, fp);
        write_vector_to_file(bp_pos_vec, fp);
        write_vector_to_file(cm_pos_vec, fp);
        write_vector_to_file(slope_vec, fp);
    }

    void
    read_from_file(BGZF *fp)
    {
        read_element_from_file(chrom_id, fp);
        read_vector_from_file(bp_pos_vec, fp);
        read_vector_from_file(cm_pos_vec, fp);
        read_vector_from_file(slope_vec, fp);
    }

    bool
    is_equal(GeneticMap &other)
    {
        return chrom_id == other.chrom_id && bp_pos_vec == other.bp_pos_vec
               && cm_pos_vec == other.cm_pos_vec && slope_vec == other.slope_vec;
    }

    void
    print()
    {
        for (size_t i = 0; i < bp_pos_vec.size(); i++) {
            std::cout << bp_pos_vec[i] << '\t' << cm_pos_vec[i] << '\t' << slope_vec[i]
                      << '\n';
        }
    }

  private:
    // add a poisition and update slope;
    void
    add_position(uint32_t bp_pos, float cm_pos)
    {
        bp_pos_vec.push_back(bp_pos);
        cm_pos_vec.push_back(cm_pos);
        size_t sz = bp_pos_vec.size();
        slope_vec.push_back((cm_pos_vec[sz - 1] - cm_pos_vec[sz - 2])
                            / (bp_pos_vec[sz - 1] - bp_pos_vec[sz - 2]));
    }

    void
    add_final_slope()
    {
        assert(slope_vec.size() + 1 == bp_pos_vec.size()
               && "The GeneticMap might already have added the final slope");
        // extrapolate the last slope
        slope_vec.push_back(slope_vec.back());
    }
};

#endif
