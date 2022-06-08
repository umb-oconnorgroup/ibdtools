#ifndef __gmap_hpp__
#define __gmap_hpp__

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "common.hpp"

class GeneticMap
{
    int chrom_id;
    std::vector<size_t> bp_pos_vec;
    std::vector<long double> cm_pos_vec; // internal using unit 10^(-8) cM
    std::vector<long double> slope_vec;  // cm per bp

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
        size_t bp_pos;
        long double cm_pos;
        size_t line_counter = 0;
        while (std::getline(ifs, line, '\n')) {
            if (line_counter < 10) {
                verify((line.npos == line.find_first_of('\t'))
                       && "Error in parsing plink map. Found tab but should use space "
                          "as column delimiter");
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

    size_t
    get_bp(long double cm)
    {
        auto id = std::distance(cm_pos_vec.begin(),
            std::upper_bound(cm_pos_vec.begin(), cm_pos_vec.end(), cm));
        id -= 1; // move to map point before the position
        auto slope = slope_vec[id];

        // if (llround((cm - cm_pos_vec[id]) / slope + bp_pos_vec[id]) > bp_pos_vec[id +
        // 1]
        //     && id < bp_pos_vec.size() - 2) {
        //     print_range_info(id);
        //     std::cout << "bp: "
        //               << llround((cm - cm_pos_vec[id]) / slope + bp_pos_vec[id])
        //               << " cm: " << cm << "\n---------\n";
        // }
        return llround((cm - cm_pos_vec[id]) / slope + bp_pos_vec[id]);
    }

    // Like this region
    // 10 . 82.484688 58212057
    // 10 . 82.484688 58212123
    // 10 . 82.484688 58212538
    // 10 . 82.484688 58212655
    bool
    isin_horizontal_region(long double cm)
    {
        return std::distance(std::lower_bound(cm_pos_vec.begin(), cm_pos_vec.end(), cm),
                   std::upper_bound(cm_pos_vec.begin(), cm_pos_vec.end(), cm))
               > 1;
    }

    size_t
    get_lower_id(size_t bp)
    {
        auto id = std::distance(bp_pos_vec.begin(),
            std::upper_bound(bp_pos_vec.begin(), bp_pos_vec.end(), bp));
        return id - 1;
    }

    size_t
    get_lower_id(long double cm)
    {
        auto id = std::distance(cm_pos_vec.begin(),
            std::upper_bound(cm_pos_vec.begin(), cm_pos_vec.end(), cm));
        return id - 1;
    }

    void
    print_range_info(size_t lower_id)
    {
        verify(lower_id >= 0);
        if (lower_id >= bp_pos_vec.size() - 1) {
            std::cout << std::setprecision(10) << "bp [" << bp_pos_vec.back() << ", Inf)"
                      << " cm [" << cm_pos_vec.back()
                      << ", Inf) with a slope: " << slope_vec.back() << '\n';
        } else {
            std::cout << "bp [" << bp_pos_vec[lower_id] << ", "
                      << bp_pos_vec[lower_id + 1] << ")"
                      << " cm [" << cm_pos_vec[lower_id] << ", "
                      << cm_pos_vec[lower_id + 1]
                      << ") with a slope: " << slope_vec[lower_id] << '\n';
        }
    }

    long double
    get_cm(size_t bp)
    {
        auto id = std::distance(bp_pos_vec.begin(),
            std::upper_bound(bp_pos_vec.begin(), bp_pos_vec.end(), bp));
        id -= 1; // move to map point before the position
        auto slope = slope_vec[id];

        // debug
        // if (((bp - bp_pos_vec[id]) * slope + cm_pos_vec[id])
        //     > cm_pos_vec[id + 1] + 0.0003) {
        //     print_range_info(id);
        //     std::cout << "bp: " << bp
        //               << " cm: " << ((bp - bp_pos_vec[id]) * slope + cm_pos_vec[id])
        //               << "\n---------\n";
        // }

        return ((bp - bp_pos_vec[id]) * slope + cm_pos_vec[id]);
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

    size_t
    get_first_nonzero_bp()
    {
        return bp_pos_vec[1];
    }
    size_t
    get_last_bp()
    {
        return bp_pos_vec.back();
    }
    long double
    get_first_nonzero_cm()
    {
        return cm_pos_vec[1];
    }
    long double
    get_last_cm()
    {
        return cm_pos_vec.back();
    }

  private:
    // add a poisition and update slope;
    void
    add_position(size_t bp_pos, long double cm_pos)
    {
        bp_pos_vec.push_back(bp_pos);
        cm_pos_vec.push_back(cm_pos);
        size_t sz = bp_pos_vec.size();
        long double slope = (cm_pos_vec[sz - 1] - cm_pos_vec[sz - 2])
                            / (bp_pos_vec[sz - 1] - bp_pos_vec[sz - 2]);
        slope_vec.push_back(slope);
    }

    void
    add_final_slope()
    {
        verify(slope_vec.size() + 1 == bp_pos_vec.size()
               && "The GeneticMap might already have added the final slope");
        // extrapolate the last slope
        slope_vec.push_back(slope_vec.back());
    }
};

#endif
