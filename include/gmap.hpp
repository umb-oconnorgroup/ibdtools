#ifndef __gmap_hpp__
#define __gmap_hpp__

#include <cmath>
#include <htslib/hts.h>
#include <vector>

void exit_on_false(bool condition, const char *message, const char *file, int lineno);

class GeneticMap
{
    int chrom_id;
    std::vector<size_t> bp_pos_vec;
    std::vector<long double> cm_pos_vec; // internal using unit 10^(-8) cM
    std::vector<long double> slope_vec;  // cm per bp

  public:
    GeneticMap() {}
    GeneticMap(int chrom_id_, const char *gmap_fn);
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

    void print_range_info(size_t lower_id);

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

    void write_to_file(BGZF *fp);

    void read_from_file(BGZF *fp);

    bool
    is_equal(GeneticMap &other)
    {
        return chrom_id == other.chrom_id && bp_pos_vec == other.bp_pos_vec
               && cm_pos_vec == other.cm_pos_vec && slope_vec == other.slope_vec;
    }

    void print();
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
        exit_on_false(slope_vec.size() + 1 == bp_pos_vec.size(),
            "The GeneticMap might already have added the final slope", __FILE__,
            __LINE__);
        // extrapolate the last slope
        slope_vec.push_back(slope_vec.back());
    }
};

#endif
