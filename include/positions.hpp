#ifndef __positions_hpp__
#define __positions_hpp__

#include "common.hpp"
#include <unordered_map>

// Positions for a given chromosome
class Positions
{
    int chrom_id;
    std::vector<uint32_t> pos_bp_vec;
    std::vector<float> pos_cm_vec;
    // for find pos id
    std::unordered_map<uint32_t, uint32_t> bp2id_map;

  public:
    // set chrom_id_ to -1 if don't care
    Positions(int chrom_id_) : chrom_id(chrom_id_) {}
    Positions() {}

    // @return false if pos_bp is duplicated; true otherwise
    bool
    add(uint32_t pos_bp, float pos_cm)
    {
        if (bp2id_map.find(pos_bp) != bp2id_map.end()) {
            std::cerr << "==> Warning: Duplicated position; only the first will be used!"
                      << " duplicated bp: " << pos_bp
                      << " duplicated position id: " << bp2id_map[pos_bp] << '\n';
            return false;
        }
        size_t next_id = bp2id_map.size();
        bp2id_map.insert({ pos_bp, next_id });
        // bp2id_map[pos_bp] = next_id;
        pos_bp_vec.push_back(pos_bp);
        pos_cm_vec.push_back(pos_cm);
        return true;
    }

    size_t
    get_size() const
    {
        return pos_bp_vec.size();
    }

    int
    get_chrom_id() const
    {
        return chrom_id;
    }

    uint32_t
    get_id(uint32_t bp_pos)
    {
        return bp2id_map.at(bp_pos);
    }

    uint32_t
    get_upper_bound_id(uint32_t bp_pos)
    {
        return std::distance(pos_bp_vec.begin(),
            std::upper_bound(pos_bp_vec.begin(), pos_bp_vec.end(), bp_pos));
    }

    uint32_t
    get_bp(uint32_t pid)
    {
        return pos_bp_vec[pid];
    }

    float
    get_cm(uint32_t pid)
    {
        // if (pid >= pos_cm_vec.size()) {
        //    std::cerr << "Stoped at pid: " << pid << '\n';
        //    print();
        //}
        return pos_cm_vec[pid];
    }
    std::vector<uint32_t> get_window_counts(float window_in_cm = 2.0);

    std::vector<region_label_t> get_gap_vector(
        float window_in_cm = 2.0, uint32_t min_snp_per_window = 2);

    // for debugging
    void print();

    void write_to_file(BGZF *fp);
    void read_from_file(BGZF *fp);

    bool
    is_equal(Positions &other)
    {
        return pos_bp_vec == other.pos_bp_vec && pos_cm_vec == other.pos_cm_vec
               && bp2id_map == other.bp2id_map;
    }
};

#endif
