#ifndef __positions_hpp
#define __positions_hpp

#include "common.hpp"

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

    int
    get_chrom_id()
    {
        return chrom_id;
    }

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
    get_bp(uint32_t pid)
    {
        return pos_bp_vec[pid];
    }

    float
    get_cm(uint32_t pid)
    {
        return pos_cm_vec[pid];
    }

    // for debugging
    void
    print()
    {
        for (size_t i = 0; i < pos_bp_vec.size(); i++) {
            std::cout << "bp and cM: " << pos_bp_vec[i] << '\t' << pos_cm_vec[i] << '\n';
        }
        for (auto m : bp2id_map) {
            std::cout << "bp2id_map " << m.first << '\t' << m.second << '\n';
        }
        std::cout << pos_bp_vec.size() << ' ' << bp2id_map.size() << '\n';
    }

    void
    write_to_file(BGZF *fp)
    {
        write_vector_to_file(pos_bp_vec, fp);
        write_vector_to_file(pos_cm_vec, fp);
    }
    void
    read_from_file(BGZF *fp)
    {
        read_vector_from_file(pos_bp_vec, fp);
        read_vector_from_file(pos_cm_vec, fp);

        // fill the bp2id_map
        bp2id_map.clear();
        for (size_t i = 0; i < pos_bp_vec.size(); i++) {
            size_t next_id = bp2id_map.size();
            bp2id_map[pos_bp_vec[i]] = next_id;
        }
    }

    bool
    is_equal(Positions &other)
    {
        return pos_bp_vec == other.pos_bp_vec && pos_cm_vec == other.pos_cm_vec
               && bp2id_map == other.bp2id_map;
    }
};

#endif
