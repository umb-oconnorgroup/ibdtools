#ifndef __positions_hpp__
#define __positions_hpp__

#include "common.hpp"
#include <algorithm>
#include <math.h>

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
    std::vector<uint32_t>
    get_window_counts(float window_in_cm = 2.0)
    {
        size_t max_cM = ceil(pos_cm_vec.back());
        // get multiples
        max_cM = ceil(max_cM / window_in_cm) * window_in_cm;

        // init vector
        std::vector<uint32_t> count_per_window;
        count_per_window.resize(max_cM / window_in_cm);
        std::fill(count_per_window.begin(), count_per_window.end(), 0);

        // count
        for (auto cm : pos_cm_vec)
            count_per_window[cm / window_in_cm] += 1;

        return count_per_window;
    }

    std::vector<region_label_t>
    get_gap_vector(float window_in_cm = 2.0, int min_snp_per_window = 2)
    {
        std::vector<uint32_t> count_per_2cm = get_window_counts(window_in_cm);

        std::vector<region_label_t> label_vec;

        size_t prev_label;
        for (uint32_t i = 0; i < count_per_2cm.size(); ++i) {

            uint32_t label = (count_per_2cm[i] >= min_snp_per_window);
            if (i == 0) {
                // first label and start
                label_vec.push_back({ 0, label });
                prev_label = label;
            } else {
                // extend
                if (label == prev_label)
                    continue;
                else {
                    // new label and start
                    float start_cm = window_in_cm * i;
                    uint32_t pid_s = std::distance(
                        pos_cm_vec.begin(), std::upper_bound(pos_cm_vec.begin(),
                                                pos_cm_vec.end(), start_cm));

                    //                                       pid_s = 52633
                    //                                       \
		    //                                       v
                    //
                    // ------x---x-----x---------------------x----x-----------------
                    //
                    //   |     keep      |    delete     |    keep       |
                    //   64              66              68              70
                    //
                    //   66 and 68 will map to the pid_s
                    //
                    // Need to change pid_s for regions maked for deletion
                    if (label == 0)
                        pid_s = pid_s - 1;

                    // std::cout << "i: " << i << "\t window_in_cm" << window_in_cm
                    //           << "\t start_cm " << start_cm << "\t pid_s: " << pid_s
                    //           << '\n';
                    //
                    label_vec.push_back({ pid_s, label });
                    prev_label = label;
                }
            }
        }

        return label_vec;
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
        write_element_to_file(chrom_id, fp);
        write_vector_to_file(pos_bp_vec, fp);
        write_vector_to_file(pos_cm_vec, fp);
    }
    void
    read_from_file(BGZF *fp)
    {
        read_element_from_file(chrom_id, fp);
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
