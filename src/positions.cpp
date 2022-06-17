
#include "positions.hpp"
#include "common.hpp"
#include <math.h>

std::vector<uint32_t>
Positions::get_window_counts(float window_in_cm)
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
Positions::get_gap_vector(float window_in_cm, uint32_t min_snp_per_window)
{
    std::vector<uint32_t> count_per_2cm = get_window_counts(window_in_cm);

    std::vector<region_label_t> label_vec;

    size_t prev_label;
    for (uint32_t i = 0; i < count_per_2cm.size(); ++i) {

        uint32_t label = (size_t) (count_per_2cm[i] >= min_snp_per_window);
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
                uint32_t pid_s = std::distance(pos_cm_vec.begin(),
                    std::upper_bound(pos_cm_vec.begin(), pos_cm_vec.end(), start_cm));

                /*                                       pid_s = 52633
                                                         \
                                                 v

                   ------x---x-----x---------------------x----x-----------------

                     |     keep      |    delete     |    keep       |
                     64              66              68              70

                     66 and 68 will map to the pid_s
                */
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
Positions::print()
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
Positions::write_to_file(BGZF *fp)
{
    write_element_to_file(chrom_id, fp);
    write_vector_to_file(pos_bp_vec, fp);
    write_vector_to_file(pos_cm_vec, fp);
}
void
Positions::read_from_file(BGZF *fp)
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