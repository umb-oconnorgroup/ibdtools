#ifndef __ibdcoverage_hpp__
#define __ibdcoverage_hpp__
#include "common.hpp"
#include "ibdfile.hpp"
#include <cstdint>
#include <filesystem>
#include <iostream>

class IbdCoverage
{
    // output variable
    std::vector<float> cm_vec;
    std::vector<size_t> count_vec; // each element corresponds to one element in cm_vec
    float window_in_cM;
    size_t total_rec_processed{ 0 };

    // inputs
    IbdFile in;
    MetaFile meta;

    // for parallelization
    using Iter = std::vector<ibd_rec1_t>::iterator;
    struct group_t {
        Iter first, last;
        std::vector<size_t> grp_cnt_vec;
    };
    size_t max_groups{ 64 };
    std::vector<group_t> grps_vec;

  public:
    IbdCoverage() {}
    IbdCoverage(const char *ibd_fn, const char *meta_fn, float win_in_cM = 1.0,
        size_t max_rec_ram = 1 * 1024 * 1024)
        : window_in_cM(win_in_cM)
    {

        // read meta file
        BGZF *fp = bgzf_open(meta_fn, "r");
        assert(fp != NULL);
        meta.read_from_file(fp);
        bgzf_close(fp);

        // prepare ibdfile object
        in = IbdFile(ibd_fn, &meta, max_rec_ram);

        // prepare cm_vec
        auto &pos = meta.get_positions();
        float max_cM = pos.get_cm(pos.get_size() - 1);
        size_t num_win = max_cM / win_in_cM + 1;
        cm_vec.reserve(num_win);

        for (size_t i = 0; i < num_win; i++)
            cm_vec.push_back(win_in_cM * i);

        // init count vec;
        count_vec.resize(cm_vec.size(), 0);

        // grp_vec allocation
        grps_vec.reserve(max_groups);

        // grp_counts_vec allocation
        for (auto &grp : grps_vec) {
            grp.grp_cnt_vec.resize(cm_vec.size(), 0);
        }
    }

    // 1. read chunk to memory
    // 2. sort
    // 3. divide by sample_pair group
    // 4. use parallele algorithm to calculate covarage.
    void
    calculate_coverage()
    {
        in.open("r");
        auto &vec = in.get_vec();
        bool read_full;
        do {
            // read ibd into memory
            in.read_from_file();

            // divide into groups for parallelizating
            divide_to_groups();

            // calculate group coverage (can use parallelizing algorithm
            for_each(grps_vec.begin(), grps_vec.end(),
                std::bind(&IbdCoverage::calculate_grp_coverage, this));

            // add each groups results to count_vec
            for (auto &grp : grps_vec)
                // this can be parallelized
                std::transform(count_vec.begin(), count_vec.end(),
                    grp.grp_cnt_vec.begin(), count_vec.begin(), std::plus<size_t>());

            total_rec_processed += vec.size();

        } while (read_full);

        in.close();
    }

  private:
    void
    divide_to_groups()
    {
        auto vec = in.get_vec();
        size_t step = vec.size() / max_groups;

        for (size_t i = 0; i < max_groups; i++) {
            grps_vec[i].first = vec.begin() + step * i;
            // first
            // find group boundaries
            while (grps_vec[i].first != vec.begin()
                   && grps_vec[i].first->is_same_pair(*(grps_vec[i].first - 1)))
                grps_vec[i].first--;

            // last
            if (i > 1) {
                grps_vec[i - 1].last = grps_vec[i].first;
            }
            if (i == max_groups - 1)
                grps_vec[i].last = vec.end();
        }
    }

    void
    calculate_grp_coverage(group_t &grp)
    {
        // clear counts to zero
        grp.grp_cnt_vec.resize(cm_vec.size(), 0);

        auto &pos = meta.get_positions();
        uint32_t start_cm, end_cm;

        for (Iter it = grp.first; it < grp.last; it++) {
            start_cm = pos.get_cm(it->get_pid1());
            end_cm = pos.get_cm(it->get_pid2());

            auto cm_vec_id = std::distance(cm_vec.begin(),
                                 std::upper_bound(cm_vec.begin(), cm_vec.end(),
                                     pos.get_cm(it->get_pid1())))
                             - 1;
            for (; cm_vec_id < cm_vec.size(); cm_vec_id++)
                grp.grp_cnt_vec[cm_vec_id]++;
        }
    }

  public:
    void
    write_to_file(const char *ibd_cov_file)
    {
        BGZF *fp = bgzf_open(ibd_cov_file, "w");
        assert(fp != NULL);
        bgzf_mt(fp, 10, 256);
        bgzf_index_build_init(fp);

        write_element_to_file(window_in_cM, fp);
        write_element_to_file(total_rec_processed, fp);
        write_vector_to_file(cm_vec, fp);
        write_vector_to_file(count_vec, fp);

        assert(0 == bgzf_index_dump(fp, ibd_cov_file, ".gzi"));
        bgzf_close(fp);
    }

    void
    read_from_file(const char *ibd_cov_file)
    {
        BGZF *fp = bgzf_open(ibd_cov_file, "r");
        assert(fp != NULL);
        bgzf_mt(fp, 10, 256);
        std::string gzi = ibd_cov_file;
        if (std::filesystem::exists(gzi + ".gzi")) {
            assert(0 == bgzf_index_load(fp, ibd_cov_file, ".gzi"));
        }

        read_element_from_file(window_in_cM, fp);
        read_element_from_file(total_rec_processed, fp);
        read_vector_from_file(cm_vec, fp);
        read_vector_from_file(count_vec, fp);

        bgzf_close(fp);
    }
};
#endif
