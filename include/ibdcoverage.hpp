#ifndef __ibdcoverage_hpp__
#define __ibdcoverage_hpp__
#include "common.hpp"
#include "ibdfile.hpp"
#include <cstdint>
#include <filesystem>
#include <fstream>
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
    std::vector<uint8_t> subpop_v; // vector of 0's and 1's , 1 means sample is of a
                                   // subpopultion of interest

    // for parallelization
    using Iter = std::vector<ibd_rec1_t>::iterator;
    struct group_t {
        Iter first, last;
        std::vector<size_t> grp_cnt_vec;
    };

    // this is can be larger when used for parallelization
    size_t max_groups{ 1 };
    std::vector<group_t> grps_vec;

  public:
    IbdCoverage() {}
    IbdCoverage(const char *ibd_fn, const char *meta_fn, float win_in_cM = 1.0,
        size_t max_rec_ram = 1 * 1024 * 1024, const char *subpop_fn = NULL)
        : window_in_cM(win_in_cM)
    {
        // read meta file
        BGZF *fp = bgzf_open(meta_fn, "r");
        assert(fp != NULL);
        meta.read_from_file(fp);
        bgzf_close(fp);

        // parse subpop
        if (subpop_fn != NULL) {
            auto &samples = meta.get_samples();
            subpop_v.resize(samples.get_num_samples(), 0);
            std::ifstream ifs(subpop_fn);
            std::string line;
            while (getline(ifs, line, '\n')) {
                uint32_t sid = samples.get_id(line);
                subpop_v[sid] = 1;
            }
        }

        // prepare ibdfile object
        in = IbdFile(ibd_fn, &meta, max_rec_ram);

        // prepare cm_vec
        auto &pos = meta.get_positions();
        float max_cM = pos.get_cm(pos.get_size() - 1);
        size_t num_win = max_cM / win_in_cM + 1;
        cm_vec.reserve(num_win);

        for (size_t i = 0; i <= num_win; i++)
            cm_vec.push_back(win_in_cM * i);

        // init count vec;
        count_vec.resize(cm_vec.size(), 0);

        // grp_vec allocation
        grps_vec.resize(max_groups);

        // grp_counts_vec allocation
        for (auto &grp : grps_vec) {
            grp.grp_cnt_vec.resize(cm_vec.size(), 0);
        }

        assert(subpop_v.size() == meta.get_samples().get_num_samples()
               || subpop_v.size() == 0);
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
            read_full = in.read_from_file();

            // divide into groups for parallelizating
            divide_to_groups();

            // calculate group coverage (can use parallelizing algorithm
            for_each(grps_vec.begin(), grps_vec.end(),
                [this](group_t &grp) mutable { this->calculate_grp_coverage(grp); });

            // add each groups results to count_vec
            for (auto &grp : grps_vec)
                // this can be parallelized
                std::transform(count_vec.begin(), count_vec.end(),
                    grp.grp_cnt_vec.begin(), count_vec.begin(), std::plus<size_t>());

            total_rec_processed += vec.size();

            // debug
            // summary(std::cout, 300);
            // print_group_state();

            // exit(1);

        } while (read_full);

        in.close();
    }

  private:
    void
    divide_to_groups()
    {
        auto &vec = in.get_vec();
        size_t step = vec.size() / max_groups;

        for (size_t i = 0; i < max_groups; i++) {
            // first
            grps_vec[i].first = vec.begin() + step * i;

            // last
            if (i > 0) {
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
        grp.grp_cnt_vec.clear();
        grp.grp_cnt_vec.resize(cm_vec.size(), 0);

        auto &pos = meta.get_positions();

        for (Iter it = grp.first; it < grp.last; it++) {

            // when the subpop vector is in use, if any of the sample in the ibd rec
            // is not in subpop, skip counting coverage for this vector
            if (subpop_v.size() != 0
                && (subpop_v[it->get_sid1()] == 0 || subpop_v[it->get_sid2()] == 0))
                continue;

            //  upper_bound:    a[i-1] <= x < a[i]
            //  lower_bound:    a[i-1] < x <= a[i]

            auto first = std::distance(cm_vec.begin(),
                             std::upper_bound(cm_vec.begin(), cm_vec.end(),
                                 pos.get_cm(it->get_pid1())))
                         - 1;
            auto last = std::distance(
                cm_vec.begin(), std::lower_bound(cm_vec.begin() + first, cm_vec.end(),
                                    pos.get_cm(it->get_pid2())));

            for (; first < last; first++)
                grp.grp_cnt_vec[first]++;
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

    // if print_n = 0, print all records
    void
    summary(std::ostream &os, size_t print_n = 10)
    {
        os << "# Window size (cM): " << window_in_cM << '\n';
        os << "# Records processed: " << total_rec_processed << '\n';
        os << "# Coverage:\n";
        os << "cM\tCount\n";
        if (print_n == 0)
            print_n = cm_vec.size();
        for (size_t i = 0; i < print_n && i < cm_vec.size(); i++)
            os << cm_vec[i] << '\t' << count_vec[i] << '\n';
        if (print_n < cm_vec.size())
            os << "...\n";
    }

    void
    print_group_state()
    {
        std::cout << "Group state: \n";
        for (auto &grp : grps_vec) {
            std::cout << std::distance(in.get_vec().begin(), grp.first) << ", "
                      << std::distance(in.get_vec().begin(), grp.last) << '\n';
        }
    }

    std::vector<float> &
    get_cm_vec()
    {
        return cm_vec;
    }

    std::vector<size_t> &
    get_count_vec()
    {
        return count_vec;
    }
};
#endif
