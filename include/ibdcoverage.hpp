#ifndef __ibdcoverage_hpp__
#define __ibdcoverage_hpp__
#include "common.hpp"
#include "metafile.hpp"
#include <iterator>
#include <string>
#include <vector>

class IbdFile;
class MetaFile;

class IbdCoverage
{
    // output variable
    std::vector<float> cm_vec;
    std::vector<size_t> count_vec; // each element corresponds to one element in cm_vec
    float window_in_cM;
    size_t total_rec_processed{ 0 };

    // inputs
    std::unique_ptr<IbdFile> in;
    std::unique_ptr<MetaFile> meta;
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
    IbdCoverage();
    IbdCoverage(const char *ibd_fn, const char *meta_fn, float win_in_cM = 1.0,
        size_t max_rec_ram = 1 * 1024 * 1024, const char *subpop_fn = NULL);

    // 1. read chunk to memory
    // 2. sort
    // 3. divide by sample_pair group
    // 4. use parallele algorithm to calculate covarage.
    void calculate_coverage();

  private:
    void divide_to_groups();

    void calculate_grp_coverage(group_t &grp);

  public:
    void write_to_file(const char *ibd_cov_file);

    void read_from_file(const char *ibd_cov_file);

    // if print_n = 0, print all records
    void summary(std::ostream &os, size_t print_n = 10);

    void print_group_state();

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
