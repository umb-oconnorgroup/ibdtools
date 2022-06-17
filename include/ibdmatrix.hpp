#ifndef __ibdmatrix_hpp__
#define __ibdmatrix_hpp__
#include "common.hpp"
#include "ibdfile.hpp"
#include "metafile.hpp"
#include "positions.hpp"
#include <filesystem>
#include <iostream>

class IbdMatrix
{
    std::vector<uint16_t> cm10x_vec; // convert float to int16_t by round(10xcM);
    uint32_t d; // d = nsam if is not hap_pair matrix; otherwise d = 2 * nsam;
    uint8_t is_hap_pair_m;
    std::vector<uint32_t> subpop_ids;

  public:
    IbdMatrix() {}

    /*Construct mat from a vector*/
    IbdMatrix(std::vector<uint16_t> &in_cm10x_vec, uint8_t in_is_hap_pair);

    // Lower trangular  row > col;
    //  -  -  -  -  -
    //  0  -  -  -  -
    //  1  2  -  -  -
    //  3  4  5  -  -
    //  6  7  8  9  -
    //
    size_t
    get_array_size()
    {
        size_t sz = d;
        return sz * (sz - 1) / 2;
    }

    size_t
    get_arr_index(size_t row, size_t col)
    {
        return row * (row - 1) / 2 + col;
    }

    uint32_t
    get_num_samples()
    {
        return d;
    }

    bool
    is_hap_pair()
    {
        return is_hap_pair_m != 0;
    }

    void
    get_matrix_index(size_t arr_index, uint32_t &row, uint32_t &col)
    {
        // r* (r - 1) / 2 + c = y;
        // 0 <= c <= r-1
        // sqrt(2*y + 0.25) + 0.5 >= r >= sqrt(2*y + 2.25) - 0.5
        //

        row = sqrt(2 * arr_index + 0.25) + 0.5;

        // overflow if uint32_t * uint32_t. Need to convert to size_t
        // col = arr_index -  row * (row - 1) / 2;
        size_t tmp = row;
        col = arr_index - tmp * (tmp - 1) / 2;
    }

    uint16_t &
    at(size_t row, size_t col)
    {
        return cm10x_vec[row * (row - 1) / 2 + col];
    }

    uint16_t &
    at(size_t arr_index)
    {
        return cm10x_vec[arr_index];
    }

    void read_matrix_file(const char *matrix_fn);

    void write_matrix_file(const char *matrix_fn);

    void add_matrix(const IbdMatrix &mat2);

    void calculate_total_from_ibdfile(IbdFile &ibdfile, bool use_hap_pair = false,
        const char *subpop_fn = NULL, float min_cM = 2.0, float max_cM = -1.0);
    // @win_size_in_10th_cm: 1 means 0.1cm window; 2 means 0.2cm ...
    void get_histogram(std::vector<size_t> &count_vec, uint16_t win_size_in_10th_cm = 1,
        MetaFile *meta = NULL, const char *subpop_fn = NULL, bool use_hap_pair = false);
    void subset_matrix(const std::vector<uint32_t> &in_subpop_ids);

    // keep any element x if low <= x < upper and set other element to 0
    void filter_matrix(uint16_t low, uint16_t upper);

    void print_to_ostream(std::ostream &os, uint32_t r_first = 0, uint32_t r_last = 10,
        uint32_t c_first = 0, uint32_t c_last = 10);

    void print_subpop_ids_to_ostream(std::ostream &os);
};

#endif
