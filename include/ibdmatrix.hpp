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
    IbdMatrix(std::vector<uint16_t> &in_cm10x_vec, uint8_t in_is_hap_pair)
    {
        size_t array_size = in_cm10x_vec.size();
        // d^2 - d  = 2 *s;
        // d = sqrt(2 * s + 0.25) + 0.5
        this->d = (uint32_t) (sqrt(2 * array_size + 0.25) + 0.5);
        this->is_hap_pair_m = in_is_hap_pair;
        cm10x_vec.reserve(array_size);
        std::copy(
            in_cm10x_vec.begin(), in_cm10x_vec.end(), std::back_inserter(cm10x_vec));
    }

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

    void
    read_matrix_file(const char *matrix_fn)
    {
        BGZF *fp = bgzf_open(matrix_fn, "r");
        exit_on_false(fp != NULL, "", __FILE__, __LINE__);

        bgzf_mt(fp, 10, 256);

        std::string gzi(matrix_fn);
        gzi += ".gzi";
        if (std::filesystem::exists(gzi))
            exit_on_false(
                0 == bgzf_index_load(fp, matrix_fn, ".gzi"), "", __FILE__, __LINE__);

        read_element_from_file(is_hap_pair_m, fp);
        read_element_from_file(d, fp);
        read_vector_from_file(this->subpop_ids, fp);
        read_vector_from_file(cm10x_vec, fp);
        exit_on_false(get_array_size() == cm10x_vec.size(), "", __FILE__, __LINE__);

        bgzf_close(fp);
    }

    void
    write_matrix_file(const char *matrix_fn)
    {
        BGZF *fp = bgzf_open(matrix_fn, "w");
        exit_on_false(fp != NULL, "", __FILE__, __LINE__);

        bgzf_mt(fp, 10, 256);
        bgzf_index_build_init(fp);
        exit_on_false(get_array_size() == cm10x_vec.size(), "", __FILE__, __LINE__);

        write_element_to_file(is_hap_pair_m, fp);
        write_element_to_file(d, fp);
        write_vector_to_file(this->subpop_ids, fp);
        write_vector_to_file(cm10x_vec, fp);

        exit_on_false(
            0 == bgzf_index_dump(fp, matrix_fn, ".gzi"), "", __FILE__, __LINE__);
        bgzf_close(fp);
    }

    void
    add_matrix(const IbdMatrix &mat2)
    {
        std::cerr << "mat.d: " << d << " mat2.d: " << mat2.d << '\n';
        exit_on_false(d == mat2.d, "", __FILE__, __LINE__);
        std::transform(cm10x_vec.begin(), cm10x_vec.end(), mat2.cm10x_vec.begin(),
            cm10x_vec.begin(), std::plus<uint16_t>());
    }

    void
    calculate_total_from_ibdfile(IbdFile &ibdfile, bool use_hap_pair = false,
        const char *subpop_fn = NULL, float min_cM = 2.0, float max_cM = -1.0)
    {
        if (max_cM < min_cM) {
            max_cM = std::numeric_limits<float>::max();
        }
        MetaFile *meta = ibdfile.get_meta();
        exit_on_false(meta != NULL, "calculate total need info from the meta file",
            __FILE__, __LINE__);

        is_hap_pair_m = use_hap_pair;

        // for get cM from pos_id
        Positions &pos = meta->get_positions();

        // for allocation
        d = meta->get_samples().get_num_samples();
        if (use_hap_pair)
            d *= 2;

        cm10x_vec.resize(get_array_size(), 0);

        bool read_full;
        auto &in_vec = ibdfile.get_vec();
        uint32_t row, col, pid1, pid2;
        float cM;

        if (is_hap_pair()) {
            do {
                read_full = ibdfile.read_from_file();
                for (auto rec : in_vec) {
                    row = rec.get_sid1();
                    row <<= 1;
                    row += rec.get_hid1();
                    col = rec.get_sid2();
                    col <<= 1;
                    col += rec.get_hid2();
                    pid1 = rec.get_pid1();
                    pid2 = rec.get_pid2();
                    cM = pos.get_cm(pid2) - pos.get_cm(pid1);
                    if (cM < min_cM || cM >= max_cM)
                        continue;
                    at(row, col) += lround(cM * 10);
                }

            } while (read_full);

        } else {
            do {
                read_full = ibdfile.read_from_file();

                for (auto rec : in_vec) {
                    row = rec.get_sid1();
                    col = rec.get_sid2();

                    pid1 = rec.get_pid1();
                    pid2 = rec.get_pid2();
                    cM = pos.get_cm(pid2) - pos.get_cm(pid1);
                    if (cM < min_cM || cM >= max_cM)
                        continue;
                    at(row, col) += lround(cM * 10);
                }
            } while (read_full);
        }
    }

    // @win_size_in_10th_cm: 1 means 0.1cm window; 2 means 0.2cm ...
    void
    get_histogram(std::vector<size_t> &count_vec, uint16_t win_size_in_10th_cm = 1,
        MetaFile *meta = NULL, const char *subpop_fn = NULL, bool use_hap_pair = false)
    {
        auto max_element = std::max_element(cm10x_vec.begin(), cm10x_vec.end());
        uint16_t max_val = *max_element;
        size_t sz = max_val / win_size_in_10th_cm
                    + ((max_val % win_size_in_10th_cm) ? 1 : 0) + 1;

        // for subpopulation
        std::vector<uint8_t> subpop_v; // vector of 0's and 1's , 1 means sample is of a
                                       // subpopultion of interest
        if (subpop_fn != NULL) {
            exit_on_false(meta != NULL, "", __FILE__, __LINE__);
            meta->get_samples().get_subpop_vector(subpop_fn, subpop_v);

            // debgu
            // size_t sz = 0;
            // for (auto x : subpop_v)
            //     if (x != 0)
            //         sz++;
            // std::cerr << "subpop size: " << sz << '\n';
        }

        count_vec.resize(sz, 0);

        uint32_t sid1, sid2;
        for (size_t row = 1; row < d; row++)
            for (size_t col = 0; col < row; col++) {
                sid1 = use_hap_pair ? (row >> 1) : row;
                sid2 = use_hap_pair ? (col >> 1) : col;

                if (subpop_v.size() != 0 && (subpop_v[sid1] == 0 || subpop_v[sid2] == 0))
                    continue;

                count_vec[at(row, col) / win_size_in_10th_cm]++;
            }
    }

    void
    subset_matrix(const std::vector<uint32_t> &in_subpop_ids)
    {
        // check if already subsetted
        bool not_subsetted = this->subpop_ids.size() == 0;
        bool input_equal_to_member = this->subpop_ids == in_subpop_ids;
        bool is_valid = (not_subsetted || input_equal_to_member);
        exit_on_false(is_valid,
            "Error: subsetting already performed; new subsetting not allowed", __FILE__,
            __LINE__);

        // copy args to member variable
        this->subpop_ids.reserve(in_subpop_ids.size());
        std::copy(in_subpop_ids.begin(), in_subpop_ids.end(),
            std::back_inserter(this->subpop_ids));

        // sort member variable
        std::sort(this->subpop_ids.begin(), this->subpop_ids.end());

        // validate subpop_ids
        auto new_d = this->subpop_ids.size();
        auto new_array_size = (new_d - 1) * new_d / 2;
        bool new_d_lt2 = new_d < 2;
        bool has_repeat
            = std::adjacent_find(this->subpop_ids.begin(), this->subpop_ids.end())
              != this->subpop_ids.end();
        bool too_large = this->subpop_ids.back() >= d;
        exit_on_false(!new_d_lt2, "subpop_ids should be at least 2", __FILE__, __LINE__);
        exit_on_false(!has_repeat, "subpop_ids has repeated ids", __FILE__, __LINE__);
        exit_on_false(!has_repeat, "subpop_ids has repeated ids", __FILE__, __LINE__);
        exit_on_false(!too_large, "subpop_ids id out of range", __FILE__, __LINE__);

        // consolidate the matrix
        for (size_t i = 0; i < new_array_size; i++) {
            // get row, col index in new array
            size_t row = sqrt(2 * i + 0.25) + 0.5;
            size_t col = i - row * (row - 1) / 2;
            // get row, col index in old array
            size_t row_old = this->subpop_ids[row];
            size_t col_old = this->subpop_ids[col];
            // copy
            this->at(row, col) = this->at(row_old, col_old);
        }

        // update d and shrink array
        this->d = new_d;
        this->cm10x_vec.resize(new_array_size);
        this->cm10x_vec.shrink_to_fit();
    }

    // keep any element x if low <= x < upper and set other element to 0
    void
    filter_matrix(uint16_t low, uint16_t upper)
    {
        for_each(cm10x_vec.begin(), cm10x_vec.end(), [low, upper](auto &x) {
            if (x < low || x >= upper)
                x = 0;
        });
    }

    void
    print_to_ostream(std::ostream &os, uint32_t r_first = 0, uint32_t r_last = 10,
        uint32_t c_first = 0, uint32_t c_last = 10)
    {
        exit_on_false(
            0 <= r_first && r_first < r_last && r_last < d, "", __FILE__, __LINE__);
        exit_on_false(
            0 <= c_first && c_first < c_last && c_last < d, "", __FILE__, __LINE__);

        for (uint32_t row = r_first; row <= r_last; row++) {
            for (uint32_t col = c_first; col <= c_last; col++) {
                if (row <= col)
                    os << "--\t";
                else
                    os << at(row, col) / 10.0 << '\t';
            }
            os << '\n';
        }
    }

    void
    print_subpop_ids_to_ostream(std::ostream &os)
    {
        for (auto id : this->subpop_ids) {
            os << id << '\n';
        }
    }
};

#endif
