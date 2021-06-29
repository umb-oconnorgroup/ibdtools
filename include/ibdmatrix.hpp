#ifndef __ibdmatrix_hpp__
#define __ibdmatrix_hpp__
#include "common.hpp"
#include "ibdfile.hpp"
#include "metafile.hpp"
#include "positions.hpp"
#include <filesystem>

class IbdMatrix
{
    std::vector<uint16_t> cm10x_vec; // convert float to int16_t by round(10xcM);
    uint32_t nsam;

  public:
    IbdMatrix() {}

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
        size_t sz = nsam;
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
        return nsam;
    }

    void
    get_matrix_index(size_t arr_index, uint32_t &row, uint32_t &col)
    {
        // r* (r - 1) / 2 + c = y;
        // 0 <= c <= r-1
        // sqrt(2*y + 0.25) + 0.5 >= r >= sqrt(2*y + 2.25) - 0.5
        //

        row = sqrtl(2 * arr_index + 0.25) + 0.5;

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
        assert(fp != NULL);

        bgzf_mt(fp, 10, 256);

        std::string gzi(matrix_fn);
        gzi += ".gzi";
        std::filesystem::exists(gzi);
        assert(0 == bgzf_index_load(fp, matrix_fn, ".gzi"));

        read_element_from_file(nsam, fp);
        read_vector_from_file(cm10x_vec, fp);
        assert(get_array_size() == cm10x_vec.size());

        bgzf_close(fp);
    }

    void
    write_matrix_file(const char *matrix_fn)
    {
        BGZF *fp = bgzf_open(matrix_fn, "w");
        assert(fp != NULL);

        bgzf_mt(fp, 10, 256);

        bgzf_index_build_init(fp);

        assert(get_array_size() == cm10x_vec.size());
        write_element_to_file(nsam, fp);
        write_vector_to_file(cm10x_vec, fp);

        assert(0 == bgzf_index_dump(fp, matrix_fn, ".gzi"));
        bgzf_close(fp);
    }

    void
    calculate_total_from_ibdfile(IbdFile &ibdfile)
    {
        MetaFile *meta = ibdfile.get_meta();
        assert(meta != NULL && "calculate total need info from the meta file");

        // for get cM from pos_id
        Positions &pos = meta->get_positions();

        // for allocation
        nsam = meta->get_samples().get_num_samples();
        cm10x_vec.resize(get_array_size());

        bool read_full;
        auto &in_vec = ibdfile.get_vec();
        uint32_t row, col, pid1, pid2;
        do {
            read_full = ibdfile.read_from_file();
            for (auto rec : in_vec) {
                row = rec.get_sid1();
                col = rec.get_sid2();
                pid1 = rec.get_pid1();
                pid2 = rec.get_pid2();
                at(row, col) += lround((pos.get_cm(pid2) - pos.get_cm(pid1)) * 10);
            }
        } while (read_full);
    }

    void
    print_to_ostream(std::ostream &os, uint32_t r_first = 0, uint32_t r_last = 10,
        uint32_t c_first = 0, uint32_t c_last = 10)
    {
        assert(0 <= r_first && r_first < r_last && r_last < nsam);
        assert(0 <= c_first && c_first < c_last && c_last < nsam);

        for (uint32_t row = r_first; row < r_last; row++) {
            for (uint32_t col = c_first; col < c_last; col++) {
                if (row <= col)
                    os << "--\t";
                else
                    os << at(row, col) / 10.0 << '\t';
            }
            os << '\n';
        }
    }
};

#endif
