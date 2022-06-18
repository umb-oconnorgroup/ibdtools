#ifndef __common_hpp__
#define __common_hpp__

#include <chrono>
#include <htslib/hts.h>
#include <iostream>
#include <vector>

// from htslib bgzf.c
typedef struct {
    uint64_t uaddr; // offset w.r.t. uncompressed data
    uint64_t caddr; // offset w.r.t. compressed data
} bgzidx1_t;

// from htslib bgzf.c
// struct bgzidx_t {
//     int noffs, moffs;     // the size of the index, n:used, m:allocated
//     bgzidx1_t *offs;      // offsets
//     uint64_t ublock_addr; // offset of the current block (uncompressed data)
// };

void exit_on_false(bool condition, const char *message, const char *file, int lineno);
void exit_with_message(const char *message);

struct region_label_t {
    uint32_t pid_s : 24; // region start position id
    uint32_t label : 8;
    friend bool
    operator==(const region_label_t &x, const region_label_t &y)
    {
        return x.pid_s == y.pid_s && x.label == y.label;
    }
};

void add_exclusion_range(
    std::vector<region_label_t> &label_v, uint32_t pid_left, uint32_t pid_right);

struct ibd_rec2_t {
    // data member
    uint32_t sid1;
    uint32_t sid2;
    uint32_t pid1;
    uint32_t pid2;
    uint8_t hid1;
    uint8_t hid2;

    // getter
    uint32_t
    get_sid1() const
    {
        return sid1;
    }
    uint8_t
    get_hid1() const

    {
        return hid1;
    }
    uint32_t
    get_sid2() const
    {
        return sid2;
    }
    uint8_t
    get_hid2() const
    {
        return hid2;
    }
    uint32_t
    get_pid1() const
    {
        return pid1;
    }
    uint32_t
    get_pid2() const
    {
        return pid2;
    }

    // setter
    void
    set_sid1(uint32_t sid1_32)
    {
        sid1 = sid1_32;
    }
    void
    set_hid1(uint8_t hid1_32)
    {
        hid1 = hid1_32;
    }
    void
    set_sid2(uint32_t sid2_32)
    {
        sid2 = sid2_32;
    }
    void
    set_hid2(uint8_t hid2_32)
    {
        hid2 = hid2_32;
    }
    void
    set_pid1(uint32_t pid1_32)
    {
        pid1 = pid1_32;
    }
    void
    set_pid2(uint32_t pid2_32)
    {
        pid2 = pid2_32;
    }

    void
    maximize()
    {
        sid1 = 0xffffffff;
        sid2 = 0xffffffff;
        pid1 = 0xffffffff;
        pid2 = 0xffffffff;
        hid1 = 0xff;
        hid2 = 0xff;
    }

    // less than operator
    friend bool
    operator<(const ibd_rec2_t &rec1, const ibd_rec2_t &rec2)
    {
        uint32_t a, b;
        a = rec1.get_sid1();
        b = rec2.get_sid1();
        if (a != b)
            return a < b;
        a = rec1.get_sid2();
        b = rec2.get_sid2();
        if (a != b)
            return a < b;
        a = rec1.get_pid1();
        b = rec2.get_pid1();
        if (a != b)
            return a < b;
        a = rec1.get_pid2();
        b = rec2.get_pid2();
        if (a != b)
            return a < b;
        a = rec1.get_hid1();
        b = rec2.get_hid1();
        if (a != b)
            return a < b;
        a = rec1.get_hid2();
        b = rec2.get_hid2();
        return a < b;
    }

    // for debug
    void
    print()
    {
        // add + to uint8_t to make std::cout treat it as integer instead of char
        std::cout << sid1 << '\t' << +hid1 << '\t' << sid2 << '\t' << +hid2 << '\t'
                  << pid1 << '\t' << pid2 << '\n';
    }

    // equal operator
    friend bool
    operator==(const ibd_rec2_t &rec1, const ibd_rec2_t &rec2)
    {
        return rec1.sid1 == rec2.sid1 && rec1.sid2 == rec2.sid2 && rec1.pid1 == rec2.pid1
               && rec1.pid2 == rec2.pid2 && rec1.hid1 == rec2.hid1
               && rec1.hid2 == rec2.hid2;
    }

    bool
    is_same_pair(const ibd_rec2_t &rec2)
    {
        return sid1 == rec2.sid1 && sid2 == rec2.sid2;
    }
};

// packed version of ibd_rec2_t
struct __attribute__((packed)) ibd_rec1_t {

    // data member
    // left -> right: low to high bits
    // up -> down, low bytes to high bytes
    uint16_t hid2 : 1, hid1 : 1, pid2_l : 14;
    uint64_t pid2_h : 6, pid1 : 20, sid2 : 19, sid1 : 19;

    ibd_rec1_t() {}

    ibd_rec1_t(ibd_rec2_t &r2)
    {
        // this works but need to make sure that the assigned value will not overflow the
        // bit
        // fields
        *((uint16_t *) this) = (r2.pid2 & 0x3fff) << 2;
        *((uint16_t *) this) |= ((r2.hid1 << 1) + r2.hid2);

        *(uint64_t *) ((uint16_t *) this + 1) = r2.sid1;
        *(uint64_t *) ((uint16_t *) this + 1) <<= 19;
        *(uint64_t *) ((uint16_t *) this + 1) |= r2.sid2;
        *(uint64_t *) ((uint16_t *) this + 1) <<= 20;
        *(uint64_t *) ((uint16_t *) this + 1) |= r2.pid1;
        *(uint64_t *) ((uint16_t *) this + 1) <<= 6;
        *(uint64_t *) ((uint16_t *) this + 1) |= (r2.pid2 >> 14);

        /*
        // this works but slower
        sid1 = r2.sid1;
        hid1 = r2.hid1;
        sid2 = r2.sid2;
        hid2 = r2.hid2;
        pid1 = r2.pid1;
        pid2_l = (r2.pid2 & 0x3fff);
        pid2_h = (r2.pid2 >> 14);
    */

        /*  // for debug
        if(!is_equal(*this, r2))
        {
            this->print();
            r2.print();
        }
        exit_on_false(is_equal(*this, r2), "", __FILE__, __LINE__);
        */
    }
    // getter
    uint32_t
    get_sid1() const
    {
        return sid1;
    }
    uint8_t
    get_hid1() const
    {
        return hid1;
    }
    uint32_t
    get_sid2() const
    {
        return sid2;
    }
    uint8_t
    get_hid2() const
    {
        return hid2;
    }
    uint32_t
    get_pid1() const
    {
        return pid1;
    }
    uint32_t
    get_pid2() const
    {
        return ((uint32_t) pid2_h << 14) + pid2_l;
    }

    // setter
    void
    set_sid1(uint32_t sid1_32)
    {
        sid1 = sid1_32;
    }
    void
    set_hid1(uint8_t hid1_32)
    {
        hid1 = hid1_32;
    }
    void
    set_sid2(uint32_t sid2_32)
    {
        sid2 = sid2_32;
    }
    void
    set_hid2(uint8_t hid2_32)
    {
        hid2 = hid2_32;
    }
    void
    set_pid1(uint32_t pid1_32)
    {
        pid1 = pid1_32;
    }
    void
    set_pid2(uint32_t pid2_32)
    {
        pid2_h = pid2_32 >> 14;
        pid2_l = pid2_32 & 0x3fff;
    }

    void
    maximize()
    {
        *((uint64_t *) this) = ~(uint64_t) 0;
        ((uint16_t *) this)[4] = ~(uint16_t) 0;
    }

    // less than operator
    // Sort by fields in the order of: sid1, sid2, pid1, pid2, hid1, hid2
    friend bool
    operator<(const ibd_rec1_t &rec1, const ibd_rec1_t &rec2)
    {
        if (*((uint64_t *) ((uint16_t *) &rec1 + 1))
            != *((uint64_t *) ((uint16_t *) &rec2 + 1)))
            return *((uint64_t *) ((uint16_t *) &rec1 + 1))
                   < *((uint64_t *) ((uint16_t *) &rec2 + 1));

        return ((uint16_t *) &rec1)[0] < ((uint16_t *) &rec2)[0];
    }

    // for debug
    friend bool
    is_less(const ibd_rec1_t &rec1, const ibd_rec1_t &rec2)
    {

        uint32_t a, b;
        a = rec1.get_sid1();
        b = rec2.get_sid1();
        if (a != b)
            return a < b;
        a = rec1.get_sid2();
        b = rec2.get_sid2();
        if (a != b)
            return a < b;
        a = rec1.get_pid1();
        b = rec2.get_pid1();
        if (a != b)
            return a < b;
        a = rec1.get_pid2();
        b = rec2.get_pid2();
        if (a != b)
            return a < b;
        a = rec1.get_hid1();
        b = rec2.get_hid1();
        if (a != b)
            return a < b;
        a = rec1.get_hid2();
        b = rec2.get_hid2();
        return a < b;
    }

    // equal operator
    friend bool
    operator==(const ibd_rec1_t &rec1, const ibd_rec1_t &rec2)
    {
        return (*((uint64_t *) &rec1) == *((uint64_t *) &rec2))
               && ((uint16_t *) &rec1)[4] == ((uint16_t *) &rec2)[4];
    }

    // for debug
    friend bool
    is_equal(const ibd_rec1_t &rec1, const ibd_rec1_t &rec2)
    {

        return rec1.get_sid1() == rec2.get_sid1() && rec1.get_sid2() == rec2.get_sid2()
               && rec1.get_pid1() == rec2.get_pid1()
               && rec1.get_pid2() == rec2.get_pid2()
               && rec1.get_hid1() == rec2.get_hid1()
               && rec1.get_hid2() == rec2.get_hid2();
    }

    friend bool
    is_equal(const ibd_rec1_t &rec1, const ibd_rec2_t &rec2)
    {

        return rec1.get_sid1() == rec2.get_sid1() && rec1.get_sid2() == rec2.get_sid2()
               && rec1.get_pid1() == rec2.get_pid1()
               && rec1.get_pid2() == rec2.get_pid2()
               && rec1.get_hid1() == rec2.get_hid1()
               && rec1.get_hid2() == rec2.get_hid2();
    }

    bool
    is_same_pair(const ibd_rec1_t &rec2)
    {
        return sid1 == rec2.sid1 && sid2 == rec2.sid2;
    }

    // for debug
    void
    print()
    {
        std::cout << sid1 << '\t' << hid1 << '\t' << sid2 << '\t' << hid2 << '\t' << pid1
                  << '\t' << ((pid2_h << 14) + pid2_l) << '\n';
    }
};

using ibd_rec_cmp_t = std::function<bool(const ibd_rec1_t &, const ibd_rec1_t &)>;

/*
struct IbdComparatorHapPair {
    bool
    operator()(const ibd_rec1_t &r1, const ibd_rec1_t &r2)
    {
        if (r1.get_sid1() != r2.get_sid1())
            return r1.get_sid1() < r2.get_sid1();
        if (r1.get_hid1() != r2.get_hid1())
            return r1.get_hid1() < r2.get_hid1();
        if (r1.get_sid2() != r2.get_sid2())
            return r1.get_sid2() < r2.get_sid2();
        if (r1.get_hid2() != r2.get_hid2())
            return r1.get_hid2() < r2.get_hid2();
        if (r1.get_pid1() != r2.get_pid1())
            return r1.get_pid1() < r2.get_pid1();
        return r1.get_pid2() < r2.get_pid2();
    }
};
*/

// Adapted from cpp_high_perforamnce_2e_bjorn_andrist
class ScopedTimer
{
    using ClockType = std::chrono::steady_clock;
    bool human_readable{ false };

  private:
    const char *function_name_{};
    const ClockType::time_point start_{};

  public:
    ScopedTimer(const char *func, bool human_readable = false);
    ScopedTimer(ScopedTimer &&) = delete;
    ScopedTimer &operator=(const ScopedTimer &) = delete;
    ScopedTimer &operator=(ScopedTimer &&) = delete;
    ~ScopedTimer();
};

template <typename T> void write_vector_to_file(std::vector<T> &v, BGZF *fp);
template <typename T> void read_vector_from_file(std::vector<T> &v, BGZF *fp);
template <typename T> void write_element_to_file(T &v, BGZF *fp);
template <typename T> void read_element_from_file(T &v, BGZF *fp);

std::vector<std::string> read_lines_from_file(const char *fn);

struct bgzf_index_t {
    uint64_t caddr;
    uint64_t uaddr;
};

struct GziFile {
    std::vector<bgzf_index_t> vec;
    // buffer
    std::string block_buffer;
    GziFile(const char *fn);
};

// Split a string into fields seprated by delimiters.
//
// The StringViewSplitter object should be instantilized out of loop scope to avoid
// repeated allocation. Two types of usages
//  1. split + get: for heterogenous type of data. This involves interal storage of
//  positions (string_view)
//  2. split_to_vector: this is more efficient for homogenous type of data.

class StringViewSplitter
{
    std::vector<std::string_view> sv_vec;
    std::string delim;
    std::string buffer;

  public:
    StringViewSplitter(const char *delim_, size_t reserved_size = 100);

    // @ max_fields_to_extract: 0 means all, >0 mean the number of fields try to extract.
    // @ return: void. Use get method to get each field.
    void split(std::string_view line_sv, size_t max_fields_to_extract = 0);

    template <typename T> void get(size_t which, T &val);

    // @ max_fields_to_extract: 0 means all, >0 mean the number of fields try to
    // extract.
    // @ return: void. The output is updated in `vec` parameter
    // @ internal empty field is parsed, the trailing delimiter will not parsed as an
    // empty field
    template <typename T>
    void split_to_vector(std::string_view line_sv, std::vector<T> &vec,
        bool clear_vec = true, size_t max_fields_to_extract = 0);
};

template <typename T> class TournamentTree
{
  public:
    // senitel value
    T max_val;
    size_t size_max;

  private:
    // tree structure
    struct Node {
        size_t parent;
        size_t leaf;
        T val;
    };
    std::vector<Node> node_vec;
    // size and position
    size_t N, k;
    size_t layer_N_start;

    Node &
    root()
    {
        return node_vec[0];
    }

  public:
    TournamentTree(int k_, T max_val);

    [[nodiscard]] T &init_run(std::vector<T> initial_vals, size_t &winner_id);

    // when no T value is specified, the function the flush out values filled in the
    // tournament tree
    T replace_run(size_t winner_id);

    // @ winner_id receives the next winner. (output)
    // @ T val is the input number to replace the winner's slot.
    // @ return the the winner's value.
    T &replace_run(T val, size_t &winner_id);
};

#endif
