#ifndef __common_hpp__
#define __common_hpp__

#include <cstdint>
#include <cstdlib>
#include <htslib/bgzf.h>
#include <htslib/hts.h>

#include <algorithm>
#include <cassert>
#include <charconv>
#include <chrono>
#include <filesystem>
#include <fmt/core.h>
#include <fstream>
#include <iostream>
#include <istream>
#include <iterator>
#include <list>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

// from htslib bgzf.c
typedef struct {
    uint64_t uaddr; // offset w.r.t. uncompressed data
    uint64_t caddr; // offset w.r.t. compressed data
} bgzidx1_t;

// from htslib bgzf.c
struct bgzidx_t {
    int noffs, moffs;     // the size of the index, n:used, m:allocated
    bgzidx1_t *offs;      // offsets
    uint64_t ublock_addr; // offset of the current block (uncompressed data)
};

struct region_label_t {
    uint32_t pid_s : 24; // region start position id
    uint32_t label : 8;
    friend bool
    operator==(const region_label_t &x, const region_label_t &y)
    {
        return x.pid_s == y.pid_s && x.label == y.label;
    }
};

inline void
add_exclusion_range(
    std::vector<region_label_t> &label_v, uint32_t pid_left, uint32_t pid_right)
{
    assert(label_v.size() > 0 && label_v[0].pid_s <= pid_left);
    assert(pid_left < pid_right);

    using namespace std;
    size_t id1, id2;
    region_label_t tmp;
    bool found_pid_left, found_pid_right;
    uint8_t left_label, right_label;
    region_label_t sentinel = { 0xffffff, 0xff };

    auto Cmp = [](const region_label_t &x, const region_label_t &y) -> bool {
        return x.pid_s < y.pid_s;
    };

    tmp = { pid_left };

    auto it1 = upper_bound(label_v.begin(), label_v.end(), tmp, Cmp);
    found_pid_left = (it1 - 1)->pid_s == tmp.pid_s;
    left_label = (it1 - 1)->label;
    tmp.label = left_label; // ensure first insertion do not change label in front the
                            // next insertion element before the second insertion

    // add left pid if needed
    if (found_pid_left) {
        id1 = distance(label_v.begin(), it1 - 1);
    } else {
        id1 = distance(label_v.begin(), it1);
        label_v.insert(it1, tmp);
    }

    tmp = { pid_right };
    auto it2 = upper_bound(label_v.begin(), label_v.end(), tmp, Cmp);
    found_pid_right = (it2 - 1)->pid_s == tmp.pid_s;
    right_label = (it2 - 1)->label;

    // add right pid if needed
    if (found_pid_right) {
        id2 = distance(label_v.begin(), it2 - 1);
    } else {
        id2 = distance(label_v.begin(), it2);
        label_v.insert(it2, tmp);
    }

    label_v[id2].label = right_label;

    // then update labels for id1, remove element from id1 and id2;
    label_v[id1].label = 0;

    // remove those need to be removed
    size_t i, j;
    i = id1 == 0 ? 0 : id1 - 1;
    for (j = i + 1; j <= id2; j++) {
        if (j > id1 && j < id2)
            label_v[j] = sentinel;
        else if (label_v[i].label == label_v[j].label) {
            // two nearby effective element using the same label
            label_v[j] = sentinel;
        } else {
            i = j;
        }
    }

    label_v.erase(remove(label_v.begin(), label_v.end(), sentinel), label_v.end());
}

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
        assert(is_equal(*this, r2));
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

    void
    random()
    {
        int a = rand();
        const size_t sz = sizeof(int);
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
               && ((uint16_t *) &rec1)[4] == ((uint16_t *) &rec1)[4];
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

// Adapted from cpp_high_perforamnce_2e_bjorn_andrist
class ScopedTimer
{
    using ClockType = std::chrono::steady_clock;
    bool human_readable{ false };

  private:
    const char *function_name_{};
    const ClockType::time_point start_{};

  public:
    ScopedTimer(const char *func, bool human_readable = false)
        : function_name_{ func }, start_{ ClockType::now() },
          human_readable(human_readable)
    {
    }
    ScopedTimer(ScopedTimer &&) = delete;
    ScopedTimer &operator=(const ScopedTimer &) = delete;
    ScopedTimer &operator=(ScopedTimer &&) = delete;

    ~ScopedTimer()
    {
        using namespace std::chrono;
        using days = duration<int, std::ratio<86400> >;

        auto stop = ClockType::now();
        auto duration = (stop - start_);
        auto ms = duration_cast<milliseconds>(duration);
        uint8_t flags = 0;
        if (human_readable) {
            auto d = duration_cast<days>(duration);
            duration -= d;
            auto h = duration_cast<hours>(duration);
            duration -= h;
            auto m = duration_cast<minutes>(duration);
            duration -= m;
            auto s = duration_cast<seconds>(duration);
            duration -= s;
            auto ms = duration_cast<milliseconds>(duration);
            if (d.count() > 0)
                std::cerr << d.count() << "days ";
            if (d.count() || h.count())
                std::cerr << h.count() << "hrs ";
            if (d.count() || h.count() || m.count())
                std::cerr << m.count() << "min ";
            if (d.count() || h.count() || m.count() || s.count())
                std::cerr << s.count() << "sec ";
            std::cerr << ms.count() << "ms\t" << function_name_ << '\n';
        } else {
            std::cerr << ms.count() << "ms\t" << function_name_ << '\n';
        }
    }
};

template <typename T>
void
write_vector_to_file(std::vector<T> &v, BGZF *fp)
{
    if constexpr (std::is_same_v<T, std::string>) {
        size_t total_bytes = 0;
        size_t bytes_written = 0;
        std::string buffer_str;

        // calculate total size
        for (auto &s : v) {
            assert(std::find(s.begin(), s.end(), '\n') == s.end()
                   && "string should not contains '\n'");
            total_bytes += s.size() + 1; // add 1 for '\n'
        }

        // write size info
        assert(bgzf_write(fp, &total_bytes, sizeof(total_bytes)) == sizeof(total_bytes));

        // write contents via buffer
        if (total_bytes > 0) {
            buffer_str.reserve(total_bytes);
            for (auto &s : v) {
                buffer_str += s;
                buffer_str += '\n';
            }
            assert(bgzf_write(fp, buffer_str.c_str(), total_bytes) == total_bytes
                   && "write_vector_to_file error for string type");
        }

    } else if constexpr (
        std::is_arithmetic_v<
            T> || std::is_same_v<T, ibd_rec1_t> || std::is_same_v<T, ibd_rec2_t>) {
        // write size info
        size_t total_bytes = v.size() * sizeof(T);
        assert(bgzf_write(fp, &total_bytes, sizeof(total_bytes)) == sizeof total_bytes);
        // write contents
        assert(bgzf_write(fp, &v[0], total_bytes) == total_bytes);

    } else {
        assert(false && "write_vector_to_file is not implemented for this Type\n");
    }
}

template <typename T>
void
read_vector_from_file(std::vector<T> &v, BGZF *fp)
{
    if constexpr (std::is_same_v<T, std::string>) {
        size_t total_bytes;
        std::string buffer_str;

        // get size info from file
        assert(bgzf_read(fp, &total_bytes, sizeof(total_bytes)) == sizeof(total_bytes)
               && "read_vector_from_file size info reading error for string type");

        if (total_bytes > 0) {
            // get content from file
            buffer_str.resize(total_bytes);
            assert(bgzf_read(fp, &buffer_str[0], total_bytes) == total_bytes
                   && "read_vector_from_file content reading error for string type");

            // fill the vector
            std::istringstream iss(buffer_str);
            std::string line;
            v.clear();
            while (std::getline(iss, line, '\n')) {
                v.push_back(line);
            }
        }

    } else if constexpr (
        std::is_arithmetic_v<
            T> || std::is_same_v<T, ibd_rec1_t> || std::is_same_v<T, ibd_rec2_t>) {
        size_t total_bytes;
        size_t vector_sz;

        // get size info from file
        assert(bgzf_read(fp, &total_bytes, sizeof(total_bytes)) == sizeof(total_bytes)
               && "read_vector_from_file size info reading error for "
                  "arithemtic/ibd_rec types");

        // get content from file
        vector_sz = total_bytes / sizeof(T);
        v.resize(vector_sz);
        assert(bgzf_read(fp, &v[0], total_bytes) == total_bytes
               && "read_vector_from_file content reading error for arithemtic/ibd_rec "
                  "types");

    } else {
        assert(false && "read_vector_to_file is not implemented for this Type\n");
    }
}

template <typename T>
void
write_element_to_file(T &v, BGZF *fp)
{
    if constexpr (
        std::is_arithmetic_v<
            T> || std::is_same_v<T, ibd_rec1_t> || std::is_same_v<T, ibd_rec2_t>) {
        // get content from file
        assert(
            bgzf_write(fp, &v, sizeof v) == sizeof v && "write_element_from_file error");

    } else {
        assert(false && "read_element_to_file is not implemented for this Type\n");
    }
}
template <typename T>
void
read_element_from_file(T &v, BGZF *fp)
{
    if constexpr (
        std::is_arithmetic_v<
            T> || std::is_same_v<T, ibd_rec1_t> || std::is_same_v<T, ibd_rec2_t>) {
        // get content from file
        assert(
            bgzf_read(fp, &v, sizeof v) == sizeof v && "read_element_from_file error");

    } else {
        assert(false && "read_element_to_file is not implemented for this Type\n");
    }
}

inline std::vector<std::string>
read_lines_from_file(const char *fn)
{
    std::vector<std::string> lines;
    BGZF *fp = bgzf_open(fn, "r");
    assert(fp != NULL && "bgzf_open error");
    kstring_t kstr = { 0 };
    while (bgzf_getline(fp, '\n', &kstr) >= 0) {
        lines.push_back(std::string(ks_c_str(&kstr)));
    }
    ks_free(&kstr);
    bgzf_close(fp);
    return lines;
}

struct bgzf_index_t {
    uint64_t caddr;
    uint64_t uaddr;
};

struct GziFile {
    std::vector<bgzf_index_t> vec;
    // buffer
    std::string block_buffer;

    GziFile(const char *fn)
    {
        vec.clear();
        FILE *fp = fopen(fn, "r");
        uint64_t sz;
        assert(fread(&sz, sizeof(sz), 1, fp) == 1);
        vec.resize(sz);
        assert(fread(&vec[0], sizeof(decltype(vec)::value_type), sz, fp) == sz);
        fclose(fp);

        block_buffer.resize(64 * 1024);
    }
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
    StringViewSplitter(const char *delim_, size_t reserved_size = 100) : delim{ delim_ }
    {
        if (reserved_size > 0)
            sv_vec.reserve(reserved_size);
    }

    // @ max_fields_to_extract: 0 means all, >0 mean the number of fields try to extract.
    // @ return: void. Use get method to get each field.
    void
    split(std::string_view line_sv, size_t max_fields_to_extract = 0)
    {

        sv_vec.clear();
        size_t s, e;
        size_t ncol = 0;
        size_t back = line_sv.size() - 1;
        for (s = 0, e = 0; ncol < max_fields_to_extract || max_fields_to_extract == 0;
             ncol++) {
            // Not using following because it skips empty field.
            //     s = line_sv.find_first_not_of(delim, e);
            //     if (s == line_sv.npos)
            //        break;
            //
            if (e >= back)
                break;
            if (e != 0)
                s = e + 1;
            e = line_sv.find_first_of(delim, s);
            sv_vec.push_back(line_sv.substr(s, e - s));

            /*
            std::cerr << "s: " << s << " e: " << e << " sv: " << line_sv.substr(s, e - s)
                      << " e - s: " << e - s
                      << " sv.size: " << line_sv.substr(s, e - s).size()
                      << " sv_vec.size: " << sv_vec.size() << '\n';
                      */
        }
    }

    template <typename T>
    void
    get(size_t which, T &val)
    {
        assert(which < sv_vec.size() && "out of range");
        if constexpr (std::is_same_v<T, std::string>) {
            val.clear();
            val += sv_vec[which];
            /*
            std::cout << "type: string\t";
            */
        } else if constexpr (std::is_same_v<T, std::string_view>) {
            val = sv_vec[which];
            /*
            std::cout << "type: string_view\t";
            */
        } else if constexpr (std::is_same_v<T, double>) {
            // Note important to set a default value, because from_chars will not
            // update the val if the first=back (length=0)
            auto sv = sv_vec[which];
            val = std::strtod(sv.begin(), NULL);
            /*
            std::cout << "type: arithmetic\t";
            */
        } else if constexpr (std::is_arithmetic_v<T>) {
            // Note important to set a default value, because from_chars will not
            // update the val if the first=back (length=0)
            val = 0;
            auto sv = sv_vec[which];
            std::from_chars(sv.begin(), sv.end(), val);
            /*
            std::cout << "type: arithmetic\t";
            */
        } else {
            assert(false && "get method Not implimented for this type");
        }
        /*
        std::cout << " field: " << which << " value: " << val << '\n';
        */
    }

    // @ max_fields_to_extract: 0 means all, >0 mean the number of fields try to
    // extract.
    // @ return: void. The output is updated in `vec` parameter
    // @ internal empty field is parsed, the trailing delimiter will not parsed as an
    // empty field
    template <typename T>
    void
    split_to_vector(std::string_view line_sv, std::vector<T> &vec, bool clear_vec = true,
        size_t max_fields_to_extract = 0)
    {
        // clear the output vector for each, otherwise just append.
        if (clear_vec)
            vec.clear();

        size_t s, e;
        size_t back = line_sv.size() - 1;
        size_t ncol = 0;
        T val;
        std::string_view sv;

        for (s = 0, e = 0; ncol < max_fields_to_extract || max_fields_to_extract == 0;
             ncol++) {
            if (e >= back)
                break;
            if (e != 0)
                s = e + 1;
            e = line_sv.find_first_of(delim, s);
            sv = line_sv.substr(s, e - s);

            /*
            std::cerr << "s: " << s << " e: " << e << " sv: " << sv
                      << " e - s: " << e - s << " sv.size: " << sv
                      << " sv_vec.size: " << sv_vec.size() << '\n';
                      */

            if constexpr (std::is_same_v<T, std::string>) {
                buffer.clear();
                buffer += sv;
                /*
                std::cout << "type: string\t";
                */
                vec.push_back(buffer);
            } else if constexpr (std::is_same_v<T, std::string_view>) {
                vec.push_back(sv);
                /*
                std::cout << "type: string_view\t";
                */
            } else if constexpr (std::is_arithmetic_v<T>) {
                // Note important to set a default value, because from_chars will not
                // update the val if the first=back (length=0)
                val = 0;
                std::from_chars(sv.begin(), sv.end(), val);
                vec.push_back(val);
                /*
                std::cout << "type: arithmetic\t";
                */
            } else {
                assert(false && "get method Not implimented for this type");
            }
            /*
            std::cout << " field: " << vec.size() - 1 << " value: " << vec.back()
                      << '\n';
                      */
        }
    }
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
    TournamentTree(int k_, T max_val) : max_val(max_val), k(k_)
    {
        assert(k >= 1 && "the number of node needs to be greater than 1");

        for (N = 1; (1 << N) < k; N++)
            ;
        // std::cout << "k = " << k << " N: " << N << " 2^n = " << (1 << N) << '\n';
        size_t sz = 0;
        for (size_t layer = 0; layer <= N; layer++)
            sz += (1 << layer);
        // std::cout << "size = " << sz << '\n';
        node_vec.resize(sz);

        size_max = std::numeric_limits<size_t>::max();
    }

    /*
     *                       0
     *           1                       2
     *     3           4           5           6
     *  7     8     9     10    11    12    13    14
     */
    // @ winner_id receives the next winner. (output)
    // @ initial_val: is a vector values to initial the K slots.
    // @ return the the winner's value.
    [[nodiscard]] T &
    init_run(std::vector<T> initial_vals, size_t &winner_id)
    {
        assert(k == initial_vals.size());
        std::vector<size_t> layer_start_vec;
        for (size_t layer = 0, count = 0; layer <= N; count += (1 << layer), layer++) {
            layer_start_vec.push_back(count);
            // std::cout << "layer: " << layer << "  start_id: " << count << '\n';
        }
        layer_N_start = layer_start_vec.back();

        for (size_t layer = N; layer != size_max; layer--) {
            // std::cout << "for layer: " << layer << '\n';
            for (size_t i = 0; i < (1 << layer); i++) {
                size_t this_node_id = layer_start_vec[layer] + i;
                Node &node = node_vec[this_node_id];
                // initialize parent
                if (layer > 0) {
                    size_t prev_layer_start = layer_start_vec[layer - 1];
                    node.parent = prev_layer_start + i / 2;
                } else {
                    node.parent = size_max;
                }
                // initialize leaf and initialize value_id
                if (layer == N) {
                    node.leaf = i;
                    if (i < initial_vals.size())
                        node.val = initial_vals[i];
                    else {
                        node.val = max_val;
                    }
                } else {
                    size_t next_layer_start = layer_start_vec[layer + 1];
                    size_t left_child_node_id = next_layer_start + 2 * i;
                    size_t right_child_node_id = left_child_node_id + 1;
                    Node &left = node_vec[left_child_node_id];
                    Node &right = node_vec[right_child_node_id];
                    if (left.val < right.val) {
                        node.leaf = left.leaf;
                        node.val = left.val;
                    } else {
                        node.leaf = right.leaf;
                        node.val = right.val;
                    }
                }

                // std::cout << "node: " << this_node_id << "\t i: " << i
                //           << "\t leaf: " << node.leaf << "\t val: " << node.val
                //           << "\t parent: " << node.parent << '\n';
            }
        }
        winner_id = root().leaf;
        return root().val;
    }

    // when no T value is specified, the function the flush out values filled in the
    // tournament tree
    size_t
    replace_run(size_t winner_id)
    {
        return replace_run(max_val, winner_id);
    }

    // @ winner_id receives the next winner. (output)
    // @ T val is the input number to replace the winner's slot.
    // @ return the the winner's value.
    T &
    replace_run(T val, size_t &winner_id)
    {
        size_t this_node_id = layer_N_start + root().leaf;
        size_t parent_id, sister_id;

        // update this node's value
        node_vec[this_node_id].val = val;

        // update parents
        for (size_t layer = N; layer != 0; layer--) {
            Node &this_node = node_vec[this_node_id];
            Node &parent = node_vec[this_node.parent];
            size_t sister_id
                = (this_node_id & 1) ? (this_node_id + 1) : (this_node_id - 1);
            Node &sister = node_vec[sister_id];
            parent.leaf = this_node.val < sister.val ? this_node.leaf : sister.leaf;
            parent.val = this_node.val < sister.val ? this_node.val : sister.val;

            // std::cout << "this_id: " << this_node_id << " this_leaf: " <<
            // this_node.leaf
            //           << "  sister_id: " << sister_id << "  sister_leaf: " <<
            //           sister.leaf
            //           << " parent_id: " << this_node.parent
            //           << " parent_val: " << parent.val << " parent.leaf: " <<
            //           parent.leaf
            //           << '\n';

            this_node_id = this_node.parent;
        }
        winner_id = root().leaf;
        return root().val;
    }
};

#endif
