
#include "common.hpp"
#include <htslib/bgzf.h>
#include <sstream>

template void write_vector_to_file(std::vector<std::string> &v, BGZF *fp);
template void write_vector_to_file(std::vector<int> &v, BGZF *fp);
template void write_vector_to_file(std::vector<uint32_t> &v, BGZF *fp);
template void write_vector_to_file(std::vector<float> &v, BGZF *fp);
template void write_vector_to_file(std::vector<size_t> &v, BGZF *fp);
template void write_vector_to_file(std::vector<long double> &v, BGZF *fp);
template void write_vector_to_file(std::vector<uint8_t> &v, BGZF *fp);
template void write_vector_to_file(std::vector<uint16_t> &v, BGZF *fp);

template <typename T>
void
write_vector_to_file(std::vector<T> &v, BGZF *fp)
{
    if constexpr (std::is_same_v<T, std::string>) {
        uint64_t total_bytes = 0;
        bool condition = false;
        std::string buffer_str;
        ssize_t ret = 0;

        // calculate total size
        for (auto &s : v) {
            exit_on_false(std::find(s.begin(), s.end(), '\n') == s.end(),
                "string should not contains '\n'", __FILE__, __LINE__);
            total_bytes += s.size() + 1; // add 1 for '\n'
        }

        // write size info
        ret = bgzf_write(fp, &total_bytes, sizeof(total_bytes));
        condition = ret >= 0 && ((size_t) ret) == sizeof(total_bytes);
        exit_on_false(condition, "", __FILE__, __LINE__);

        // write contents via buffer
        if (total_bytes > 0) {
            buffer_str.reserve(total_bytes);
            for (auto &s : v) {
                buffer_str += s;
                buffer_str += '\n';
            }
            ret = bgzf_write(fp, buffer_str.c_str(), total_bytes);
            condition = ret >= 0 && ((size_t) ret) == total_bytes;
            const char *msg = "write_vector_to_file error for string type";
            exit_on_false(condition, msg, __FILE__, __LINE__);
        }

    } else if constexpr (
        std::is_arithmetic_v<
            T> || std::is_same_v<T, ibd_rec1_t> || std::is_same_v<T, ibd_rec2_t>) {
        // write size info
        ssize_t ret = 0;
        uint64_t total_bytes = v.size() * sizeof(T);
        bool condition = false;

        ret = bgzf_write(fp, &total_bytes, sizeof(total_bytes));
        condition = ret >= 0 && ((size_t) ret) == sizeof(total_bytes);
        exit_on_false(condition, "", __FILE__, __LINE__);

        // write contents
        ret = bgzf_write(fp, &v[0], total_bytes);
        condition = ret >= 0 && ((size_t) ret) == total_bytes;
        const char *msg = "write_vector_to_file error for string type";
        exit_on_false(condition, msg, __FILE__, __LINE__);

    } else {
        std::cerr << "write_vector_to_file is not implemented for this Type\n";
        exit(-1);
    }
}

template void read_vector_from_file(std::vector<std::string> &v, BGZF *fp);
template void read_vector_from_file(std::vector<int> &v, BGZF *fp);
template void read_vector_from_file(std::vector<uint32_t> &v, BGZF *fp);
template void read_vector_from_file(std::vector<float> &v, BGZF *fp);
template void read_vector_from_file(std::vector<size_t> &v, BGZF *fp);
template void read_vector_from_file(std::vector<long double> &v, BGZF *fp);
template void read_vector_from_file(std::vector<uint8_t> &v, BGZF *fp);
template void read_vector_from_file(std::vector<uint16_t> &v, BGZF *fp);

template <typename T>
void
read_vector_from_file(std::vector<T> &v, BGZF *fp)
{
    if constexpr (std::is_same_v<T, std::string>) {
        uint64_t total_bytes = 0;
        ssize_t ret = 0;
        std::string buffer_str;

        // get size info from file
        ret = bgzf_read(fp, &total_bytes, sizeof(total_bytes));
        exit_on_false(
            ret >= 0 && ((size_t) ret) == sizeof(total_bytes), "", __FILE__, __LINE__);

        if (total_bytes > 0) {
            // get content from file
            buffer_str.resize(total_bytes);
            ret = bgzf_read(fp, &buffer_str[0], total_bytes);
            exit_on_false(ret >= 0 && ((size_t) ret) == total_bytes,
                "read_vector_from_file content reading error for string type", __FILE__,
                __LINE__);

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
        uint64_t total_bytes = 0;
        size_t vector_sz = 0;
        ssize_t ret = 0;

        // get size info from file
        ret = bgzf_read(fp, &total_bytes, sizeof(total_bytes));
        exit_on_false(ret >= 0 && ((size_t) ret) == sizeof(total_bytes),
            "read_vector_from_file size info reading error for arithemtic/ibd_rec types",
            __FILE__, __LINE__);

        // get content from file
        vector_sz = total_bytes / sizeof(T);
        v.resize(vector_sz);
        ret = bgzf_read(fp, &v[0], total_bytes);
        exit_on_false(ret >= 0 && ((size_t) ret) == total_bytes,
            "read_vector_from_file content reading error for arithemtic/ibd_rec types",
            __FILE__, __LINE__);

    } else {
        exit_with_message("read_vector_to_file is not implemented for this Type\n");
    }
}

template void write_element_to_file(int &v, BGZF *fp);
template void write_element_to_file(size_t &v, BGZF *fp);
template void write_element_to_file(uint8_t &v, BGZF *fp);
template void write_element_to_file(uint32_t &v, BGZF *fp);
template void write_element_to_file(float &v, BGZF *fp);

template <typename T>
void
write_element_to_file(T &v, BGZF *fp)
{
    if constexpr (
        ssize_t ret = 0;
        std::is_arithmetic_v<
            T> || std::is_same_v<T, ibd_rec1_t> || std::is_same_v<T, ibd_rec2_t>) {
        // get content from file
        ret = bgzf_write(fp, &v, sizeof v);
        exit_on_false(ret > 0 && ((size_t) ret) == sizeof v,
            "write_element_from_file error", __FILE__, __LINE__);

    } else {
        exit_with_message("read_element_to_file is not implemented for this Type\n");
    }
}

template void read_element_from_file(int &v, BGZF *fp);
template void read_element_from_file(size_t &v, BGZF *fp);
template void read_element_from_file(uint8_t &v, BGZF *fp);
template void read_element_from_file(uint32_t &v, BGZF *fp);
template void read_element_from_file(float &v, BGZF *fp);

template <typename T>
void
read_element_from_file(T &v, BGZF *fp)
{
    if constexpr (
        ssize_t ret = 0;
        std::is_arithmetic_v<
            T> || std::is_same_v<T, ibd_rec1_t> || std::is_same_v<T, ibd_rec2_t>) {
        // get content from file
        ret = bgzf_read(fp, &v, sizeof v);
        exit_on_false(ret > 0 && ((size_t) ret) == sizeof v,
            "read_element_from_file error", __FILE__, __LINE__);

    } else {
        exit_with_message("read_element_to_file is not implemented for this Type\n");
    }
}

std::vector<std::string>
read_lines_from_file(const char *fn)
{
    std::vector<std::string> lines;
    BGZF *fp = bgzf_open(fn, "r");
    exit_on_false(fp != NULL, "bgzf_open error", __FILE__, __LINE__);
    kstring_t kstr = { 0 };
    while (bgzf_getline(fp, '\n', &kstr) >= 0) {
        lines.push_back(std::string(ks_c_str(&kstr)));
    }
    ks_free(&kstr);
    bgzf_close(fp);
    return lines;
}

GziFile::GziFile(const char *fn)
{
    vec.clear();
    ssize_t ret = 0;
    FILE *fp = fopen(fn, "r");
    uint64_t sz = 0;
    ret = fread(&sz, sizeof(sz), 1, fp);
    exit_on_false(ret >= 0 && ((size_t) ret) == 1, "", __FILE__, __LINE__);
    vec.resize(sz);
    ret = fread(&vec[0], sizeof(decltype(vec)::value_type), sz, fp);
    exit_on_false(ret >= 0 && ((size_t) ret) == sz, "", __FILE__, __LINE__);
    fclose(fp);

    block_buffer.resize(64 * 1024);
}
