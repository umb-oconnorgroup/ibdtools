#include "common.hpp"
#include <charconv>

StringViewSplitter::StringViewSplitter(const char *delim_, size_t reserved_size)
    : delim{ delim_ }
{
    if (reserved_size > 0)
        sv_vec.reserve(reserved_size);
}

void
StringViewSplitter::split(std::string_view line_sv, size_t max_fields_to_extract)
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

template void StringViewSplitter::get(size_t which, int &val);
template void StringViewSplitter::get(size_t which, double &val);
template void StringViewSplitter::get(size_t which, uint32_t &val);
template void StringViewSplitter::get(size_t which, std::string &val);
template void StringViewSplitter::get(size_t which, std::string_view &val);

template <typename T>
void
StringViewSplitter::get(size_t which, T &val)
{
    exit_on_false(which < sv_vec.size(), "out of range", __FILE__, __LINE__);
    if constexpr (std::is_same_v<T, std::string>) {
        val.clear();
        val += sv_vec[which];
    } else if constexpr (std::is_same_v<T, std::string_view>) {
        val = sv_vec[which];
    } else if constexpr (std::is_same_v<T, double>) {
        // Note important to set a default value, because from_chars will not
        // update the val if the first=back (length=0)
        auto sv = sv_vec[which];
        val = std::strtod(sv.begin(), NULL);
    } else if constexpr (std::is_arithmetic_v<T>) {
        // Note important to set a default value, because from_chars will not
        // update the val if the first=back (length=0)
        val = 0;
        auto sv = sv_vec[which];
        std::from_chars(sv.begin(), sv.end(), val);
    } else {
        exit_with_message("get method Not implimented for this type");
    }
}

template void StringViewSplitter::split_to_vector(
    std::string_view, std::vector<int> &, bool, size_t);

template <typename T>
void
StringViewSplitter::split_to_vector(std::string_view line_sv, std::vector<T> &vec,
    bool clear_vec, size_t max_fields_to_extract)
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

        if constexpr (std::is_same_v<T, std::string>) {
            buffer.clear();
            buffer += sv;
            vec.push_back(buffer);
        } else if constexpr (std::is_same_v<T, std::string_view>) {
            vec.push_back(sv);
        } else if constexpr (std::is_arithmetic_v<T>) {
            // Note important to set a default value, because from_chars will not
            // update the val if the first=back (length=0)
            val = 0;
            std::from_chars(sv.begin(), sv.end(), val);
            vec.push_back(val);
        } else {
            exit_with_message("get method Not implimented for this type");
        }
    }
}
