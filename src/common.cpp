#include "common.hpp"
#include <iostream>
void
exit_on_false(bool condition, const char *message, const char *file, int lineno)
{
    if (!condition) {
        std::cerr << "Error from " << file << ":" << lineno << ": " << message << "\n";
        exit(-1);
    }
}

void
exit_with_message(const char *message)
{
    std::cerr << "Error: " << message << "\n";
    exit(-1);
}

void
add_exclusion_range(
    std::vector<region_label_t> &label_v, uint32_t pid_left, uint32_t pid_right)
{
    exit_on_false(label_v.size() > 0, "", __FILE__, __LINE__);
    exit_on_false(label_v[0].pid_s <= pid_left, "", __FILE__, __LINE__);
    exit_on_false(pid_left < pid_right, "", __FILE__, __LINE__);

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

ScopedTimer::ScopedTimer(const char *func, bool in_human_readable)
    : human_readable(in_human_readable), function_name_{ func }, start_{
          ClockType::now()
      }
{
}
ScopedTimer::~ScopedTimer()
{
    using namespace std::chrono;
    using days = duration<int, std::ratio<86400> >;

    auto stop = ClockType::now();
    auto duration = (stop - start_);
    auto ms = duration_cast<milliseconds>(duration);
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

template <typename T>
TournamentTree<T>::TournamentTree(int k_, T max_val) : max_val(max_val), k(k_)
{
    exit_on_false(
        k >= 1, "the number of node needs to be greater than 1", __FILE__, __LINE__);

    for (N = 1; (1UL << N) < (size_t) k; N++)
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
template <typename T>
[[nodiscard]] T &
TournamentTree<T>::init_run(std::vector<T> initial_vals, size_t &winner_id)
{
    exit_on_false(k == initial_vals.size(), "", __FILE__, __LINE__);
    std::vector<size_t> layer_start_vec;
    for (size_t layer = 0, count = 0; layer <= N; count += (1 << layer), layer++) {
        layer_start_vec.push_back(count);
        // std::cout << "layer: " << layer << "  start_id: " << count << '\n';
    }
    layer_N_start = layer_start_vec.back();

    for (size_t layer = N; layer != size_max; layer--) {
        // std::cout << "for layer: " << layer << '\n';
        for (size_t i = 0; i < (1UL << layer); i++) {
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

// @ winner_id receives the next winner. (output)
// @ T val is the input number to replace the winner's slot.
// @ return the the winner's value.
template <typename T>
T &
TournamentTree<T>::replace_run(T val, size_t &winner_id)
{
    size_t this_node_id = layer_N_start + root().leaf;
    // size_t parent_id, sister_id;

    // update this node's value
    node_vec[this_node_id].val = val;

    // update parents
    for (size_t layer = N; layer != 0; layer--) {
        Node &this_node = node_vec[this_node_id];
        Node &parent = node_vec[this_node.parent];
        size_t sister_id = (this_node_id & 1) ? (this_node_id + 1) : (this_node_id - 1);
        Node &sister = node_vec[sister_id];
        parent.leaf = this_node.val < sister.val ? this_node.leaf : sister.leaf;
        parent.val = this_node.val < sister.val ? this_node.val : sister.val;

        this_node_id = this_node.parent;
    }
    winner_id = root().leaf;
    return root().val;
}

template <typename T>
T
TournamentTree<T>::replace_run(size_t winner_id)
{
    return replace_run(max_val, winner_id);
}

template class TournamentTree<ibd_rec1_t>;
template class TournamentTree<uint16_t>;
