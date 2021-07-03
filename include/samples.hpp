#ifndef __samples_hpp__
#define __samples_hpp__

#include <cassert>
#include <unordered_map>

#include "common.hpp"

class Samples
{
    std::vector<std::string> names_vec;
    std::unordered_map<std::string, uint32_t> name2id_map;

  public:
    void
    add(std::string name)
    {
        if (name2id_map.find(name) == name2id_map.end()) {
            size_t next_id = name2id_map.size();
            name2id_map[name] = next_id;
            names_vec.push_back(name);
        }
    }

    std::string &
    get_name(uint32_t sid)
    {
        return names_vec[sid];
    }

    uint32_t
    get_id(const std::string name)
    {
        return name2id_map.at(name);
    }

    // for debugging
    void
    print()
    {
        for (auto &str : names_vec) {
            std::cout << str << '\t';
        }
        std::cout << '\n';
    }

    void
    write_to_file(BGZF *fp)
    {
        write_vector_to_file(names_vec, fp);
    }

    void
    read_from_file(BGZF *fp)
    {
        read_vector_from_file(names_vec, fp);
        name2id_map.clear();
        for (auto &n : names_vec) {
            size_t next_id = name2id_map.size();
            name2id_map[n] = next_id;
        }
    }

    bool
    is_equal(Samples &other)
    {
        return names_vec == other.names_vec && name2id_map == other.name2id_map;
    }

    uint32_t
    get_num_samples()
    {
        return names_vec.size();
    }
};

#endif
