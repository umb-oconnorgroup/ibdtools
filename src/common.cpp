#include "common.hpp"
#include <iostream>

void
exit_with_message(const char *message, const char *file, const int lineno)
{
    std::cerr << "Error: " << (message) << " (from " << (file) << " + " << (lineno)
              << ") \n";
    exit(-1);
}

void
add_exclusion_range(
    std::vector<region_label_t> &label_v, uint32_t pid_left, uint32_t pid_right)
{
    my_assert(label_v.size() > 0, "");
    my_assert(label_v[0].pid_s <= pid_left, "");
    my_assert(pid_left < pid_right, "");

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
