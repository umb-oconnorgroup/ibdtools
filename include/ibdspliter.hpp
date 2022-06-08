#ifndef __ibdsplitter_hpp__
#define __ibdsplitter_hpp__

#include "./common.hpp"
#include "./ibdfile.hpp"
#include "metafile.hpp"
#include <algorithm>

class IbdSplitter
{
    std::vector<region_label_t> &labels;
    IbdFile in;
    std::vector<IbdFile> out_files;
    MetaFile meta;

    float cm_threshold;

  public:
    IbdSplitter(const char *in_fn, const char *out_fn_prefix, const char *meta_fn,
        std::vector<region_label_t> &labels_, uint32_t max_label = 1,
        float cm_threshold = 2.0, size_t in_ibdrec_vec_capacity = 1 * 1024 * 1024,
        size_t out_ibdrec_vec_capacity = 1 * 1024 * 1024)
        : labels(labels_), cm_threshold(cm_threshold)
    {
        BGZF *fp = bgzf_open(meta_fn, "r");
        verify(fp != NULL);
        meta.read_from_file(fp);
        bgzf_close(fp);

        std::string prefix = out_fn_prefix;
        in = IbdFile(in_fn, NULL, in_ibdrec_vec_capacity);

        // push an IbdFile for each label ( 0, 1, 2, ..., max_label)
        for (uint32_t i = 0; i <= max_label; i++) {
            std::string out_fn = out_fn_prefix + std::to_string(i);
            out_files.push_back(IbdFile(out_fn.c_str(), NULL, out_ibdrec_vec_capacity));
        }
    }

    void
    print_labels()
    {
        for (auto label : labels) {
            std::cout << "pid_s: " << label.pid_s << "\t label: " << label.label << '\n';
        }
    }

    // driving code
    void
    split()
    {
        in.open("r");
        for (auto &f : out_files)
            f.open("w");

        auto &in_vec = in.get_vec();

        bool read_full;
        do {
            read_full = in.read_from_file(false); // flush in info from input file

            for (auto rec : in_vec) {

                split_ibd_rec(rec); // split each record
            }

        } while (read_full); // if read_full is false, the input file reachs eof.

        in.close();

        // flush out ibd_vec to file
        for (auto &f : out_files) {
            if (f.get_vec().size() > 0) {
                f.write_to_file();
            }
            f.close();
        }
    }

    void
    split_ibd_rec(ibd_rec1_t rec)
    {
        using namespace std;
        uint32_t rid_s, rid_e;
        uint32_t pid_l, pid_r, rid;

        // ----------------
        // get a range of region id
        // ----------------
        //
        auto cmp_start_id = [](region_label_t x, region_label_t y) -> bool {
            return x.pid_s < y.pid_s;
        };
        auto lbeg = labels.begin();
        auto lend = labels.end();
        auto qleft = region_label_t{ rec.get_pid1(), 0 };
        auto qright = region_label_t{ rec.get_pid2(), 0 };

        rid_s = distance(lbeg, upper_bound(lbeg, lend, qleft, cmp_start_id));
        rid_e = distance(lbeg, upper_bound(lbeg + rid_s, lend, qright, cmp_start_id));

        ibd_rec1_t tmp_rec;
        // ----------------
        // Split by regions
        // ----------------

        /*                        rid_e
                                  |
                                  rid_s
                                  |

            ------------++++++++++----------++++++++++-----------++++++++++
                          ^----^

                        /       \
            rec.pid1        rec.pid2
        */
        // if rec is within a single split region
        if (rid_s == rid_e) {
            rid = rid_s - 1;
            pid_l = rec.get_pid1();
            pid_r = rec.get_pid2();

            if (is_longer_than_2cm(pid_l, pid_r))
                transfer(rec, rid);

        } else {
            /*                        s                                        e
                ------------++++++++++----------++++++++++-----------++++++++++
                                 ^-----------------------------------------^

                                 |    |
                                 l    r
            */
            // for the left side, need to change pid2
            rid = rid_s - 1;
            pid_l = rec.get_pid1();
            pid_r = labels[rid + 1].pid_s;

            tmp_rec = rec;
            tmp_rec.set_pid2(pid_r);
            if (is_longer_than_2cm(pid_l, pid_r))
                transfer(tmp_rec, rid);
            /*                        s                                        e
                ------------++++++++++----------++++++++++-----------++++++++++
                                 ^-----------------------------------------^

                                      |         |
                                      l         r
            */
            // in the middle, need to change both pid1 and pid2
            for (rid = rid_s; rid < rid_e - 1; rid++) {
                pid_l = labels[rid].pid_s;
                pid_r = labels[rid + 1].pid_s;

                tmp_rec = rec;
                tmp_rec.set_pid1(pid_l);
                tmp_rec.set_pid2(pid_r);
                if (is_longer_than_2cm(pid_l, pid_r))
                    transfer(tmp_rec, rid);
            }
            /*  ------------++++++++++----------++++++++++-----------++++++++++
                                 ^-----------------------------------------^

                                                                     |     |
                                                                     l     r
            */
            // for the righ side, only need to change pid1
            rid = rid_e - 1;
            pid_l = labels[rid].pid_s;
            pid_r = rec.get_pid2();

            tmp_rec = rec;
            tmp_rec.set_pid1(pid_l);
            if (is_longer_than_2cm(pid_l, pid_r))
                transfer(tmp_rec, rid);
        }
    }

    void
    transfer(ibd_rec1_t rec, uint32_t rid)
    {
        uint32_t label = labels[rid].label;

        auto &f = out_files[label];
        auto &out_vec = f.get_vec();

        out_vec.push_back(rec);

        if (out_vec.capacity() <= out_vec.size()) {
            f.write_to_file(); // flush out info to output file
            out_vec.clear();
        }
    }

  private:
    bool
    is_longer_than_2cm(uint32_t pid_l, uint32_t pid_r)
    {
        auto &pos = meta.get_positions();
        return pos.get_cm(pid_r) - pos.get_cm(pid_l) >= cm_threshold;
    }
};

#endif
