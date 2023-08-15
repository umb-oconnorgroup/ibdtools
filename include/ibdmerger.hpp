#ifndef __ibdmerger_hpp__
#define __ibdmerger_hpp__

#include "./common.hpp"
#include "./ibdfile.hpp"
#include "./metafile.hpp"
#include "genotypes.hpp"
#include "gmap.hpp"
#include "positions.hpp"

class IbdMerger
{
    IbdFile in;
    IbdFile out;

    MetaFile meta;

    size_t max_rec_allowed_by_ram;
    std::string out_mode;

    struct range_t {
        size_t first, last;
    };

    std::vector<range_t> range_vec;

    // threshold
    float max_cm;
    size_t max_snp;

    std::vector<ibd_rec1_t> &vec{ in.get_vec() };

  public:
    IbdMerger(const char *in_fn, const char *out_fn, const char *out_mode,
        const char *meta_fn, size_t max_rec_allowed_by_ram, size_t max_snp = 1,
        float max_cm = 0.6)
        : max_rec_allowed_by_ram(max_rec_allowed_by_ram), out_mode(out_mode),
          max_cm(max_cm), max_snp(max_snp)
    {
        // MetaFile object can not be a local variable in this function. make it a
        // member variable
        //
        // read meta
        BGZF *fp = bgzf_open(meta_fn, "r");
        exit_on_false(fp != NULL, "", __FILE__, __LINE__);
        meta.read_from_file(fp);
        bgzf_close(fp);

        in = IbdFile(in_fn, &meta, this->max_rec_allowed_by_ram);
        out = IbdFile(out_fn, NULL, 0);
    }

    void
    merge()
    {
        in.open("r");
        out.open(out_mode.c_str());
        auto &vec = in.get_vec();

        bool read_full;
        do {
            // read and append to in's vector
            read_full = in.read_from_file(true);

            size_t last = 0;
            if (vec.size() > 0) {
                last = find_group_start(vec.size() - 1);
            }

            // This is can split to many regions for paralellization
            merge_range(0, last);

            // squish
            auto it = std::remove_if(vec.begin(), vec.begin() + last,
                [](auto x) -> bool { return x.sid1 == 0; });

            // write
            in.write_to_file(out, distance(vec.begin(), it));

            // clear the vector
            vec.erase(vec.begin(), vec.begin() + last);
        } while (read_full);

        // merge left over records and output
        merge_range(0, vec.size());
        auto it = std::remove_if(
            vec.begin(), vec.end(), [](auto x) -> bool { return x.sid1 == 0; });
        in.write_to_file(out, distance(vec.begin(), it));

        in.close();
        out.close();
    }

  private:
    size_t
    find_group_start(size_t id)
    {
        // size_t id_o = id;

        auto &vec = in.get_vec();
        exit_on_false(id >= 0 && id < vec.size(), "", __FILE__, __LINE__);

        // if the same as the one above, move above
        while (id != 0 && vec[id - 1].is_same_pair(vec[id]))
            id--;
        exit_on_false(id >= 0 && id < vec.size(), "", __FILE__, __LINE__);

        // std::cout << "id in: " << id_o;
        // std::cout << " id out: " << id << '\n';
        // if(id_o != id && id!=0)
        // {
        // 	in.get_vec()[id-1].print();
        // 	in.get_vec()[id].print();
        // }

        return id;
    }

    bool
    is_mergeable(size_t id1, size_t id2)
    {
        ibd_rec1_t &r1 = in.get_vec()[id1];
        ibd_rec1_t &r2 = in.get_vec()[id2];

        // not the same pair
        if (!(r1.is_same_pair(r2)))
            return false;

        // overlapping
        if (r1.get_pid2() >= r2.get_pid1()) {

            // std::cout << "===\n";
            // r1.print();
            // r2.print();
            return true;
        }

        //
        // non overlapping
        //
        Positions &pos = in.get_meta()->get_positions();
        Genotypes &gt = in.get_meta()->get_genotypes();
        uint32_t sid1 = r1.get_sid1();
        uint32_t sid2 = r1.get_sid2();
        uint32_t prev_pid2 = r1.get_pid2();
        uint32_t next_pid1 = r2.get_pid1();

        // too far
        if (pos.get_cm(next_pid1) - pos.get_cm(prev_pid2) > max_cm)
            return false;

        // too many discordant homozygotes
        size_t num_discordant_homezygotes = 0;
        uint8_t a, b, c, d;
        for (uint32_t id = prev_pid2 + 1;
             id < next_pid1 && num_discordant_homezygotes <= max_snp; id++) {
            //    std::cout << " id " << id << " sid1 " << sid1 << " sid2 " << sid2 <<
            //    '\n'
            //	  ;

            a = gt.at(id, sid1);
            b = a >> 2;
            a &= 0b11;
            c = gt.at(id, sid2);
            d = c >> 2;
            c &= 0b11;

            exit_on_false(a <= 1 && b <= 1 && c <= 1 && d <= 1, "", __FILE__, __LINE__);

            // sample1 is homozgyote and sample 2 is homozygote but their genotype is
            // different, then add 1;
            // if (a == b && a <= 1 && c == d && c <= 1 && a != c)
            if (a == b && c == d && a != c)
                num_discordant_homezygotes++;
        }

        if (num_discordant_homezygotes > max_snp)
            return false;
        return true;
    }

    void
    merge_range(size_t first, size_t last)
    {
        auto &vec = in.get_vec();

        // id1 points to the last survived record; id2 points to the rec to be checked
        size_t id1 = first, id2 = first + 1;

        while (id2 < last) {
            if (is_mergeable(id1, id2)) {
                // merge but extend the end of IBD segment
                vec[id1].set_pid2(std::max(vec[id1].get_pid2(), vec[id2].get_pid2()));
                // set 0 to mark the record that will not be present in out
                vec[id2].sid1 = 0;
                id2++;
            } else {
                // before move to the next, set the hap bits to zeros
                vec[id1].set_hid1(0);
                vec[id1].set_hid2(0);

                // std::cout << "keep : ";
                // vec[id1].print();

                id1 = id2;
                id2++;
            }
        }
        // the last rec needs also to reset hap bits
        vec[id1].set_hid1(0);
        vec[id1].set_hid2(0);

        // std::cout << "keep : ";
        // vec[id1].print();
    }
};

#endif
