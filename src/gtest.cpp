#include "../include/ibdcoverage.hpp"
#include "../include/ibdfile.hpp"
#include "../include/ibdmatrix.hpp"
#include "../include/ibdmerger.hpp"
#include "../include/ibdsorter.hpp"
#include "../include/ibdspliter.hpp"
#include "../include/metafile.hpp"
#include <algorithm>
#include <cstdint>
#include <gtest/gtest.h>
#include <iterator>
#include <memory>
#include <numeric>
#include <random>
#include <stdio.h>

const char *map_fn = "data/example.map";
const char *vcf_fn = "data/example.bcf.gz";
const char *haps_fn = "data/example_haplotypes.txt.gz";
const char *ibd_txt_fn = "data/example.ibd.gz";
const char *browning_merged = "data/merged_browning.ibd.gz";
const char *temp_file1 = "tmp1.gz";
const char *temp_file2 = "tmp2.gz";
const char *temp_file3 = "tmp3.gz";
const char *temp_file4 = "tmp4.gz";
const char *temp_file5 = "tmp5.gz";

// 1. ksplit seems to alloc new ids buffer everytime. Bad;
// 2. ksplit_core seems to realloc ids if necessary. Good. The document was bad. It
// returns the number of field (n_fields) found. The buffer should be finally freed by
// user.
// 3. kstring_t is reused in getline. Similar to ids in ksplit_core
TEST(htslib, k_split_k_split_core)
{
    auto get_lines = []() -> std::vector<std::string> {
        std::string short_line = "1\t2\t3\t4\t";
        std::string long_line = short_line;
        long_line += long_line;
        long_line += long_line;
        long_line += long_line;
        long_line += long_line;

        std::vector<std::string> lines;
        lines.push_back(short_line);
        lines.push_back(short_line);
        lines.push_back(long_line);
        lines.push_back(short_line);
        return lines;
    };

    auto lines = get_lines();
    int *buf;
    int buf_sz;

    auto does_ksplit_cause_reallocation = [&](std::string line) mutable {
        kstring_t str;
        str.s = &line[0];
        str.l = line.size();
        str.m = line.capacity();

        int *old_buf = buf;
        buf = ksplit(&str, '\t', &buf_sz);
        return buf != old_buf;
    };

    auto does_ksplit_core_cause_reallocation = [&](std::string line) mutable {
        int *old_buf = buf;
        ksplit_core(&line[0], '\t', &buf_sz, &buf);
        return buf != old_buf;
    };

    buf = NULL;
    buf_sz = 0;
    // ksplit always realloc for offsets
    EXPECT_EQ(does_ksplit_cause_reallocation(lines[0]), true);
    EXPECT_EQ(does_ksplit_cause_reallocation(lines[1]), true);
    EXPECT_EQ(does_ksplit_cause_reallocation(lines[2]), true);
    EXPECT_EQ(does_ksplit_cause_reallocation(lines[3]), true);
    free(buf);

    buf = NULL;
    buf_sz = 0;
    // ksplit_core only realloc for offsets if larger size is needed
    EXPECT_EQ(does_ksplit_core_cause_reallocation(lines[0]), true);
    // EXPECT_EQ(does_ksplit_core_cause_reallocation(lines[1]), false);
    // EXPECT_EQ(does_ksplit_core_cause_reallocation(lines[2]), true);
    EXPECT_EQ(does_ksplit_core_cause_reallocation(lines[3]), false);
    free(buf);
}

// Much attention to horizontal region where bp increases without change in cM.
TEST(ibdtools, GeneticMap)
{
    GeneticMap gmap(-1, map_fn);

    size_t first_bp = gmap.get_first_nonzero_bp();
    size_t last_bp = gmap.get_last_bp();
    // long double first_cm = gmap.get_first_nonzero_cm();
    // long double last_cm = gmap.get_last_cm();

    // std::cout << first_bp << ',' << last_bp << ',' << first_cm << ',' << last_cm <<
    // '\n';

    size_t bp;
    long double cm;

    bp = first_bp;
    cm = gmap.get_cm(bp);
    if (!gmap.isin_horizontal_region(cm)) {
        EXPECT_EQ(bp, gmap.get_bp(cm));
    }

    bp = last_bp;
    cm = gmap.get_cm(bp);
    if (!gmap.isin_horizontal_region(cm)) {
        EXPECT_EQ(bp, gmap.get_bp(cm));
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> unif(0);

    for (size_t i = 0; i < 20; i++) {

        bp = unif(gen) % first_bp;
        cm = gmap.get_cm(bp);
        if (!gmap.isin_horizontal_region(cm)) {
            EXPECT_EQ(bp, gmap.get_bp(cm));
        }

        bp = unif(gen) % (last_bp - first_bp) + first_bp;
        cm = gmap.get_cm(bp);
        if (!gmap.isin_horizontal_region(cm)) {
            EXPECT_EQ(bp, gmap.get_bp(cm));
        }

        while ((bp = unif(gen) % 300 * 1024 * 1024) <= last_bp)
            ;
        cm = gmap.get_cm(bp);
        if (!gmap.isin_horizontal_region(cm)) {
            EXPECT_EQ(bp, gmap.get_bp(cm));
        }
    }
}

TEST(ibdtools, Meta_read_write)
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);

    Samples &sam1 = meta.get_samples();
    const char *temp_file = "tmp.gz";
    BGZF *fp = bgzf_open(temp_file, "w");
    sam1.write_to_file(fp);
    bgzf_close(fp);

    Samples sam2;
    fp = bgzf_open(temp_file, "r");
    sam2.read_from_file(fp);
    bgzf_close(fp);

    EXPECT_EQ(true, meta.get_samples().is_equal(sam2));

    Chromosomes chrs1;
    fp = bgzf_open(temp_file, "w");
    chrs1.write_to_file(fp);
    bgzf_close(fp);

    Chromosomes chrs2;
    fp = bgzf_open(temp_file, "r");
    chrs2.read_from_file(fp);
    bgzf_close(fp);

    EXPECT_EQ(true, chrs1.is_equal(chrs2));

    Positions &pos1 = meta.get_positions();
    fp = bgzf_open(temp_file, "w");
    pos1.write_to_file(fp);
    bgzf_close(fp);

    Positions pos2(-1);
    fp = bgzf_open(temp_file, "r");
    pos2.read_from_file(fp);
    bgzf_close(fp);

    EXPECT_EQ(true, meta.get_positions().is_equal(pos2));

    fp = bgzf_open(temp_file1, "w");
    bgzf_mt(fp, 10, 256);
    meta.write_to_file(fp);
    bgzf_close(fp);

    MetaFile meta2;
    fp = bgzf_open(temp_file1, "r");
    meta2.read_from_file(fp);
    bgzf_close(fp);

    EXPECT_EQ(true, meta.is_equal(meta2));
}

TEST(ibdtools, Genotypes_get_haplotypes)
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);
    auto haps1 = meta.get_genotypes().get_haplotypes();
    auto haps2 = read_lines_from_file(haps_fn);

    EXPECT_EQ(haps1.size(), haps2.size());
    EXPECT_EQ(true, (haps1 == haps2));
}

TEST(ibdtools, IbdFile_encode_raw_ibd)
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);

    IbdFile ibdfile1(temp_file1, &meta);
    ibdfile1.open("w");
    ibdfile1.from_raw_ibd(ibd_txt_fn);
    ibdfile1.close();

    IbdFile ibdfile2(temp_file1, &meta);
    ibdfile2.open("r");
    ibdfile2.read_from_file();
    ibdfile2.close();

    EXPECT_EQ(ibdfile1.has_equal_vector(ibdfile2), true);
}

TEST(ibdtools, IbdFile_decode_packed_ibd)
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true, "5");

    IbdFile ibdfile1(temp_file1, &meta);
    ibdfile1.open("w");
    ibdfile1.from_raw_ibd(ibd_txt_fn); // encode from raw ibd file to temp_file
    ibdfile1.close();

    IbdFile ibdfile2(temp_file1, &meta);
    ibdfile2.open("r");
    ibdfile2.to_raw_ibd(temp_file2); // decode to temp2_file
    ibdfile2.close();

    IbdFile ibdfile3(temp_file1, &meta);
    ibdfile3.open("w");
    ibdfile3.from_raw_ibd(temp_file2); // encode again
    ibdfile3.close();

    EXPECT_EQ(true, ibdfile1.has_equal_vector(ibdfile3));
}

TEST(ibdtools, Positions_get_bp_get_id)
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);
    auto &pos = meta.get_positions();

    size_t bp = pos.get_bp(52400);
    size_t id = pos.get_id(bp);
    EXPECT_EQ(id, 52400);

    bp = pos.get_bp(52398);
    id = pos.get_id(bp);
    EXPECT_EQ(id, 52398);
}

TEST(ibdtools, GziFile_caddr_and_uaddr)
{
    IbdFile ibdfile(temp_file1, NULL, 100 * 1024);
    auto &vec = ibdfile.get_vec();
    vec.resize(vec.capacity());

    std::string fn = temp_file1;
    fn += ".gzi";

    GziFile gzif(fn.c_str());
    for (auto &e : gzif.vec)
        EXPECT_TRUE(e.caddr <= e.uaddr && e.uaddr <= 64 * 1024);
}

TEST(ibdtools, StringViewSplitter)
{
    const char *line1 = "1\t2\t\t3\t4";
    const char *line2 = "0\t1\t2\t\t5\t4\t";
    const char *delim = "\t";
    StringViewSplitter svs(delim);

    svs.split(line1);

    std::string field;
    std::string_view field_sv;
    int field_int;
    svs.get(2, field);
    svs.get(3, field_sv);
    svs.get(0, field_int);

    EXPECT_EQ(field, "");
    EXPECT_EQ(field_sv, "3");
    EXPECT_EQ(field_int, 1);

    std::vector<int> vec;
    svs.split_to_vector(line1, vec);
    EXPECT_TRUE((vec == std::vector<int>{ 1, 2, 0, 3, 4 }));

    svs.split_to_vector(line2, vec);
    EXPECT_TRUE((vec == std::vector<int>{ 0, 1, 2, 0, 5, 4 }));

    svs.split_to_vector(line1, vec, true);  // clear
    svs.split_to_vector(line2, vec, false); // append
    EXPECT_TRUE((vec == std::vector<int>{ 1, 2, 0, 3, 4, 0, 1, 2, 0, 5, 4 }));

    svs.split_to_vector(line1, vec, true, 3);  // clear, only extract 3 columns
    svs.split_to_vector(line2, vec, false, 3); // append, only extract 3 columns
    EXPECT_TRUE((vec == std::vector<int>{ 1, 2, 0, 0, 1, 2 }));
}

TEST(cpp, test_from_chars_zero_length)
{
    std::string_view sv("1\t2\t3\t4");
    std::string_view sv2 = sv.substr(1, 0);
    int val = 100;
    std::from_chars(sv2.begin(), sv2.end(), val);
    EXPECT_FALSE((val == 1));
    EXPECT_FALSE((val == 0));
    EXPECT_TRUE((val == 100));
}

TEST(ibdtools, TournamentTree)
{

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint16_t> unif(0);

    std::vector<std::vector<uint16_t> > channels;
    std::vector<uint16_t> initials, out, all_nums;

    channels.resize(5);

    TournamentTree<uint16_t> ttree(5, 255);

    for (auto &v : channels) {
        for (size_t i = 0; i < 20; i++) {
            v.push_back(unif(gen) % 199); // fill in randome number less than max
            all_nums.push_back(v.back()); // save this for comparison
        }
        std::sort(v.begin(), v.end()); // each source input itself is sorted
        EXPECT_TRUE(v.end() == std::find(v.begin(), v.end(), 255));
        v.push_back(255); // add an sentinel

        // std::cout << "v: ";
        // std::copy(v.begin(), v.end(), std::ostream_iterator<uint16_t>(std::cout, "
        // ")); std::cout << "\n";
    }

    for (auto &v : channels) {
        initials.push_back(v.front());
        v.erase(v.begin());
    }

    // std::cout << "initials: ";
    // std::copy(initials.begin(), initials.end(),
    //     std::ostream_iterator<uint16_t>(std::cout, " "));
    // std::cout << "\n ";

    size_t winner_id;
    out.push_back(ttree.init_run(initials, winner_id));

    // std::cout << "resutls back: " << ((uint16_t) results.back()) << "\n";

    while (out.back() != 255) {
        auto &input = channels[winner_id];
        out.push_back(ttree.replace_run(input.front(), winner_id));

        // std::cout << "resutls back: " << ((uint16_t) results.back()) << "\n";
        //
        input.erase(input.begin());
    }

    out.pop_back(); // remove the last sentinel;

    // the out should be sorted version of all_numbers
    std::sort(all_nums.begin(), all_nums.end());
    EXPECT_EQ(out.size(), all_nums.size());
    EXPECT_TRUE((out == all_nums));
}

// Two times difference. Disabled to save big io cost
TEST(htslib, DISABLED_bgzf_write_uncompressed)
{
    std::vector<int> vec;
    size_t sz = 100 * 1024 * 1024;
    vec.resize(sz);
    ssize_t res;

    BGZF *bg_fp = bgzf_open("tmp_write_test.gz", "wu");
    res = bgzf_write(bg_fp, &vec[0], sz * sizeof(int));
    EXPECT_TRUE(res > 0 && (size_t) (res) == sz * sizeof(int));
    bgzf_close(bg_fp);

    BGZF *fp = bgzf_open("tmp_write_test.gz", "ru");
    assert(sz * sizeof(int) == bgzf_read(fp, &vec[0], sz * sizeof(int)));
    bgzf_close(fp);
}

// Two times difference. Disabled to save big io cost
TEST(htslib, DISABLED_bgzf_write_uncompressed_vs_fwrite)
{
    std::vector<int> vec;
    size_t sz = 100 * 1024 * 1014;
    vec.resize(sz);
    ssize_t res;

    FILE *bg_fp = fopen("tmp_write_test.gz", "wu");
    res = fwrite(&vec[0], sizeof(int), sz, bg_fp);
    EXPECT_TRUE((res > 0) && ((size_t) res) == sz);
    fclose(bg_fp);

    FILE *fp = fopen("tmp_write_test.gz", "ru");
    res = fread(&vec[0], sizeof(int), sz, fp);
    EXPECT_TRUE((res > 0) && ((size_t) res) == sz);
    fclose(fp);
}

TEST(cpp, default_move_constructor)
{
    class B
    {
      public:
        const char *str;
        std::vector<int> vec;
        FILE *fp;
    };

    class A
    {
      public:
        const char *str;
        std::vector<int> vec;
        FILE *fp;
        B *b;
    };

    B b;
    b.str = "b_str";
    b.vec.push_back(123);
    b.vec.push_back(456);

    A a1;
    a1.b = &b;
    a1.str = "A1_str";
    a1.vec.push_back(9999);

    FILE *fp = fopen(vcf_fn, "r");
    a1.fp = fp;

    A a2(std::move(a1));

    EXPECT_TRUE(strcmp(a2.b->str, "b_str") == 0);

    char buffer[1000];
    char *res_p = fgets(buffer, 1000, a2.fp);

    EXPECT_EQ(res_p, buffer);

    fclose(a2.fp);
}

TEST(ibdtools, IbdSorter)
{
    // meta
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);

    // from TXT (ibd_txt_fn) --> Binary (temp_file1)
    IbdFile ibdfile1(temp_file1, &meta);
    ibdfile1.open("wu");
    ibdfile1.from_raw_ibd(ibd_txt_fn);
    ibdfile1.close();

    // Sort: from Binary(temp_file1) to Binary(temp_file2)
    std::string out_prefix = temp_file2;
    out_prefix += "_temp_";
    IbdSorter sorter(temp_file1, temp_file2, "wu", out_prefix.c_str(), 1000);
    sorter.sort_into_chunks();
    sorter.merge_chunks(3);

    // read sorted and unsorted file
    IbdFile file_orig(temp_file1, NULL);
    IbdFile file_sort(temp_file2, NULL);

    file_orig.open("r");
    file_sort.open("r");

    file_orig.read_from_file();
    file_sort.read_from_file();

    file_orig.close();
    file_sort.close();

    // sort the unsored record in memory
    auto &orig_vec = file_orig.get_vec();
    std::sort(orig_vec.begin(), orig_vec.end());

    EXPECT_TRUE((orig_vec == file_sort.get_vec()));
}

TEST(ibdtools, IbdMerger)
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true, "10");
    BGZF *fp = bgzf_open("tmp_meta.gz", "w");
    meta.write_to_file(fp);
    bgzf_close(fp);

    IbdFile ibdfile1(temp_file1, &meta);
    ibdfile1.open("wu");
    ibdfile1.from_raw_ibd(ibd_txt_fn);
    ibdfile1.close();

    std::string out_prefix = temp_file2;
    out_prefix += "_temp_";
    IbdSorter sorter(temp_file1, temp_file2, "wu", out_prefix.c_str(), 1000);
    sorter.sort_into_chunks();
    sorter.merge_chunks(3);

    // read merged file
    MetaFile meta2;
    meta.parse_files(vcf_fn, map_fn, true, "10");
    IbdFile sorted_here(temp_file2, &meta2);
    sorted_here.open("r");
    sorted_here.read_from_file(false);
    sorted_here.close();

    IbdMerger merger(temp_file2, temp_file3, "w", "tmp_meta.gz", 100);
    merger.merge();

    // read merged file
    IbdFile merged_here(temp_file3, &meta);
    merged_here.open("r");
    merged_here.read_from_file(false);
    merged_here.close();

    // encode from brownning's merge txt file
    IbdFile merged_browning("/dev/null", &meta);
    merged_browning.open("w");
    merged_browning.from_raw_ibd(browning_merged);
    merged_browning.close();
    // set hap bits to 0
    for (auto &x : merged_browning.get_vec()) {
        x.hid1 = 0;
        x.hid2 = 0;
    }
    auto &vec = merged_browning.get_vec();
    std::sort(vec.begin(), vec.end());

    // compare
    EXPECT_TRUE(merged_here.get_vec() == merged_browning.get_vec());
}

TEST(cpp, bit_fields)
{
    struct A {
        uint64_t a : 1, b : 1, c : 1, d : 1, e : 60;
    };

    uint64_t x = ~uint64_t(0);
    A a = { 0 };

    a.a = x;

    EXPECT_EQ(a.a, 1);
    EXPECT_EQ(a.e, 0);

    auto y = a.a;
    EXPECT_EQ(sizeof(y), sizeof(uint64_t));

    struct __attribute__((packed)) B {
        uint16_t x;
        uint64_t a : 16, b : 16, c : 16, d : 16;
    };
    B b;
    b.x = 0x1;
    b.a = 0x2;
    b.b = 0x3;
    b.c = 0x4;
    b.d = 0x5;

    uint16_t *p = (uint16_t *) ((B *) &b);

    EXPECT_EQ(p[0], 1);
    EXPECT_EQ(p[1], 2);
    EXPECT_EQ(p[2], 3);
    EXPECT_EQ(p[3], 4);
    EXPECT_EQ(p[4], 5);
}

// some definition
struct seg_t {
    uint32_t pid_l;
    uint32_t pid_r;
    uint32_t label;
    friend bool
    operator==(const seg_t &s1, const seg_t &s2)
    {
        return s1.pid_l == s2.pid_l && s1.pid_r == s2.pid_r && s1.label == s2.label;
    };
    friend bool
    operator<(const seg_t &s1, const seg_t &s2)
    {
        return s1.pid_l < s2.pid_r;
    }
};

TEST(ibdtools, IbdSplitter)
{

    {
        // make meta file
        Positions pos(0);
        for (int i = 0; i <= 11; i++)
            pos.add(i, i * 1.0);

        MetaFile meta;
        meta.get_positions() = pos;
        BGZF *fp = bgzf_open("tmp_meta.gz", "w");
        meta.write_to_file(fp);
        bgzf_close(fp);

        // make region labels
        auto labels = std::vector<region_label_t>{ { 0, 1 }, { 1, 0 }, { 2, 1 },
            { 5, 0 }, { 8, 1 } };

        using SegV = std::vector<seg_t>;

        // helper function
        auto split = [&labels](seg_t seg) -> SegV {
            // make a record
            SegV seg_vec;
            ibd_rec2_t r2{ 1, 0, 1, 10, 0, 1 };
            r2.pid1 = seg.pid_l;
            r2.pid2 = seg.pid_r;

            // Make ibd binary file contains one record
            IbdFile ibd(temp_file1, NULL, 10);
            ibd.get_vec().push_back(r2);
            ibd.open("w");
            ibd.write_to_file();
            ibd.close();

            // run splitter
            IbdSplitter splitter(temp_file1, "tmp_", "tmp_meta.gz", labels);
            splitter.split();

            // read in file 0
            uint32_t label = 0;
            IbdFile ibd2("tmp_0", NULL, 100);
            ibd2.open("r");
            ibd2.read_from_file();
            ibd2.close();

            for (auto rec : ibd2.get_vec())
                seg_vec.push_back({ rec.get_pid1(), rec.get_pid2(), label });

            // read in file 1
            label = 1;
            IbdFile ibd3("tmp_1", NULL, 100);
            ibd3.open("r");
            ibd3.read_from_file();
            ibd3.close();

            for (auto rec : ibd3.get_vec())
                seg_vec.push_back(seg_t{ rec.get_pid1(), rec.get_pid2(), label });

            std::sort(seg_vec.begin(), seg_vec.end());

            return seg_vec;
        };

        // auto print_segv = [](SegV seg_vec) {
        //     std::cout << "Print_seg_start \n";
        //     for (auto seg : seg_vec) {
        //         std::cout << "Seg: (" << seg.pid_l << ", " << seg.pid_r
        //                   << "), label = " << seg.label << '\n';
        //     }
        //     std::cout << "Print_seg_end \n";
        // };

        EXPECT_TRUE((split({ 0, 1 }) == SegV{}));
        EXPECT_TRUE((split({ 2, 4 }) == SegV{ { 2, 4, 1 } }));
        EXPECT_TRUE((split({ 5, 7 }) == SegV{ { 5, 7, 0 } }));
        EXPECT_TRUE((split({ 0, 2 }) == SegV{}));
        EXPECT_TRUE((split({ 3, 4 }) == SegV{}));
        EXPECT_TRUE((split({ 0, 3 }) == SegV{}));
        // print_segv(split({ 0, 3 }));
        EXPECT_TRUE((split({ 0, 4 }) == SegV{ { 2, 4, 1 } }));
        EXPECT_TRUE((split({ 0, 5 }) == SegV{ { 2, 5, 1 } }));
        EXPECT_TRUE((split({ 3, 7 }) == SegV{ { 3, 5, 1 }, { 5, 7, 0 } }));
        EXPECT_TRUE(
            (split({ 2, 11 }) == SegV{ { 2, 5, 1 }, { 5, 8, 0 }, { 8, 11, 1 } }));
        EXPECT_TRUE((split({ 4, 9 }) == SegV{ { 5, 8, 0 } }));
    }

    {
        // just make sure the program run with real data files
        MetaFile meta;
        meta.parse_files(vcf_fn, map_fn, true, "10");
        BGZF *fp = bgzf_open("tmp_meta.gz", "w");
        meta.write_to_file(fp);
        bgzf_close(fp);

        IbdFile ibdfile1(temp_file1, &meta);
        ibdfile1.open("wu");
        ibdfile1.from_raw_ibd(ibd_txt_fn);
        ibdfile1.close();

        auto labels = meta.get_positions().get_gap_vector(2.0, 100);
        std::string prefix = "tmp__splitted__";

        IbdSplitter splitter(
            temp_file1, prefix.c_str(), "tmp_meta.gz", labels, 1, 2.0, 1000, 1000);

        splitter.split();

        EXPECT_TRUE(std::filesystem::exists(prefix + '0'));
        EXPECT_TRUE(std::filesystem::exists(prefix + '1'));

        IbdFile f1((prefix + '0').c_str(), &meta);
        f1.open("r");
        f1.to_raw_ibd((prefix + "0.txt.gz").c_str());
        f1.close();

        IbdFile f2((prefix + '1').c_str(), &meta);
        f2.open("r");
        f2.to_raw_ibd((prefix + "1.txt.gz").c_str());
        f2.close();
    }
}

TEST(math, low_triangular_index_conversion)
{
    IbdMatrix m;

    EXPECT_EQ(m.get_arr_index(1, 0), 0);
    EXPECT_EQ(m.get_arr_index(2, 0), 1);
    EXPECT_EQ(m.get_arr_index(2, 1), 2);
    EXPECT_EQ(m.get_arr_index(3, 0), 3);
    EXPECT_EQ(m.get_arr_index(3, 1), 4);
    EXPECT_EQ(m.get_arr_index(3, 2), 5);
    EXPECT_EQ(m.get_arr_index(4, 0), 6);
    EXPECT_EQ(m.get_arr_index(4, 1), 7);
    EXPECT_EQ(m.get_arr_index(4, 2), 8);
    EXPECT_EQ(m.get_arr_index(4, 3), 9);
    EXPECT_EQ(m.get_arr_index(5, 0), 10);
    EXPECT_EQ(m.get_arr_index(5, 1), 11);
    EXPECT_EQ(m.get_arr_index(5, 2), 12);
    EXPECT_EQ(m.get_arr_index(5, 3), 13);
    EXPECT_EQ(m.get_arr_index(5, 4), 14);

    std::random_device ran;
    std::mt19937 gen(ran());
    std::uniform_int_distribution<uint32_t> unif(0, 100000);

    uint32_t row, col, r, c;
    for (uint32_t i = 0; i < 10000; i++) {
        row = unif(gen);
        col = unif(gen);
        if (row == col)
            row = col + 1;
        else if (row < col)
            std::swap(row, col);

        size_t arr_index = m.get_arr_index(row, col);

        m.get_matrix_index(arr_index, r, c);
        EXPECT_EQ(r, row);
        EXPECT_EQ(c, col);
    }
}

TEST(ibdtools, IbdMatrix)
{
    // run through
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, 0, "10");

    IbdFile in(temp_file1, &meta);
    in.open("w");
    in.from_raw_ibd(ibd_txt_fn);
    in.close();

    IbdMatrix mat_sam, mat_hap;

    in.open("r");
    mat_sam.calculate_total_from_ibdfile(in);
    in.close();
    mat_sam.write_matrix_file(temp_file2);
    mat_sam.read_matrix_file(temp_file2);

    in.open("r");
    mat_hap.calculate_total_from_ibdfile(in, true);
    in.close();
    mat_hap.write_matrix_file(temp_file3);
    mat_hap.read_matrix_file(temp_file3);

    // exact compair 1
    in.open("r");
    in.read_from_file();
    in.close();
    auto &vec = in.get_vec();
    auto &pos = meta.get_positions();

    std::random_device ran;
    std::mt19937 gen(ran());
    std::uniform_int_distribution<uint32_t> unif(0, mat_sam.get_num_samples() - 1);

    uint32_t row, col;
    for (uint32_t i = 0; i < 10000; i++) {
        row = unif(gen);
        col = unif(gen);
        if (row == col) {
            continue;
        } else if (row < col)
            std::swap(row, col);
        else {
            // do nothing
        }

        uint16_t total = 0;
        uint16_t total_hap[4] = { 0, 0, 0, 0 };
        for (auto rec : vec) {
            if (rec.get_sid1() == row && rec.get_sid2() == col) {
                float cm = pos.get_cm(rec.get_pid2()) - pos.get_cm(rec.get_pid1());
                uint16_t cmx10 = lround(cm * 10);
                total += cmx10;
                if (rec.get_hid1() == 0) {
                    if (rec.get_hid2() == 0) {
                        total_hap[0] += cmx10;
                    } else {
                        total_hap[1] += cmx10;
                    }
                } else {
                    if (rec.get_hid2() == 0) {
                        total_hap[2] += cmx10;
                    } else {
                        total_hap[3] += cmx10;
                    }
                }
            }
        }

        EXPECT_EQ(mat_sam.at(row, col), total);
        EXPECT_EQ(mat_hap.at((row << 1) + 0, (col << 1) + 0), total_hap[0]);
        EXPECT_EQ(mat_hap.at((row << 1) + 0, (col << 1) + 1), total_hap[1]);
        EXPECT_EQ(mat_hap.at((row << 1) + 1, (col << 1) + 0), total_hap[2]);
        EXPECT_EQ(mat_hap.at((row << 1) + 1, (col << 1) + 1), total_hap[3]);
    }
}

TEST(ibdtools, IbdCoverage_run_thru)
{

    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true, "10");
    BGZF *fp = bgzf_open("tmp_meta.gz", "w");
    meta.write_to_file(fp);
    bgzf_close(fp);

    IbdFile in(temp_file1, &meta);
    in.open("w");
    in.from_raw_ibd(ibd_txt_fn);
    in.close();

    IbdCoverage cov(temp_file1, "tmp_meta.gz", 1, 1000);
    cov.calculate_coverage();

    cov.write_to_file(temp_file2);

    IbdCoverage cov2;
    cov2.read_from_file(temp_file2);
}

TEST(ibdtools, IbdCoverage_exact)
{

    // make meta file
    Positions pos(0);
    for (int i = 0; i <= 11; i++)
        pos.add(i, i * 1.0);

    MetaFile meta;
    meta.get_positions() = pos;
    BGZF *fp = bgzf_open("tmp_meta.gz", "w");
    meta.write_to_file(fp);
    bgzf_close(fp);

    // make ibd vec
    using namespace std;
    auto pairs
        = vector<pair<uint32_t, uint32_t> >{ { 0, 1 }, { 2, 4 }, { 5, 7 }, { 0, 2 },
              { 3, 4 }, { 0, 3 }, { 0, 4 }, { 0, 5 }, { 3, 7 }, { 2, 11 }, { 4, 9 } };

    std::vector<ibd_rec1_t> rec_vec;
    for (auto pair : pairs) {
        ibd_rec2_t r2{ 1, 0, 1, 10, 0, 1 };
        r2.set_pid1(get<0>(pair));
        r2.set_pid2(get<1>(pair));
        rec_vec.push_back(r2);
    }

    // write to ibd file
    IbdFile ibd(temp_file1, NULL, 10);
    ibd.get_vec() = rec_vec;
    ibd.open("w");
    ibd.write_to_file();
    ibd.close();

    IbdCoverage cov(temp_file1, "tmp_meta.gz");
    cov.calculate_coverage();
    cov.write_to_file(temp_file2);
    // cov.summary(std::cout);

    IbdCoverage cov2;
    cov2.read_from_file(temp_file2);

    vector<float> expected_cm_vec;
    for (int i = 0; i <= 12; i++)
        expected_cm_vec.push_back(i * 1.0);

    EXPECT_EQ(cov2.get_cm_vec(), expected_cm_vec);

    EXPECT_EQ(cov2.get_count_vec()[0], 5);
    EXPECT_EQ(cov2.get_count_vec()[1], 4);
    EXPECT_EQ(cov2.get_count_vec()[2], 5);
    EXPECT_EQ(cov2.get_count_vec()[3], 6);
    EXPECT_EQ(cov2.get_count_vec()[4], 4);
    EXPECT_EQ(cov2.get_count_vec()[5], 4);
    EXPECT_EQ(cov2.get_count_vec()[6], 4);
    EXPECT_EQ(cov2.get_count_vec()[7], 2);
    EXPECT_EQ(cov2.get_count_vec()[8], 2);
    EXPECT_EQ(cov2.get_count_vec()[9], 1);
    EXPECT_EQ(cov2.get_count_vec()[10], 1);
    EXPECT_EQ(cov2.get_count_vec()[11], 0);
    EXPECT_EQ(cov2.get_count_vec()[12], 0);
}

TEST(ibdtools, bgzidx)
{

    {
        MetaFile meta;
        meta.parse_files(vcf_fn, map_fn, true);
        BGZF *fp = bgzf_open(temp_file1, "w");
        bgzf_index_build_init(fp);
        meta.write_to_file(fp);
        assert(0 == bgzf_index_dump(fp, temp_file1, ".gzi"));
        bgzf_close(fp);
    }
    {
        BGZF *fp = bgzf_open(temp_file1, "r");
        assert(0 == bgzf_index_load(fp, temp_file1, ".gzi"));
        __used(fp);
    }
}

TEST(ibdtools, add_excluded_region)
{
    std::vector<region_label_t> label_v_orig
        = { { 1, 0 }, { 2, 1 }, { 5, 0 }, { 11, 1 }, { 15, 0 }, { 18, 1 } };

    std::vector<region_label_t> expected_v, result_v;

    result_v = label_v_orig;
    add_exclusion_range(result_v, 3, 4);
    expected_v = { { 1, 0 }, { 2, 1 }, { 3, 0 }, { 4, 1 }, { 5, 0 }, { 11, 1 },
        { 15, 0 }, { 18, 1 } };
    EXPECT_EQ(result_v, expected_v);

    result_v = label_v_orig;
    add_exclusion_range(result_v, 2, 4);
    expected_v = { { 1, 0 }, { 4, 1 }, { 5, 0 }, { 11, 1 }, { 15, 0 }, { 18, 1 } };
    EXPECT_EQ(result_v, expected_v);

    result_v = label_v_orig;
    add_exclusion_range(result_v, 6, 8);
    expected_v = label_v_orig;
    EXPECT_EQ(result_v, expected_v);

    result_v = label_v_orig;
    add_exclusion_range(result_v, 2, 13);
    expected_v = { { 1, 0 }, { 13, 1 }, { 15, 0 }, { 18, 1 } };
    EXPECT_EQ(result_v, expected_v);

    result_v = label_v_orig;
    add_exclusion_range(result_v, 3, 13);
    expected_v = { { 1, 0 }, { 2, 1 }, { 3, 0 }, { 13, 1 }, { 15, 0 }, { 18, 1 } };
    EXPECT_EQ(result_v, expected_v);

    result_v = label_v_orig;
    add_exclusion_range(result_v, 3, 7);
    expected_v = { { 1, 0 }, { 2, 1 }, { 3, 0 }, { 11, 1 }, { 15, 0 }, { 18, 1 } };
    EXPECT_EQ(result_v, expected_v);

    result_v = label_v_orig;
    add_exclusion_range(result_v, 2, 7);
    expected_v = { { 1, 0 }, { 11, 1 }, { 15, 0 }, { 18, 1 } };
    EXPECT_EQ(result_v, expected_v);

    result_v = label_v_orig;
    add_exclusion_range(result_v, 6, 15);
    expected_v = { { 1, 0 }, { 2, 1 }, { 5, 0 }, { 18, 1 } };
    EXPECT_EQ(result_v, expected_v);

    result_v = label_v_orig;
    add_exclusion_range(result_v, 6, 14);
    expected_v = { { 1, 0 }, { 2, 1 }, { 5, 0 }, { 14, 1 }, { 15, 0 }, { 18, 1 } };
    EXPECT_EQ(result_v, expected_v);

    result_v = label_v_orig;
    add_exclusion_range(result_v, 6, 11);
    expected_v = { { 1, 0 }, { 2, 1 }, { 5, 0 }, { 11, 1 }, { 15, 0 }, { 18, 1 } };
    EXPECT_EQ(result_v, expected_v);

    result_v = label_v_orig;
    add_exclusion_range(result_v, 6, 20);
    expected_v = { { 1, 0 }, { 2, 1 }, { 5, 0 }, { 20, 1 } };
    EXPECT_EQ(result_v, expected_v);
}

int
main(int argc, char **argv)
{
    bool targeted_test = false;

    if (targeted_test) {
        int myargc = 2;
        char *myargv[2];
        char arg1[] = "./gtest";
        char arg2[] = "--gtest_filter=ibdtools.add_excluded_region";
        myargv[0] = arg1;
        myargv[1] = arg2;
        ::testing::InitGoogleTest(&myargc, myargv);
    } else {
        ::testing::InitGoogleTest(&argc, argv);
    }

    return RUN_ALL_TESTS();
}
