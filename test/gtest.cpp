#include "../include/ibdfile.hpp"
#include "../include/ibdmerger.hpp"
#include "../include/ibdsorter.hpp"
#include "../include/metafile.hpp"
#include <algorithm>
#include <cstdint>
#include <gtest/gtest.h>
#include <iterator>
#include <memory>
#include <random>
#include <stdio.h>

const char *map_fn = "../data/example.map";
const char *vcf_fn = "../data/example.bcf.gz";
const char *haps_fn = "../data/example_haplotypes.txt.gz";
const char *ibd_txt_fn = "../data/example.ibd.gz";
const char *browning_merged = "../data/merged_browning.ibd.gz";
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
        int acutal_num_fields = ksplit_core(&line[0], '\t', &buf_sz, &buf);
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
    EXPECT_EQ(does_ksplit_core_cause_reallocation(lines[1]), false);
    EXPECT_EQ(does_ksplit_core_cause_reallocation(lines[2]), true);
    EXPECT_EQ(does_ksplit_core_cause_reallocation(lines[3]), false);
    free(buf);
}

// Much attention to horizontal region where bp increases without change in cM.
TEST(ibdtools, GeneticMap)
{
    GeneticMap gmap(-1, map_fn);

    size_t first_bp = gmap.get_first_nonzero_bp();
    size_t last_bp = gmap.get_last_bp();
    long double first_cm = gmap.get_first_nonzero_cm();
    long double last_cm = gmap.get_last_cm();

    // std::cout << first_bp << ',' << last_bp << ',' << first_cm << ',' << last_cm <<
    // '\n';

    size_t bp;
    long double cm;

    bp = first_bp;
    cm = gmap.get_cm(bp);
    if (!gmap.isin_horizontal_region(cm))
        EXPECT_EQ(bp, gmap.get_bp(cm));

    bp = last_bp;
    cm = gmap.get_cm(bp);
    if (!gmap.isin_horizontal_region(cm))
        EXPECT_EQ(bp, gmap.get_bp(cm));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> unif(0);

    for (size_t i = 0; i < 20; i++) {

        bp = unif(gen) % first_bp;
        cm = gmap.get_cm(bp);
        if (!gmap.isin_horizontal_region(cm))
            EXPECT_EQ(bp, gmap.get_bp(cm));

        bp = unif(gen) % (last_bp - first_bp) + first_bp;
        cm = gmap.get_cm(bp);
        if (!gmap.isin_horizontal_region(cm))
            EXPECT_EQ(bp, gmap.get_bp(cm));

        while ((bp = unif(gen) % 300 * 1024 * 1024) <= last_bp)
            ;
        cm = gmap.get_cm(bp);
        if (!gmap.isin_horizontal_region(cm))
            EXPECT_EQ(bp, gmap.get_bp(cm));
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

TEST(htslib, bgzf_write_uncompressed)
{
    std::vector<int> vec;
    size_t sz = 100 * 1024 * 1014;
    vec.resize(sz);
    ssize_t res;

    BGZF *bg_fp = bgzf_open("1.gz", "wu");
    res = bgzf_write(bg_fp, &vec[0], sz * sizeof(int));
    assert(res == sz * sizeof(int));
    bgzf_close(bg_fp);

    BGZF *fp = bgzf_open("1.gz", "ru");
    assert(sz * sizeof(int) == bgzf_read(fp, &vec[0], sz * sizeof(int)));
    bgzf_close(fp);
}

TEST(htslib, bgzf_write_uncompressed_vs_fwrite)
{
    std::vector<int> vec;
    size_t sz = 100 * 1024 * 1014;
    vec.resize(sz);
    ssize_t res;

    FILE *bg_fp = fopen("1.gz", "wu");
    res = fwrite(&vec[0], sizeof(int), sz, bg_fp);
    assert(res == sz);
    fclose(bg_fp);

    int val = 0;
    FILE *fp = fopen("1.gz", "ru");
    res = fread(&vec[0], sizeof(int), sz, fp);
    assert(res == sz);
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
    BGZF *fp = bgzf_open("meta.gz", "w");
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

    IbdMerger merger(temp_file2, temp_file3, "w", "meta.gz", 100);
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

void
test_auto_variable()
{
    class A
    {
        std::vector<int> vec;

      public:
        A() { vec.resize(10, 10); }
        std::vector<int> &
        get_vec()
        {
            return vec;
        }
    };

    A a;

    auto x = a.get_vec();
    std::fill(x.begin(), x.end(), 20);
    auto &y = a.get_vec();
    std::cout << std::addressof(x) << " " << std::addressof(y) << '\n';

    std::copy(x.begin(), x.end(), std::ostream_iterator<int>(std::cout, " "));
    std::copy(y.begin(), y.end(), std::ostream_iterator<int>(std::cout, " "));
}

void
test_bit_fields()
{
    struct A {
        uint64_t a : 1, b : 1, c : 1, d : 1, e : 60;
    };

    A a;
    uint64_t x = ~uint64_t(0);

    a.e = x;

    auto y = a.a;

    std::cout << a.a << '\n';
    std::cout << sizeof(A) << '\n';
    std::cout << a.e << '\n';

    struct __attribute__((packed)) B {
        uint16_t x;
        uint64_t a : 16, b : 16, c : 16, d : 16;
    };
    B b;
    b.x = 0x11;
    b.a = 0xa;
    b.b = 0xb;
    b.c = 0xc;
    b.d = 0xd;

    std::cout << std::hex << *((uint16_t *) &b) << '\n';
    std::cout << std::hex << *((uint16_t *) &b + 4) << '\n';
}

int
main(int argc, char **argv)
{
    // test_ibdrec();
    // test_bgzf_getline_ksplit();
    // test_gmap();
    // test_vcffile();
    // test_samples_read_write();
    // test_positions_map();
    // test_positions_read_write();
    // test_chromosomes_read_write();
    // test_get_haplotypes();
    // test_metafile_read_write();
    // test_GziFile();
    // test_bgzf_index();
    // test_string_view_splitter();
    // test_tournament_tree();
    // test_bgzf_write();
    // test_move_default_move_constructor();
    // test_ibdsorter();
    // test_auto_variable();
    // test_bit_fields();
    // test_ibdfile_decode_pakced_ibd();
    // test_ibdmerger();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
