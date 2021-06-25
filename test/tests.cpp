#include "../include/ibdfile.hpp"
#include "../include/ibdsorter.hpp"
#include "../include/metafile.hpp"
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <memory>
#include <stdio.h>

const char *map_fn = "../data/example.map";
const char *vcf_fn = "../data/example.bcf.gz";
const char *haps_fn = "../data/example_haplotypes.txt.gz";
const char *ibd_in_file = "../data/example.ibd.gz";
const char *temp_file = "tmp.gz";
const char *temp_file2 = "tmp2.gz";

// ksplit seems to alloc new ids buffer everytime. Bad;
// //
// ksplit_core seems to realloc ids if necessary. Good. The document was bad
// @ first: the char * pointer to line buffer
// @ second: the index buffer size.
// @ third: the index buffer pointer. If buffer size is 0, then this need to be NULL. If
// there are more fields then buffer size. It will realloc the buffer and update buffer
// size and buffer pointer.
// @ return: number of field (n_fields) found. Only the first n_fields of the index
// buffer are valid. Note: there is no way to only split and get only the first n
// fields. The buffer should be finally freed by user.
// //
// kstring_t is reused in getline. Similar to ids in ksplit_core
void
test_bgzf_getline_ksplit()
{
    const char *line = "1\t2\t3\t4\n";
    std::string line2 = "1\t\t2\t3\t4\t";
    line2 += line2;
    line2 += line2;
    line2.back() = '\n';

    BGZF *fp = bgzf_open(temp_file, "w");
    assert(bgzf_write(fp, line, strlen(line)) == strlen(line));
    assert(bgzf_write(fp, line, strlen(line)) == strlen(line));
    assert(bgzf_write(fp, &line2[0], line2.size()) == line2.size());
    assert(bgzf_write(fp, line, strlen(line)) == strlen(line));
    bgzf_close(fp);

    {
        int *ids = 0, count;
        int max_fix;
        kstring_t str;

        auto print_info = [&, line_no = 1]() mutable {
            ids = ksplit(&str, '\t', &count);
            std::cout << "------- ksplit Line " << line_no++ << " ---------- \t"
                      << "str.s @: " << (void *) &str.s[0] << '\t'
                      << "ids @: " << (void *) ids << "\n ";
            for (int i = 0; i < count; i++)
                std::cout << i << ":" << str.s[ids[i]] << '\t';
            std::cout << '\n';
        };

        str = { 0 };
        ids = 0;
        fp = bgzf_open(temp_file, "r");
        assert(bgzf_getline(fp, '\n', &str) == strlen(line) - 1);
        print_info();
        assert(bgzf_getline(fp, '\n', &str) == strlen(line) - 1);
        print_info();
        assert(bgzf_getline(fp, '\n', &str) == line2.size() - 1);
        print_info();
        assert(bgzf_getline(fp, '\n', &str) == strlen(line) - 1);
        print_info();
        bgzf_close(fp);
        ks_free(&str);
        free(ids);
        ids = NULL;
    }

    {
        kstring_t str = { 0 };
        int buf_sz = 5;
        int *id_buf = (int *) malloc(sizeof(int) * buf_sz);

        auto print_info2 = [&, line_no = 1](int &max) mutable {
            int n_fields = ksplit_core(str.s, '\t', &max, &id_buf);
            std::cout << "------- ksplit_core max= " << max << " Line " << line_no
                      << " ---------- \t"
                      << "str.s @: " << (void *) str.s << '\t'
                      << "ids @: " << (void *) id_buf << "\n ";
            for (int i = 0; i < n_fields; i++)
                std::cout << i << ":" << str.s[id_buf[i]] << '\t';
            std::cout << '\n';
        };
        fp = bgzf_open(temp_file, "r");
        assert(bgzf_getline(fp, '\n', &str) == strlen(line) - 1);
        print_info2(buf_sz);
        assert(bgzf_getline(fp, '\n', &str) == strlen(line) - 1);
        print_info2(buf_sz);
        assert(bgzf_getline(fp, '\n', &str) == line2.size() - 1);
        print_info2(buf_sz);
        assert(bgzf_getline(fp, '\n', &str) == strlen(line) - 1);
        print_info2(buf_sz);
        bgzf_close(fp);

        ks_free(&str);
        free(id_buf);
        id_buf = NULL;
        buf_sz = 0;
    }
}
void
test_gmap()
{
    GeneticMap gmap(-1, map_fn);
    gmap.print();

    std::cout << "bp = 50,000,000, cm = " << gmap.get_cm(1000 * 1000 * 50) << '\n';
    std::cout << "cm = 8, bp = " << gmap.get_bp(8.0) << '\n';
}

void
test_vcffile()
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);
    // meta.print();
}

void
test_samples_read_write()
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);
    std::cerr << "parsed vcf\n";

    Samples &sam1 = meta.get_samples();
    const char *temp_file = "tmp.gz";
    BGZF *fp = bgzf_open(temp_file, "w");
    std::cerr << "write sample\n";
    sam1.write_to_file(fp);
    bgzf_close(fp);

    Samples sam2;
    fp = bgzf_open(temp_file, "r");
    std::cerr << "read sample\n";
    sam2.read_from_file(fp);
    bgzf_close(fp);

    std::cout << "Sample read_write is equal: " << meta.get_samples().is_equal(sam2)
              << '\n';
}

void
test_ibdrec()
{
    // test size
    std::vector<ibd_rec1_t> vec(10);
    std::sort(vec.begin(), vec.end());
    std::cout << sizeof(ibd_rec1_t) << '\n';
    std::cout << (char *) &vec[9] - (char *) &vec[0] << '\n';

    // test speed
    std::vector<ibd_rec1_t> vec1(1000000);
    std::vector<ibd_rec2_t> vec2(1000000);
    {
        ScopedTimer timer1("rec1");
        std::sort(vec1.begin(), vec1.end());
    }
    {
        ScopedTimer timer2("rec2");
        std::sort(vec2.begin(), vec2.end());
    }
}

void
test_chromosomes_read_write()
{
    const char *temp_file = "tmp.gz";
    Chromosomes chrs1;
    BGZF *fp = bgzf_open(temp_file, "w");
    chrs1.write_to_file(fp);
    bgzf_close(fp);

    Chromosomes chrs2;
    fp = bgzf_open(temp_file, "r");
    chrs2.read_from_file(fp);
    bgzf_close(fp);

    std::cout << "chromosomes_read_write is_equal: " << chrs1.is_equal(chrs2) << '\n';
}
void
test_positions_read_write()
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);
    std::cerr << "parsed vcf\n";

    Positions &pos1 = meta.get_positions();
    const char *temp_file = "tmp.gz";
    BGZF *fp = bgzf_open(temp_file, "w");
    std::cerr << "write positions\n";
    pos1.write_to_file(fp);
    bgzf_close(fp);

    Positions pos2(-1);
    fp = bgzf_open(temp_file, "r");
    std::cerr << "read positions\n";
    pos2.read_from_file(fp);
    bgzf_close(fp);

    std::cout << "Positions read_write is equal: " << meta.get_positions().is_equal(pos2)
              << '\n';
}

void
test_metafile_read_write()
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);
    std::cerr << "parsed vcf\n";

    BGZF *fp = bgzf_open(temp_file, "w");
    bgzf_mt(fp, 10, 256);
    bgzf_index_build_init(fp);
    assert(bgzf_index_dump(fp, temp_file, ".gzi") == 0);
    std::cerr << " write metafile\n";
    meta.write_to_file(fp);
    bgzf_close(fp);

    MetaFile meta2;
    fp = bgzf_open(temp_file, "r");
    bgzf_mt(fp, 10, 256);
    assert(bgzf_index_load(fp, temp_file, ".gzi") == 0);
    std::cerr << " read metafile\n";
    meta2.read_from_file(fp);
    bgzf_close(fp);

    std::cout << "MetaFile read_write is equal: " << meta.is_equal(meta2) << '\n';
}

void
test_get_haplotypes()
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);
    auto haps1 = meta.get_genotypes().get_haplotypes();
    auto haps2 = read_lines_from_file(haps_fn);

    std::cout << "haps1 == haps2 " << (haps1 == haps2) << '\n';
    std::cout << haps1[0] << '\n';
    std::cout << haps2[0] << '\n';
}

void
test_ibdfile_encode_raw_ibd()
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);

    // IbdFile ibdfile0(meta);
    // ibdfile.from_raw_ibd(ibd_in_file, temp_file, { 1000, 0, 1, 2, 3 });

    IbdFile ibdfile1(temp_file, &meta);
    ibdfile1.open("w");
    ibdfile1.from_raw_ibd(ibd_in_file);
    ibdfile1.close();
    // debug
    // return;
    IbdFile ibdfile2(temp_file, &meta);
    ibdfile2.open("r");
    ibdfile2.read_from_file();
    ibdfile2.close();

    std::cout << "ibdfile_encde_raw_ibd: is_equal? "
              << ibdfile1.has_equal_vector(ibdfile2) << '\n';
}

void
test_positions_map()
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);

    std::cout << "meta.get_positions().get_bp(52400): "
              << meta.get_positions().get_bp(52400) << '\n';
    std::cout << "meta.get_positions().get_id(45423941): "
              << meta.get_positions().get_id(45423941) << '\n';
    std::cout << "meta.get_positions().get_id(45416789): "
              << meta.get_positions().get_id(45416789) << '\n';
    std::cout << "meta.get_positions().get_bp(52398): "
              << meta.get_positions().get_bp(52398) << '\n';
}

void
test_ibdfile_decode_pakced_ibd()
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true, "5");

    IbdFile ibdfile1(temp_file, &meta);
    ibdfile1.open("w");
    ibdfile1.from_raw_ibd(ibd_in_file); // encode from raw ibd file to temp_file
    ibdfile1.close();

    IbdFile ibdfile2(temp_file, &meta);
    ibdfile2.open("r");
    ibdfile2.to_raw_ibd(temp_file2); // decode to temp2_file
    ibdfile2.close();

    IbdFile ibdfile3(temp_file, &meta);
    ibdfile3.open("w");
    ibdfile3.from_raw_ibd(temp_file2); // encode again
    ibdfile3.close();

    std::cout << "ibdfile_decode_packed_ibd is_equal? "
              << ibdfile1.has_equal_vector(ibdfile3) << '\n';
}

void
test_GziFile()
{
    std::string fn = "../data/vcf_pos.txt.gz.gzi";
    GziFile gzif(fn.c_str());
    std::cout << "number entries: " << gzif.vec.size() << '\n';
    for (auto &e : gzif.vec)
        std::cout << "caddr: " << e.caddr << ", uaddr: " << e.uaddr << '\n';
}

void
test_bgzf_index()
{
    const size_t sz = 20 * 1020;
    char buffer[sz];
    for (size_t i = 0; i < sz; i++)
        buffer[i] = rand() % 26 + 'a';

    BGZF *fp = bgzf_open(temp_file2, "w");
    bgzf_index_build_init(fp);
    size_t bytes_written = 0;
    for (int i = 0; i < 20; i++) {
        std::cout << "blockaddr: " << fp->block_address
                  << "\t blockoffset: " << fp->block_offset
                  << "\t bytes_written: " << bytes_written << '\n';
        assert(bgzf_write(fp, buffer, sz) == sz);
        bytes_written += sz;
    }
    bgzf_close(fp);
}

void
test_string_view_splitter()
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

    std::vector<int> vec;
    svs.split_to_vector(line1, vec);
    svs.split_to_vector(line2, vec);

    std::copy(vec.begin(), vec.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << '\n';

    svs.split_to_vector(line1, vec, true);  // clear
    svs.split_to_vector(line2, vec, false); // append

    std::copy(vec.begin(), vec.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << '\n';

    svs.split_to_vector(line1, vec, true, 3);  // clear, only extract 3 columns
    svs.split_to_vector(line2, vec, false, 3); // append, only extract 3 columns

    std::copy(vec.begin(), vec.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << '\n';
}

void
test_from_chars_zero_length()
{
    std::string_view sv("1\t2\t3\t4");
    std::string_view sv2 = sv.substr(1, 0);
    int val = 100;
    std::from_chars(sv2.begin(), sv2.end(), val);
    std::cout << val << '\n';
    std::cout << sv2.front() << '\n';
}

void
test_tournament_tree()
{
    TournamentTree<int> ttree(8, 10000);

    std::vector<int> result;

    size_t winner_id;
    result.push_back(ttree.init_run({ 3, 1, 5, 6, 3, 7, 4, 0 }, winner_id));
    for (int i = 0; i < 8; i++) {
        result.push_back(ttree.replace_run(1000));
    }
    std::copy(result.begin(), result.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << '\n';
}

void
test_bgzf_write()
{
    std::vector<int> vec;
    size_t sz = 100 * 1024 * 1014;
    vec.resize(sz);
    vec[100] = 12345;

    {
        ScopedTimer timer("bgzf_write");
        BGZF *bg_fp = bgzf_open("1.gz", "wu");
        assert(sz * sizeof(int) == bgzf_write(bg_fp, &vec[0], sz * sizeof(int)));
        bgzf_close(bg_fp);
    }

    {
        ScopedTimer timer("bgzf_read");
        int val = 0;
        BGZF *fp = bgzf_open("1.gz", "ru");
        assert(bgzf_seek(fp, sizeof(int) * 100, SEEK_SET) == 0);
        assert(1 * sizeof(int) == bgzf_read(fp, &val, 1 * sizeof(int)));
        bgzf_close(fp);
        std::cout << val << '\n';
    }
}

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

void
test_move_default_move_constructor()
{
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

    std::cout << a2.b->str << '\n';
    char buffer[1000];
    fgets(buffer, 1000, a2.fp);
    std::cout << buffer << '\n';
    fclose(a2.fp);

    std::cout << a1.b->vec[0] << '\n';
}

void
test_ibdsorter()
{
    MetaFile meta;
    meta.parse_files(vcf_fn, map_fn, true);

    IbdFile ibdfile1(temp_file, &meta);
    ibdfile1.open("wu");
    ibdfile1.from_raw_ibd(ibd_in_file, 0, 1, 2, 3);
    ibdfile1.close();

    std::string out_prefix = temp_file2;
    out_prefix += "_temp_";
    IbdSorter sorter(temp_file, temp_file2, "wu", out_prefix.c_str(), 1000);
    {
        ScopedTimer timer("sort_into_chunks");
        sorter.sort_into_chunks();
    }
    {
        ScopedTimer timer("merge_chunks");
        sorter.merge_chunks(3);
    }

    IbdFile file_orig(temp_file, NULL);
    file_orig.open("r");
    IbdFile file_sort(temp_file2, NULL);
    file_sort.open("r");

    {
        ScopedTimer timer("read_from_file");
        file_orig.read_from_file();
        file_sort.read_from_file();
    }

    file_orig.close();
    file_sort.close();

    auto &orig_vec = file_orig.get_vec();
    std::sort(orig_vec.begin(), orig_vec.end());

    std::cout << "orig:  " << orig_vec.size() << '\n';
    file_orig.summary();
    std::cout << "sort:  " << file_sort.get_vec().size() << '\n';
    file_sort.summary();
    std::cout << "orig == sort? " << (orig_vec == file_sort.get_vec()) << '\n';
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
main()
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
    //  test_ibdsorter();
    // test_auto_variable();
    // test_bit_fields();
    test_ibdfile_decode_pakced_ibd();
}
