#include "../include/common.hpp"
#include "../include/ibdcoverage.hpp"
#include "../include/ibdfile.hpp"
#include "../include/ibdmatrix.hpp"
#include "../include/ibdmerger.hpp"
#include "../include/ibdsorter.hpp"
#include "../include/ibdspliter.hpp"
#include "../include/metafile.hpp"
#include "htslib/bgzf.h"
#include <argp.h>
#include <boost/program_options.hpp>
#include <cstdint>
#include <fmt/core.h>
#include <fstream>
#include <limits>
#include <numeric>

using namespace std;
using namespace boost::program_options;

int
ibdtools_encode_main(int argc, char *argv[])
{

    string ibd_in_fn, vcf_in_fn, map_in_fn, chr_name, ibd_out_fn, meta_out_fn;
    float mem = 10.0;

    options_description desc{ "ibdtools encode" };

    try {
        auto add = desc.add_options();

        add("ibd_in,i", value<string>(&ibd_in_fn)->required(), "raw txt ibd file");
        add("vcf_in,v", value<string>(&vcf_in_fn)->required(), "VCF file");
        add("gmap_in,g", value<string>(&map_in_fn)->required(),
            "recombination map file");
        add("chr_name,c", value<string>(&chr_name)->required(), "chromosome name or id");
        add("ibd_out,o", value<string>(&ibd_out_fn)->required(), "output ibd(encoded)");
        add("meta_out,m", value<string>(&meta_out_fn)->required(), "RAM to use (Gb)");
        add("mem,M1", value<float>(&mem)->default_value(10.0), "output metafile");
        add("help,h", "print help information");

        variables_map vm;
        auto options = parse_command_line(argc, argv, desc);
        store(options, vm);

        if (vm.count("help")) {
            cerr << desc << '\n';
            exit(-1);
        }

        notify(vm);

        cerr << "Options received: \n";
        cerr << "--ibd_in: " << ibd_in_fn << '\n';
        cerr << "--vcf_in: " << vcf_in_fn << '\n';
        cerr << "--gmap_in: " << map_in_fn << '\n';
        cerr << "--chr_name: " << chr_name << '\n';
        cerr << "--ibd_out: " << ibd_out_fn << '\n';
        cerr << "--meta_out: " << meta_out_fn << '\n';
        cerr << "--mem: " << mem << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cerr << desc << '\n';
        exit(-1);
    }

    // open file for save meta file
    BGZF *fp = NULL;
    fp = bgzf_open(meta_out_fn.c_str(), "w");
    assert(fp != NULL);
    bgzf_index_build_init(fp);

    MetaFile meta;
    meta.parse_files(vcf_in_fn.c_str(), map_in_fn.c_str(), true, chr_name.c_str());
    meta.write_to_file(fp);
    assert(0 == bgzf_index_dump(fp, meta_out_fn.c_str(), ".gzi"));
    bgzf_close(fp);
    fp = NULL;

    // encode ibdfile
    IbdFile ibd_in(ibd_out_fn.c_str(), &meta, 1 * 1024 * 1024 * 1024);
    ibd_in.open("w");
    ibd_in.from_raw_ibd(ibd_in_fn.c_str());
    ibd_in.close();

    return 0;
}
////////////////////////////////////////////////////////////

int
ibdtools_coverage_main(int argc, char *argv[])
{
    string ibd_in_fn, meta_in_fn, coverage_out_fn;
    float mem = 10.0;   // 10gb
    float window = 1.0; // 1cM

    options_description desc{ "ibdtools coverage" };
    try {
        auto add = desc.add_options();
        add("ibd_in,i", value<string>(&ibd_in_fn)->required(), "encoded ibd file");
        add("meta,m", value<string>(&meta_in_fn)->required(), "meta file");
        add("mem", value<float>(&mem)->default_value(10.0), "RAM to use in Gb");
        add("window", value<float>(&window)->default_value(1.0), "window size in cM");
        add("out,o", value<string>(&coverage_out_fn)->required(),
            "output file for coverage");
        add("help,h", "help information");
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cout << desc << '\n';
            exit(1);
        }

        notify(vm);
        cerr << "Options received: \n";
        cout << "--ibd_in: " << ibd_in_fn << '\n';
        cout << "--meta: " << meta_in_fn << '\n';
        cout << "--out: " << coverage_out_fn << '\n';
        cout << "--mem: " << mem << '\n';
        cout << "--window: " << window << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cout << desc << '\n';
        exit(-1);
    }

    IbdCoverage cov(
        ibd_in_fn.c_str(), meta_in_fn.c_str(), 1.0, mem / 10 * 1024 * 1024 * 1024);
    cov.calculate_coverage();
    ofstream ofs(coverage_out_fn);
    cov.summary(ofs, 3000);
    return 0;
}

////////////////////////////////////////////////////////////
int
ibdtools_split_main(int argc, char *argv[])
{
    string meta_in_fn, ibd_in_fn, out_prefix;
    vector<region_label_t> labels;
    float min_cM = 2.0;
    float mem = 10.0;

    options_description desc{ "ibdtools split" };
    try {
        auto add = desc.add_options();
        add("ibd_in,i", value<string>(&ibd_in_fn)->required(), "input ibd (encoded)");
        add("meta_in,m", value<string>(&meta_in_fn)->required(),
            "input metafile (encoded)");
        add("min_cM", value<float>(&min_cM)->default_value(2.0),
            "mininum length to keep after splitting (cM)");
        add("mem", value<float>(&mem)->default_value(10), "RAM to use");
        add("out_prefix,o", value<string>(&out_prefix)->required(),
            "output file prefix (encoded)");
        add("help,h", "print help information");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cout << desc << '\n';
            exit(1);
        }

        notify(vm);

        cerr << "Options received: \n";
        cout << "--ibd_in: " << ibd_in_fn << '\n';
        cout << "--out: " << out_prefix << '\n';
        cout << "--meta: " << meta_in_fn << '\n';
        cout << "--min_cM: " << min_cM << '\n';
        cout << "--mem: " << mem << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cout << desc << '\n';
        exit(-1);
    }

    {
        BGZF *fp = bgzf_open(meta_in_fn.c_str(), "r");
        MetaFile meta;
        meta.read_from_file(fp);
        bgzf_close(fp);
        labels = meta.get_positions().get_gap_vector(2.0, 100);

        ofstream ofs(out_prefix + "_label.txt");
        ofs << "pid\tpos_bp\tlabel\n";
        for (auto label : labels)
            ofs << label.pid_s << '\t' << meta.get_positions().get_bp(label.pid_s)
                << '\t' << label.label << '\n';
    }

    IbdSplitter splitter(ibd_in_fn.c_str(), out_prefix.c_str(), meta_in_fn.c_str(),
        labels, 1, min_cM, mem / 10 * 0.33 * 1024 * 1024 * 1024,
        mem / 10 * 0.66 * 1024 * 1024 * 1024);
    splitter.split();

    return 0;
}
////////////////////////////////////////////////////////////
int
ibdtools_sort_main(int argc, char *argv[])
{
    string ibd_in, ibd_out;
    float mem;
    size_t k_way;
    options_description desc{ "ibdtools sort" };

    try {
        auto add = desc.add_options();
        add("ibd_in,i", value<string>(&ibd_in)->required(),
            "input IBD file (encoded, unsorted)");
        add("ibd_out,o", value<string>(&ibd_out)->required(),
            "output ibd file (encoded");
        add("mem,M", value<float>(&mem)->default_value(10.0), "RAM to use (Gb)");
        add("k_way,K", value<size_t>(&k_way)->default_value(10),
            "no. way for the merging step");
        add("help,h", "print help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cerr << desc << '\n';
            exit(-2);
        }
        notify(vm);

        cerr << "Options received: \n";
        cerr << "--ibd_in: " << ibd_in << '\n';
        cerr << "--ibd_out: " << ibd_out << '\n';
        cerr << "--mem: " << mem << '\n';
        cerr << "--k_way: " << k_way << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cerr << desc << '\n';
        exit(-1);
    }
    IbdSorter sorter(ibd_in.c_str(), ibd_out.c_str(), "w", (ibd_out + "_tmp_").c_str(),
        mem / 10 * 1024 * 1024 * 1024);
    sorter.sort_into_chunks();
    sorter.merge_chunks(k_way);

    return 0;
}

////////////////////////////////////////////////////////////
int
ibdtools_merge_main(int argc, char *argv[])
{
    string ibd_in, ibd_out, meta_in;
    float mem;
    int max_snp;
    float max_cm;

    options_description desc{ "ibdtools merge" };
    try {
        auto add = desc.add_options();
        add("ibd_in,i", value<string>(&ibd_in)->required(),
            "input IBD file (from ibdtools sort)");
        add("meta_in,m", value<string>(&meta_in)->required(),
            "input meta file (from ibdtools encode)");
        add("ibd_out,o", value<string>(&ibd_out)->required(),
            "output IBD file (encoded)");
        add("max_snp,S", value<int>(&max_snp)->default_value(1),
            "max no. of discordant homozygote SNPs in the gap of two nearby IBD "
            "segments belonging to the sample sample pair");
        add("max_cm,C", value<float>(&max_cm)->default_value(0.6),
            "max distance in cM between two nearby IBD segments belonging the same "
            "sample pair that can be merged");
        add("mem,M", value<float>(&mem)->default_value(10.0), "RAM to use (Gb)");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cerr << desc << '\n';
            exit(-2);
        }
        notify(vm);

        cerr << "Options received: \n";
        cerr << "--ibd_in: " << ibd_in << '\n';
        cerr << "--ibd_out: " << ibd_out << '\n';
        cerr << "--max_snp: " << max_snp << '\n';
        cerr << "--max_cm: " << max_cm << '\n';
        cerr << "--mem: " << mem << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cerr << desc << '\n';
        exit(-1);
    }

    IbdMerger merger(ibd_in.c_str(), ibd_out.c_str(), "w", meta_in.c_str(),
        mem / 10 * 1024 * 1024 * 1024, max_snp, max_cm);
    merger.merge();
    return 0;
}

////////////////////////////////////////////////////////////
int
ibdtools_matrix_main(int argc, char *argv[])
{
    string ibd_in, meta_in, out_prefix;
    float mem;
    int use_hap_pair;
    float filt_lower_cm, filt_upper_cm;
    uint16_t filt_lower_cm10x, filt_upper_cm10x;
    auto filt_max = numeric_limits<uint16_t>::max();

    options_description desc{ "ibdtools matrix" };
    try {
        auto add = desc.add_options();
        add("ibd_in,i", value<string>(&ibd_in)->required(), "input ibdfile (encoded)");
        add("meta_in,m", value<string>(&meta_in)->required(), "input metafile");
        add("out_prefix,o", value<string>(&out_prefix)->required(),
            "output file prefix");
        add("use_hap_pair,H", value<int>(&use_hap_pair)->default_value(0),
            "1: calculate haplotype-pair total IBD; 0: calculate sample-pair total IBD");
        add("mem,M", value<float>(&mem)->default_value(10.0), "RAM to use");
        add("filt_lower_cm,L", value<float>(&filt_lower_cm)->default_value(0),
            "element of total IBD matrix with value less than filt_low is set to 0");
        add("filt_upper_cm,U",
            value<float>(&filt_upper_cm)->default_value(filt_max / 10.0),
            "element of total IBD matrix with value less than filt_low is set to 0");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cerr << desc << '\n';
            exit(-1);
        }

        notify(vm);

        cerr << "Options received: \n";
        cerr << "--ibd_in: " << ibd_in << '\n';
        cerr << "--meta_in: " << meta_in << '\n';
        cerr << "--out_prefix: " << out_prefix << '\n';
        cerr << "--use_hap_pair: " << use_hap_pair << '\n';
        cerr << "--mem: " << mem << '\n';
        cerr << "--filt_lower_cm: " << filt_lower_cm << '\n';
        cerr << "--filt_upper_cm: " << filt_upper_cm << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cerr << desc << '\n';
        exit(-1);
    }

    // prepare meta and ibd file objects
    MetaFile meta;
    BGZF *fp = bgzf_open(meta_in.c_str(), "r");
    assert(fp != NULL);
    meta.read_from_file(fp);
    bgzf_close(fp);
    fp = NULL;
    IbdFile ibdfile(ibd_in.c_str(), &meta, mem / 10 * 1024 * 1024 * 1024);

    // matrix: calculate total, filter, save and get histogram
    IbdMatrix mat;
    ibdfile.open("r");
    mat.calculate_total_from_ibdfile(ibdfile, use_hap_pair != 0);
    ibdfile.close();

    // filter
    filt_lower_cm10x = 10 * filt_lower_cm;
    filt_upper_cm10x = 10 * filt_upper_cm;
    assert(filt_lower_cm10x < filt_upper_cm10x);
    string mat_fn;
    if (filt_lower_cm10x != 0 || filt_upper_cm10x != filt_max) {
        mat.filter_matrix(filt_lower_cm10x, filt_upper_cm10x);
        mat_fn = out_prefix + "_" + to_string(filt_lower_cm) + "_"
                 + to_string(filt_upper_cm) + ".mat";
    } else {
        mat_fn = out_prefix + ".mat";
    }
    // save mat file
    mat.write_matrix_file(mat_fn.c_str());

    // save the histogram
    vector<size_t> hist;
    uint16_t win_10th_cM = 1;
    mat.get_histogram(hist, win_10th_cM);
    size_t total_element = 0, total_nonzeros = 0;
    for (size_t i = 0; i < hist.size(); i++) {
        if (i != 0)
            total_nonzeros += hist[i];
        total_element += hist[i];
    }

    ofstream hist_ofs(out_prefix + "_hist.txt");
    hist_ofs << "# window_size (cM): " << win_10th_cM / 10.0 << '\n';
    hist_ofs << "# total no. of elements: " << total_element << '\n';
    hist_ofs << "# total no. of nonzero elements: " << total_nonzeros << '\n';
    hist_ofs << "Start_cM\tCount\n";
    for (size_t i = 0; i < hist.size(); i++)
        hist_ofs << win_10th_cM * i / 10.0 << '\t' << hist[i] << '\n';

    return 0;
}
////////////////////////////////////////////////////////////
int
ibdtools_decode_main(int argc, char *argv[])
{

    string ibd_in, meta_in, ibd_out;
    float mem = 10.0;

    options_description desc{ "ibdtools decode" };

    try {
        auto add = desc.add_options();

        add("ibd_in,i", value<string>(&ibd_in)->required(), "input ibd file (encoded)");
        add("meta_in,m", value<string>(&meta_in)->required(), "RAM to use (Gb)");
        add("mem,M", value<float>(&mem)->default_value(10.0), "output metafile");
        add("ibd_out,o", value<string>(&ibd_out)->required(), "output ibd file(txt)");
        add("help,h", "print help information");

        variables_map vm;
        auto options = parse_command_line(argc, argv, desc);
        store(options, vm);

        if (vm.count("help")) {
            cerr << desc << '\n';
            exit(-1);
        }

        notify(vm);

        cerr << "Options received: \n";
        cerr << "--ibd_in: " << ibd_in << '\n';
        cerr << "--meta_in: " << meta_in << '\n';
        cerr << "--ibd_out: " << ibd_out << '\n';
        cerr << "--mem: " << mem << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cerr << desc << '\n';
        exit(-1);
    }

    // open file for save meta file
    BGZF *fp = NULL;

    fp = bgzf_open(meta_in.c_str(), "r");
    assert(fp != NULL);
    MetaFile meta;
    meta.read_from_file(fp);
    bgzf_close(fp);
    fp = NULL;

    // encode ibdfile
    IbdFile ibdfile(ibd_in.c_str(), &meta, mem / 10.0 * 0.33 * 1024 * 1024 * 1024);
    ibdfile.open("r");
    ibdfile.to_raw_ibd(ibd_out.c_str(), mem / 10.0 * 0.66 * 1024 * 1024 * 1024);
    ibdfile.close();

    return 0;
}

////////////////////////////////////////////////////////////

int
main(int argc, char *argv[])
{
    vector<string> subcmds
        = { "encode", "coverage", "split", "sort", "merge", "matrix", "decode" };
    // "summary", "view"

    auto err_msg = [&]() {
        cerr << "Usage: ibdtools <subcommand> [options]. Avaiable subcommands: \n";
        for (auto s : subcmds)
            cerr << "\t" << s << '\n';
        cerr << "  For more information, run: ibdtools <subcommand> --help\n";
    };

    if (argc < 2) {

        err_msg();
        return -1;
    }

    string subcmd(argv[1]);
    if (subcmd == "-h" || subcmd == "--help") {
        err_msg();
        return -1;
    } else {
        ScopedTimer timer("run time", true);
        if (subcmd == "encode") {
            return ibdtools_encode_main(argc - 1, argv + 1);
        } else if (subcmd == "coverage") {
            return ibdtools_coverage_main(argc - 1, argv + 1);
        } else if (subcmd == "split") {
            return ibdtools_split_main(argc - 1, argv + 1);
        } else if (subcmd == "sort") {
            return ibdtools_sort_main(argc - 1, argv + 1);
        } else if (subcmd == "merge") {
            return ibdtools_merge_main(argc - 1, argv + 1);
        } else if (subcmd == "matrix") {
            return ibdtools_matrix_main(argc - 1, argv + 1);
        } else if (subcmd == "decode") {
            return ibdtools_decode_main(argc - 1, argv + 1);
        } else {
            cerr << "Subcommand `" << subcmd << "` NOT Implemented!\n";
            exit(-1);
        }
    }
}