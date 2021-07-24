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
#include <sys/types.h>

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
        add("meta_out,m", value<string>(&meta_out_fn)->required(), "output metafile");
        add("mem,M1", value<float>(&mem)->default_value(10.0), "RAM to use (Gb)");
        add("help,h", "print help information");

        variables_map vm;
        auto options = parse_command_line(argc, argv, desc);
        store(options, vm);

        if (vm.count("help")) {
            cerr << desc << '\n';
            exit(-1);
        }

        notify(vm);

        cerr << "ibdtools encode options received: \n";
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
ibdtools_snpdens_main(int argc, char *argv[])
{

    string meta_in_fn, snp_density_fn;
    float window_cM;

    options_description desc{ "ibdtools snpdens" };

    try {
        auto add = desc.add_options();

        add("meta_in,m", value<string>(&meta_in_fn)->required(), "input meta file");
        add("window_cM,W", value<float>(&window_cM)->default_value(2.0),
            "window size in cM");
        add("snp_density_out,o", value<string>(&snp_density_fn)->required(),
            "output snp desnity (txt file)");
        add("help,h", "print help information");

        variables_map vm;
        auto options = parse_command_line(argc, argv, desc);
        store(options, vm);

        if (vm.count("help")) {
            cerr << desc << '\n';
            exit(-1);
        }

        notify(vm);

        cerr << "ibdtools snpdens options received: \n";
        cerr << "--meta_in: " << meta_in_fn << '\n';
        cerr << "--window_cM: " << window_cM << '\n';
        cerr << "--snp_density_out: " << snp_density_fn << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cerr << desc << '\n';
        exit(-1);
    }

    // open file for save meta file
    BGZF *fp = NULL;
    fp = bgzf_open(meta_in_fn.c_str(), "r");
    assert(fp != NULL);

    MetaFile meta;
    meta.read_from_file(fp, false);
    bgzf_close(fp);
    fp = NULL;

    // output snp density file
    ofstream ofs(snp_density_fn);
    auto counts = meta.get_positions().get_window_counts(window_cM);
    ofs << "# Chromosome: "
        << meta.get_chromosomes().get_name(meta.get_positions().get_chrom_id()) << '\n';
    ofs << "# Window size (cM) : " << window_cM << '\n';
    ofs << "# No. of windows: " << counts.size() << '\n';
    ofs << "cM\tCount\n";
    for (size_t i = 0; i < counts.size(); i++)
        ofs << window_cM * i << '\t' << counts[i] << '\n';

    return 0;
}
//
//
//

int
ibdtools_coverage_main(int argc, char *argv[])
{
    string ibd_in_fn, meta_in_fn, coverage_out_fn, subpop_fn;
    float mem = 10.0;   // 10gb
    float window = 1.0; // 1cM

    options_description desc{ "ibdtools coverage" };
    try {
        auto add = desc.add_options();
        add("ibd_in,i", value<string>(&ibd_in_fn)->required(), "encoded ibd file");
        add("meta,m", value<string>(&meta_in_fn)->required(), "meta file");
        add("subpop_samples,P", value<string>(&subpop_fn)->default_value("(None)"),
            "file contains a list subpopulation samples");
        add("mem", value<float>(&mem)->default_value(10.0), "RAM to use in Gb");
        add("window", value<float>(&window)->default_value(1.0), "window size in cM");
        add("out,o", value<string>(&coverage_out_fn)->required(),
            "output file for coverage");
        add("help,h", "help information");
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cerr << desc << '\n';
            exit(1);
        }

        notify(vm);
        cerr << "ibdtools coverage options received: \n";
        cerr << "--ibd_in: " << ibd_in_fn << '\n';
        cerr << "--meta: " << meta_in_fn << '\n';
        cerr << "--subpop_samples: " << subpop_fn << '\n';
        cerr << "--out: " << coverage_out_fn << '\n';
        cerr << "--mem: " << mem << '\n';
        cerr << "--window: " << window << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cerr << desc << '\n';
        exit(-1);
    }

    const char *subpop_fn_char = NULL;
    if (subpop_fn != "(None)")
        subpop_fn_char = subpop_fn.c_str();

    IbdCoverage cov(ibd_in_fn.c_str(), meta_in_fn.c_str(), window,
        mem / 10 * 1024 * 1024 * 1024, subpop_fn_char);
    cov.calculate_coverage();
    ofstream ofs(coverage_out_fn);
    cov.summary(ofs, 0);
    return 0;
}

////////////////////////////////////////////////////////////
int
ibdtools_split_main(int argc, char *argv[])
{
    string meta_in_fn, ibd_in_fn, out_prefix, exclusion_range_fn;
    vector<region_label_t> labels;
    int out_range_label_only = 0;
    float window_cM = 2.0;
    size_t min_snp_in_window = 100;
    float min_cM = 2.0;
    float mem = 10.0;
    struct range_t {
        uint32_t pid1, pid2;
    };

    options_description desc{ "ibdtools split" };
    try {
        auto add = desc.add_options();
        add("ibd_in,i", value<string>(&ibd_in_fn)->required(), "input ibd (encoded)");
        add("meta_in,m", value<string>(&meta_in_fn)->required(),
            "input metafile (encoded)");
        add("window_cM,W", value<float>(&window_cM)->default_value(2.0),
            "window to count SNPs (cM)");
        add("min_snp_in_window,S", value<size_t>(&min_snp_in_window)->default_value(10),
            "min snp per windows for exclusion");
        add("min_cM", value<float>(&min_cM)->default_value(2.0),
            "min length to keep after splitting (cM)");
        add("exclusion_range_fn",
            value<string>(&exclusion_range_fn)->default_value("(None)"),
            "additional range to exclude in addition to low SNP density region");
        add("mem", value<float>(&mem)->default_value(10), "RAM to use");
        add("out_prefix,o", value<string>(&out_prefix)->required(),
            "output file prefix (encoded)");
        add("out_range_label_only,R",
            value<int>(&out_range_label_only)->default_value(0),
            "if nonzero, only output the range_labels for splitting but not run split");
        add("help,h", "print help information");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cerr << desc << '\n';
            exit(1);
        }

        notify(vm);

        cerr << "ibdtools split options received: \n";
        cerr << "--ibd_in: " << ibd_in_fn << '\n';
        cerr << "--meta_in: " << meta_in_fn << '\n';
        cerr << "--window_cM: " << window_cM << '\n';
        cerr << "--min_snp_in_window: " << min_snp_in_window << '\n';
        cerr << "--min_cM: " << min_cM << '\n';
        cerr << "--exclusion_range_fn: " << exclusion_range_fn << '\n';
        cerr << "--mem: " << mem << '\n';
        cerr << "--out_prefix: " << out_prefix << '\n';
        cerr << "--out_range_label_only: " << out_range_label_only << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cerr << desc << '\n';
        exit(-1);
    }

    {
        // ranges by SNP density
        BGZF *fp = bgzf_open(meta_in_fn.c_str(), "r");
        MetaFile meta;
        meta.read_from_file(fp);
        bgzf_close(fp);
        labels = meta.get_positions().get_gap_vector(window_cM, min_snp_in_window);

        // Add ranges to exclude to the labels vector.
        if (exclusion_range_fn != "(None)") {
            ifstream ifs(exclusion_range_fn);
            uint32_t left, right;
            uint32_t pid1, pid2;
            string line;
            StringViewSplitter spliter("\t");
            while (getline(ifs, line, '\n')) {
                cerr << "line" << line << '\n';
                spliter.split(line, 2);
                spliter.get(0, left);
                spliter.get(1, right);
                auto &pos = meta.get_positions();
                pid1 = pos.get_upper_bound_id(left) - 1;
                pid2 = pos.get_upper_bound_id(right);

                add_exclusion_range(labels, pid1, pid2);
            }
        }

        // Output to file
        ofstream ofs(out_prefix + "_label.txt");
        ofs << "pid\tpos_bp\tcM\tlabel\n";
        for (auto label : labels)
            ofs << label.pid_s << '\t' << meta.get_positions().get_bp(label.pid_s)
                << '\t' << meta.get_positions().get_cm(label.pid_s) << '\t'
                << label.label << '\n';
    }

    if (out_range_label_only == 0) {
        IbdSplitter splitter(ibd_in_fn.c_str(), out_prefix.c_str(), meta_in_fn.c_str(),
            labels, 1, min_cM, mem / 10 * 0.33 * 1024 * 1024 * 1024,
            mem / 10 * 0.66 * 1024 * 1024 * 1024);
        splitter.split();
    }

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

        cerr << "ibdtools sort options received: \n";
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

        cerr << "ibdtools merge options received: \n";
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
    string ibd_in, meta_in, out_prefix, subpop_fn;
    vector<string> matrices_in;
    float mem, hist_win_cM;
    int use_hap_pair;
    float filt_lower_cm, filt_upper_cm;
    uint16_t filt_lower_cm10x, filt_upper_cm10x;
    auto filt_max = numeric_limits<uint16_t>::max();

    options_description desc{ "ibdtools matrix" };
    try {
        auto add = desc.add_options();
        add("ibd_in,i", value<string>(&ibd_in), "input ibdfile (encoded)");
        add("matrices_in", value<vector<string> >(&matrices_in)->multitoken(),
            "matrices to be combined. When this is given, --meta_in can be the metafile "
            "of any chromsome as it only uses sample information from it which is "
            "shared among chromosmes");
        add("meta_in,m", value<string>(&meta_in)->required(), "input metafile");
        add("hist_win_cM,W", value<float>(&hist_win_cM)->default_value(1),
            "window size (cM) for histogram");
        add("subpop_samples,P", value<string>(&subpop_fn)->default_value("(None)"),
            "file contains a list subpopulation samples");
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

        cerr << "ibdtools matrix options received: \n";
        cerr << "--ibd_in: " << ibd_in << '\n';
        cerr << "--matrices_in: ";
        for (auto s : matrices_in)
            cerr << s << ' ';
        cerr << '\n';
        cerr << "--meta_in: " << meta_in << '\n';
        cerr << "--hist_win_cM: " << hist_win_cM << '\n';
        cerr << "--subpop_samples: " << subpop_fn << '\n';
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

    // matrix: calculate total, filter, save and get histogram
    IbdMatrix mat;

    if (ibd_in != "") {
        IbdFile ibdfile(ibd_in.c_str(), &meta, mem / 10 * 1024 * 1024 * 1024);
        ibdfile.open("r");

        assert(matrices_in.size() == 0
               && "--ibd_in and --matrices_in are mutual exclusive");
        mat.calculate_total_from_ibdfile(ibdfile, use_hap_pair != 0);

        ibdfile.close();
    } else {
        assert(
            matrices_in.size() > 1 && "--ibd_in and --matrices_in are mutual exclusive");
        mat.read_matrix_file(matrices_in[0].c_str());
        cerr << "read matrix: " << matrices_in[0] << '\n';

        IbdMatrix mat2;
        for (size_t i = 1; i < matrices_in.size(); i++) {
            mat2.read_matrix_file(matrices_in[i].c_str());
            cerr << "read matrix: " << matrices_in[i] << '\n';
            mat.add_matrix(mat2);
            cerr << "add matrix: " << matrices_in[i] << '\n';
        }
    }

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
    const char *subpop_fn_char = NULL;
    if (subpop_fn != "(None)")
        subpop_fn_char = subpop_fn.c_str();

    uint16_t win_10th_cM = round(hist_win_cM * 10);
    mat.get_histogram(hist, win_10th_cM, &meta, subpop_fn_char, use_hap_pair != 0);
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

    string ibd_in, meta_in, ibd_out, subpop_fn;
    float mem = 10.0;

    options_description desc{ "ibdtools decode" };

    try {
        auto add = desc.add_options();

        add("ibd_in,i", value<string>(&ibd_in)->required(), "input ibd file (encoded)");
        add("meta_in,m", value<string>(&meta_in)->required(), "RAM to use (Gb)");
        add("subpop_samples,P", value<string>(&subpop_fn)->default_value("(None)"),
            "file contains a list subpopulation samples");
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

        cerr << "ibdtools decode options received: \n";
        cerr << "--ibd_in: " << ibd_in << '\n';
        cerr << "--meta_in: " << meta_in << '\n';
        cerr << "--subpop_samples: " << subpop_fn << '\n';
        cerr << "--mem: " << mem << '\n';
        cerr << "--ibd_out: " << ibd_out << '\n';

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

    const char *subpop_fn_char = NULL;
    if (subpop_fn != "(None)")
        subpop_fn_char = subpop_fn.c_str();

    // decode ibdfile
    IbdFile ibdfile(ibd_in.c_str(), &meta, mem / 10.0 * 0.33 * 1024 * 1024 * 1024);
    ibdfile.open("r");
    ibdfile.to_raw_ibd(
        ibd_out.c_str(), mem / 10.0 * 0.66 * 1024 * 1024 * 1024, subpop_fn_char);
    ibdfile.close();

    return 0;
}

////////////////////////////////////////////////////////////
int
ibdtools_view_main(int argc, char *argv[])
{
    string ibd_in, meta_in;
    uint32_t sid1, sid2;
    string sample1, sample2;
    float mem;
    options_description desc{ "ibdtools view" };

    try {
        auto add = desc.add_options();
        add("ibd_in,i", value<string>(&ibd_in)->required(),
            "input IBD file (encoded, unsorted)");
        add("meta_in,m", value<string>(&meta_in)->required(),
            "output ibd file (encoded");
        add("sid1,1", value<uint32_t>(&sid1)->default_value(1), "id for sample 1");
        add("sid2,2", value<uint32_t>(&sid2)->default_value(0), "id for sample 2");
        add("sample1", value<string>(&sample1)->default_value(""),
            "name for sample 1. If not empty will override --sid1");
        add("sample2", value<string>(&sample2)->default_value(""),
            "name for sample 2. If not empty will override --sid2");
        add("mem,M", value<float>(&mem)->default_value(10.0), "RAM to use (Gb)");
        add("help,h", "print help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cerr << desc << '\n';
            exit(-2);
        }
        notify(vm);

        cerr << "ibdtools view options received: \n";
        cerr << "--ibd_in: " << ibd_in << '\n';
        cerr << "--meta_in: " << meta_in << '\n';
        if (sample1 != "" || sample2 != "") {
            cerr << "--sample1: " << sample1 << '\n';
            cerr << "--sample2: " << sample2 << '\n';
        } else {
            cerr << "--sid1: " << sid1 << '\n';
            cerr << "--sid2: " << sid2 << '\n';
        }
        cerr << "--mem: " << mem << '\n';

    } catch (const error &ex) {
        cerr << ex.what() << '\n';
        cerr << desc << '\n';
        exit(-1);
    }

    MetaFile meta;
    BGZF *fp = bgzf_open(meta_in.c_str(), "r");
    assert(fp != NULL);
    meta.read_from_file(fp);
    bgzf_close(fp);

    IbdFile ibdfile(ibd_in.c_str(), &meta, mem / 10 * 1024 * 1024 * 1024);

    ibdfile.open("r");
    bool read_full;
    auto &pos = meta.get_positions();
    auto &samples = meta.get_samples();
    auto chrname = meta.get_chromosomes().get_name(pos.get_chrom_id());

    // samepl name overried sid
    if (sample1 != "" || sample2 != "") {
        assert(sample1 != "" && sample2 != "" && sample1 != sample2);
        sid1 = meta.get_samples().get_id(sample1);
        sid2 = meta.get_samples().get_id(sample2);
    }

    if (sid1 < sid2) {
        std::swap(sid1, sid2);
        std::swap(sample1, sample2);
    }

    std::cerr << "sid1: " << sid1 << " sample1: " << meta.get_samples().get_name(sid1)
              << " sid2: " << sid2 << " sample2: " << meta.get_samples().get_name(sid2)
              << '\n';
    do {
        read_full = ibdfile.read_from_file();
        for (auto x : ibdfile.get_vec()) {

            if (x.get_sid1() == sid1 && x.get_sid2() == sid2) {

                float cm1 = pos.get_cm(x.pid1);
                float cm2 = pos.get_cm(x.get_pid2());
                uint32_t start = pos.get_bp(x.pid1);
                uint32_t end = pos.get_bp(x.get_pid2());
                string &s1 = samples.get_name(x.sid1);
                string &s2 = samples.get_name(x.sid2);
                fmt::print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", s1, x.get_hid1(), s2,
                    x.get_hid2(), chrname, start, end, cm2 - cm1);
            }
        }
    } while (read_full);
    ibdfile.close();

    return 0;
}
////////////////////////////////////////////////////////////

int
main(int argc, char *argv[])
{
    vector<string> subcmds = { "encode", "snsdens", "coverage", "split", "sort", "merge",
        "matrix", "decode", "view" };

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
        } else if (subcmd == "snpdens") {
            return ibdtools_snpdens_main(argc - 1, argv + 1);
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
        } else if (subcmd == "view") {
            return ibdtools_view_main(argc - 1, argv + 1);
        } else {
            cerr << "Subcommand `" << subcmd << "` NOT Implemented!\n";
            exit(-1);
        }
    }
}
