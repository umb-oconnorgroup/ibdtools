#include "../include/common.hpp"
#include "../include/ibdcoverage.hpp"
#include "../include/ibdfile.hpp"
#include "../include/ibdmatrix.hpp"
#include "../include/ibdmerger.hpp"
#include "../include/ibdsorter.hpp"
#include "../include/ibdspliter.hpp"
#include "../include/ibdstat.hpp"
#include "../include/metafile.hpp"
#include "cxxopts.hpp"
#include "htslib/bgzf.h"
#include <algorithm>
#include <argp.h>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <fmt/core.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sys/types.h>

using namespace std;
using namespace cxxopts;

static void
check_required_option(Options &options, ParseResult &result, const char *opt)
{
    if (result.count(opt) < 1) {
        cerr << options.help({ "" }) << '\n';
        cerr << "Required Option is not specified: " << opt << '\n';
        exit(-1);
    }
}

int
ibdtools_encode_main(int argc, char *argv[])
{
    string ibd_in_fn, vcf_in_fn, map_in_fn, chr_name, ibd_out_fn, meta_out_fn;
    float mem = 10.0;

    // parse cmd options
    Options options(
        "ibdtools encode", "encode genotype/genetic map information into binary");
    options.set_width(80);

    try {
        auto add = options.add_options();

        add("i,ibd_in", "raw txt ibd file", value<string>(ibd_in_fn));
        add("v,vcf_in", "VCF file", value<string>(vcf_in_fn));
        add("g,gmap_in",
            "recombination map (format: plink .map; separator: single whitespace",
            value<string>(map_in_fn));
        add("c,chr_name", "chromosome name or id", value<string>(chr_name));
        add("o,ibd_out", "output ibd(encoded)", value<string>(ibd_out_fn));
        add("m,meta_out", "output metafile", value<string>(meta_out_fn));
        add("M,mem", "RAM to use (Gb)", value<float>(mem)->default_value("10.0"));
        add("h,help", "print help information");

        auto result = options.parse(argc, argv);
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "vcf_in");
        check_required_option(options, result, "gmap_in");
        check_required_option(options, result, "chr_name");
        check_required_option(options, result, "ibd_out");
        check_required_option(options, result, "meta_out");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }

        cerr << "ibdtools encode options received: \n";
        cerr << "--ibd_in:   " << ibd_in_fn << '\n';
        cerr << "--vcf_in:   " << vcf_in_fn << '\n';
        cerr << "--gmap_in:  " << map_in_fn << '\n';
        cerr << "--chr_name: " << chr_name << '\n';
        cerr << "--ibd_out:  " << ibd_out_fn << '\n';
        cerr << "--meta_out: " << meta_out_fn << '\n';
        cerr << "--mem:      " << mem << '\n';

    } catch (const OptionException &e) {
        cerr << "Error in parsing options: " << e.what() << std::endl;
        exit(-1);
    }

    // call actual functions
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
int
ibdtools_snpdens_main(int argc, char *argv[])
{
    string meta_in_fn, snp_density_fn;
    float window_cM;

    // parse argument
    Options options("ibdtools snpdens", "calculate number snp in given window size");
    options.set_width(80);

    try {
        auto add = options.add_options();

        add("m,meta_in", "input meta file", value<string>(meta_in_fn));
        add("W,window", "window size in cM",
            value<float>(window_cM)->default_value("2.0"));
        add("o,snp_density_out", "output snp desnity (txt file)",
            value<string>(snp_density_fn));
        add("h,help", "print help information");

        auto result = options.parse(argc, argv);
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "snp_density_out");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }

        cerr << "ibdtools snpdens options received: \n";
        cerr << "--meta_in:         " << meta_in_fn << '\n';
        cerr << "--window:       " << window_cM << '\n';
        cerr << "--snp_density_out: " << snp_density_fn << '\n';

    } catch (const OptionException &e) {
        cerr << e.what() << '\n';
        exit(-1);
    }

    // Calculate snp density
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
int
ibdtools_coverage_main(int argc, char *argv[])
{
    string ibd_in_fn, meta_in_fn, coverage_out_fn, subpop_fn;
    float mem = 10.0;   // 10gb
    float window = 1.0; // 1cM

    Options options("ibdtools coverage", "Caculate Snp density using encoded meta file");
    try {
        auto add = options.add_options();
        add("i,ibd_in", "encoded ibd file", value<string>(ibd_in_fn));
        add("m,meta_in", "encoded meta file", value<string>(meta_in_fn));
        add("P,subpop_samples", "file contains a list subpopulation samples",
            value<string>(subpop_fn)->default_value("(None)"));
        add("M,mem", "RAM to use in Gb", value<float>(mem)->default_value("10.0"));
        add("W,window", "window size in cM", value<float>(window)->default_value("1.0"));
        add("o, out", "output file for coverage", value<string>(coverage_out_fn));
        add("h, help", "help information");

        auto result = options.parse(argc, argv);
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "out");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(1);
        }

        cerr << "ibdtools coverage options received: \n";
        cerr << "--ibd_in:         " << ibd_in_fn << '\n';
        cerr << "--meta_in:        " << meta_in_fn << '\n';
        cerr << "--window:         " << window << '\n';
        cerr << "--subpop_samples: " << subpop_fn << '\n';
        cerr << "--out:            " << coverage_out_fn << '\n';
        cerr << "--mem:            " << mem << '\n';

    } catch (const OptionException &ex) {
        cerr << ex.what() << '\n';
        exit(-1);
    }

    // calc coverage
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
int
ibdtools_split_main(int argc, char *argv[])
{

    string meta_in_fn, ibd_in_fn, out_prefix, exclusion_range_fn;
    vector<region_label_t> labels;
    bool out_range_label_only = false;
    float window_cM = 2.0;
    size_t min_snp_in_window = 100;
    float min_cM = 2.0;
    float mem = 10.0;
    struct range_t {
        uint32_t pid1, pid2;
    };

    Options options("ibdtools split", "split and remove IBD within specified "
                                      "region or regions with low snp density");
    try {
        auto add = options.add_options();
        add("i,ibd_in", "input ibd (encoded)", value<string>(ibd_in_fn));
        add("m,meta_in", "input metafile (encoded)", value<string>(meta_in_fn));
        add("W,window", "window to count SNPs (cM)",
            value<float>(window_cM)->default_value("2.0"));
        add("S,min_snp_in_window", "min snp per windows for exclusion",
            value<size_t>(min_snp_in_window)->default_value("10"));
        add("C,min_cM", "min length to keep after splitting (cM)",
            value<float>(min_cM)->default_value("2.0"));
        add("E,exclusion_range_fn",
            "additional range to exclude in addition to low SNP density region",
            value<string>(exclusion_range_fn)->default_value("(None)"));
        add("M,mem", "RAM to use", value<float>(mem)->default_value("10"));

        add("o,out_prefix", "output file prefix (encoded)", value<string>(out_prefix));
        add("O,out_range_label_only",
            "if nonzero, only output the range_labels for splitting but not run "
            "split",
            value<bool>(out_range_label_only)
                ->default_value("false")
                ->implicit_value("True"));
        add("h,help", "print help information");

        auto result = options.parse(argc, argv);
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "out_prefix");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(1);
        }

        cerr << "ibdtools split options received: \n";
        cerr << "--ibd_in:               " << ibd_in_fn << '\n';
        cerr << "--meta_in:              " << meta_in_fn << '\n';
        cerr << "--window:               " << window_cM << '\n';
        cerr << "--min_snp_in_window:    " << min_snp_in_window << '\n';
        cerr << "--min_cM:               " << min_cM << '\n';
        cerr << "--exclusion_range_fn:   " << exclusion_range_fn << '\n';
        cerr << "--mem:                  " << mem << '\n';
        cerr << "--out_prefix:           " << out_prefix << '\n';
        cerr << "--out_range_label_only: " << out_range_label_only << '\n';

    } catch (const OptionException &ex) {
        cerr << ex.what() << '\n';
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
            double leftf, rightf;
            uint32_t pid1, pid2;
            string type;
            string line;
            StringViewSplitter spliter("\t");
            while (getline(ifs, line, '\n')) {
                cerr << "line: " << line << '\n';
                spliter.split(line, 3);
                spliter.get(0, type);
                std::transform(type.begin(), type.end(), type.begin(), ::toupper);

                auto &pos = meta.get_positions();
                auto &gmap = meta.get_genetic_map();

                assert((type == "CM" || type == "BP")
                       && "region list should be either cM or bp");

                if (type == "CM") {
                    spliter.get(1, leftf);
                    spliter.get(2, rightf);
                    uint32_t tmp = pos.get_upper_bound_id(gmap.get_bp(leftf));
                    pid1 = tmp == 0 ? 0 : tmp - 1;
                    pid2 = pos.get_upper_bound_id(gmap.get_bp(rightf));
                } else if (type == "BP") {
                    spliter.get(1, left);
                    spliter.get(2, right);
                    uint32_t tmp = pos.get_upper_bound_id(left);
                    pid1 = tmp == 0 ? 0 : tmp - 1;
                    pid2 = pos.get_upper_bound_id(right);
                    std::cerr << "pid1: " << pid1 << " pid2: " << pid2 << '\n';
                    std::cerr << " " << pos.get_bp(0) << " " << pos.get_cm(0) << '\n';
                }

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
int
ibdtools_sort_main(int argc, char *argv[])
{
    string ibd_in, ibd_out;
    float mem;
    size_t k_way;

    Options options{ "ibdtools sort" };

    try {
        auto add = options.add_options();
        add("i,ibd_in", "input IBD file (encoded, unsorted)", value<string>(ibd_in));
        add("o,ibd_out", "output ibd file (encoded", value<string>(ibd_out));
        add("M,mem", "RAM to use (Gb)", value<float>(mem)->default_value("10.0"));
        add("K,k_way", "no. way for the merging step",
            value<size_t>(k_way)->default_value("10"));
        add("h,help", "print help message");

        auto result = options.parse(argc, argv);
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "ibd_out");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }

        cerr << "ibdtools sort options received: \n";
        cerr << "--ibd_in:  " << ibd_in << '\n';
        cerr << "--ibd_out: " << ibd_out << '\n';
        cerr << "--mem:     " << mem << '\n';
        cerr << "--k_way:   " << k_way << '\n';

    } catch (const OptionException &ex) {
        cerr << ex.what() << '\n';
        exit(-1);
    }
    IbdSorter sorter(ibd_in.c_str(), ibd_out.c_str(), "w", (ibd_out + "_tmp_").c_str(),
        mem / 10 * 1024 * 1024 * 1024);
    sorter.sort_into_chunks();
    sorter.merge_chunks(k_way);

    return 0;
}
int
ibdtools_merge_main(int argc, char *argv[])
{

    string ibd_in, ibd_out, meta_in;
    float mem;
    int max_snp;
    float max_cm;

    Options options("ibdtools merge",
        "merge IBD segments if they close enough and only sperated "
        "by a few discordant sites");
    try {
        auto add = options.add_options();
        add("i,ibd_in", "input IBD file (from ibdtools sort)", value<string>(ibd_in));
        add("m,meta_in", "input meta file (from ibdtools encode)",
            value<string>(meta_in));
        add("o,ibd_out", "output IBD file (encoded)", value<string>(ibd_out));
        add("S,max_snp",
            "max no. of discordant homozygote SNPs in the gap of two nearby IBD "
            "segments belonging to the sample sample pair",
            value<int>(max_snp)->default_value("1"));
        add("C,max_cm",
            "max distance in cM between two nearby IBD segments belonging the same "
            "sample pair that can be merged",
            value<float>(max_cm)->default_value("0.6"));
        add("M,mem", "RAM to use (Gb)", value<float>(mem)->default_value("10.0"));
        add("h,help", "print help message");

        auto result = options.parse(argc, argv);
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "ibd_out");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-2);
        }

        cerr << "ibdtools merge options received: \n";
        cerr << "--ibd_in:  " << ibd_in << '\n';
        cerr << "--meta_in: " << meta_in << '\n';
        cerr << "--ibd_out: " << ibd_out << '\n';
        cerr << "--max_snp: " << max_snp << '\n';
        cerr << "--max_cm:  " << max_cm << '\n';
        cerr << "--mem:     " << mem << '\n';

    } catch (const OptionException &ex) {
        cerr << ex.what() << '\n';
        exit(-1);
    }

    IbdMerger merger(ibd_in.c_str(), ibd_out.c_str(), "w", meta_in.c_str(),
        mem / 10 * 1024 * 1024 * 1024, max_snp, max_cm);
    merger.merge();
    return 0;
}
int
ibdtools_matrix_main(int argc, char *argv[])
{
    string ibd_in, meta_in, out_prefix, subpop_fn;
    vector<string> matrices_in;
    float mem, hist_win_cM;
    int use_hap_pair;
    float min_seg_cM;
    float filt_lower_cm, filt_upper_cm;
    uint16_t filt_lower_cm10x, filt_upper_cm10x;
    auto filt_max = numeric_limits<uint16_t>::max();

    Options options{ "ibdtools matrix" };
    try {
        auto add = options.add_options();
        add("i,ibd_in", "input ibdfile (encoded)", value<string>(ibd_in));
        add("x,matrices_in",
            "matrices to be combined. When this is given, --meta_in can be the "
            "metafile of any chromsome as it only uses sample information from it "
            "which is shared among chromosmes. NOTE: Files should separated with "
            "',' without spaces",
            value<vector<string> >(matrices_in));
        add("m,meta_in", "input metafile", value<string>(meta_in));
        add("T,min_seg_cM",
            "IBD segment not shorter than min_seg_cM will be summed up to total "
            "IBD",
            value<float>(min_seg_cM)->default_value("2"));
        add("W,hist_win_cM", "window size (cM) for histogram",
            value<float>(hist_win_cM)->default_value("1"));
        add("P,subpop_samples", "file contains a list subpopulation samples",
            value<string>(subpop_fn)->default_value("(None)"));
        add("o,out_prefix", "output file prefix", value<string>(out_prefix));
        add("H,use_hap_pair",
            "1: calculate haplotype-pair total IBD; 0: calculate sample-pair total "
            "IBD",
            value<int>(use_hap_pair)->default_value("0"));
        add("M,mem", "RAM to use", value<float>(mem)->default_value("10.0"));
        add("L,filt_lower_cm",
            "element of total IBD matrix with value less than filt_low is set to 0",
            value<float>(filt_lower_cm)->default_value("0"));
        add("U,filt_upper_cm",
            "element of total IBD matrix with value less than filt_low is set to 0",
            value<float>(filt_upper_cm)->default_value(to_string(filt_max / 10.0)));

        auto result = options.parse(argc, argv);
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "matrices_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "out_prefix");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }

        cerr << "ibdtools matrix options received: \n";
        cerr << "--ibd_in: " << ibd_in << '\n';
        cerr << "--matrices_in: ";
        for (auto s : matrices_in)
            cerr << s << ' ';
        cerr << '\n';
        cerr << "--meta_in: " << meta_in << '\n';
        cerr << "--min_seg_cM: " << min_seg_cM << '\n';
        cerr << "--hist_win_cM: " << hist_win_cM << '\n';
        cerr << "--subpop_samples: " << subpop_fn << '\n';
        cerr << "--out_prefix: " << out_prefix << '\n';
        cerr << "--use_hap_pair: " << use_hap_pair << '\n';
        cerr << "--mem: " << mem << '\n';
        cerr << "--filt_lower_cm: " << filt_lower_cm << '\n';
        cerr << "--filt_upper_cm: " << filt_upper_cm << '\n';

    } catch (const OptionException &ex) {
        cerr << ex.what() << '\n';
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
        mat.calculate_total_from_ibdfile(ibdfile, use_hap_pair != 0, NULL, min_seg_cM);

        ibdfile.close();
    } else {
        assert((matrices_in.size() > 1)
               && "--ibd_in and --matrices_in are mutual exclusive");
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
int
ibdtools_decode_main(int argc, char *argv[])
{
    string ibd_in, meta_in, ibd_out, subpop_fn;
    float mem = 10.0;

    Options options("ibdtools decode", "decode binary ibd to text format");

    try {
        auto add = options.add_options();

        add("i,ibd_in", "input ibd file (encoded)", value<string>(ibd_in));
        add("m,meta_in", "RAM to use (Gb)", value<string>(meta_in));
        add("P,subpop_samples", "file contains a list subpopulation samples",
            value<string>(subpop_fn)->default_value("(None)"));
        add("M,mem", "output metafile", value<float>(mem)->default_value("10.0"));
        add("o,ibd_out", "output ibd file(txt)", value<string>(ibd_out));
        add("h,help", "print help information");

        auto result = options.parse(argc, argv);
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "ibd_out");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }

        cerr << "ibdtools decode options received: \n";
        cerr << "--ibd_in:         " << ibd_in << '\n';
        cerr << "--meta_in:        " << meta_in << '\n';
        cerr << "--subpop_samples: " << subpop_fn << '\n';
        cerr << "--mem:            " << mem << '\n';
        cerr << "--ibd_out:        " << ibd_out << '\n';

    } catch (const OptionException &ex) {
        cerr << ex.what() << '\n';
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

int
ibdtools_view_main(int argc, char *argv[])
{
    string ibd_in, meta_in;
    uint32_t sid1, sid2;
    string sample1, sample2;
    float mem;
    Options options("ibdtools view", "view ibd by sample names or ids");

    try {
        auto add = options.add_options();
        add("i,ibd_in", "input IBD file (encoded, unsorted)", value<string>(ibd_in));
        add("m,meta_in", "output ibd file (encoded", value<string>(meta_in));
        add("1,sid1", "id for sample 1", value<uint32_t>(sid1)->default_value("1"));
        add("2,sid2", "id for sample 2", value<uint32_t>(sid2)->default_value("0"));
        add("sample1", "name for sample 1. If not empty will override --sid1",
            value<string>(sample1)->default_value(""));
        add("sample2", "name for sample 2. If not empty will override --sid2",
            value<string>(sample2)->default_value(""));
        add("M,mem", "RAM to use (Gb)", value<float>(mem)->default_value("10.0"));
        add("h,help", "print help message");

        auto result = options.parse(argc, argv);
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }

        cerr << "ibdtools view options received: \n";
        cerr << "--ibd_in:  " << ibd_in << '\n';
        cerr << "--meta_in: " << meta_in << '\n';
        if (sample1 != "" || sample2 != "") {
            cerr << "--sample1: " << sample1 << '\n';
            cerr << "--sample2: " << sample2 << '\n';
        } else {
            cerr << "--sid1:    " << sid1 << '\n';
            cerr << "--sid2:    " << sid2 << '\n';
        }
        cerr << "--mem:     " << mem << '\n';

    } catch (const OptionException &ex) {
        cerr << ex.what() << '\n';
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

int
ibdtools_stat_main(int argc, char *argv[])
{

    string ibd_in, meta_in, stat_out;
    float genome_sz_cM, mem;
    Options options("ibdtools stat", "TODO");

    try {
        auto add = options.add_options();
        add("i,ibd_in", "input IBD file (encoded, unsorted)", value<string>(ibd_in));
        add("m,meta_in", "output ibd file (encoded)", value<string>(meta_in));
        add("o,out", "output file", value<string>(stat_out));
        add("G,genome_size", "genome size in cM",
            value<float>(genome_sz_cM)->default_value("3545.83"));
        add("M,mem", "RAM to use (Gb)", value<float>(mem)->default_value("10.0"));
        add("h,help", "print help message");

        auto result = options.parse(argc, argv);
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "out");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-2);
        }

        assert(stat_out != "");

        cerr << "ibdtools stat options received: \n";
        cerr << "--ibd_in:       " << ibd_in << '\n';
        cerr << "--meta_in:      " << meta_in << '\n';
        cerr << "--genome_sz_cM: " << genome_sz_cM << '\n';
        cerr << "--mem:          " << mem << '\n';
        cerr << "--out:          " << stat_out << '\n';

    } catch (const OptionException &ex) {
        cerr << ex.what() << '\n';
        exit(-1);
    }

    IbdStat stats(ibd_in.c_str(), meta_in.c_str(), genome_sz_cM,
        mem * 1024 * 1024 * 1024 / sizeof(ibd_rec1_t));

    auto counters = stats.get_stat();

    ofstream ofs(stat_out);
    ofs << "# Count of Ibd segments with a length in each window \n";
    ofs << "cM\tCount\n";
    for (size_t i = 0; i < counters.size(); i++) {
        ofs << i << '\t' << counters[i] << '\n';
    }

    return 0;
}

int
main(int argc, char *argv[])
{
    vector<string> subcmds = { "encode", "snpdens", "coverage", "split", "sort", "merge",
        "matrix", "decode", "view", "stat" };
    auto err_msg = [&]() {
        cerr << "Usage: ibdtools <subcmd> [options]. Available subcmd: \n";
        for (auto s : subcmds)
            cerr << "\t - " << s << '\n';
        cerr << " For more information, run: ibdtools <subcmd> --help\n";
    };
    if (argc < 2) {
        err_msg();
        return -1;
    } else {
        string subcmd = argv[1];

        if (subcmd == "-h" || subcmd == "--help") {
            err_msg();
            return -1;
        }
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
        } else if (subcmd == "stat") {
            return ibdtools_stat_main(argc - 1, argv + 1);
        } else {
            cerr << "Subcommand `" << subcmd << "` NOT Implemented!\n";
            exit(-1);
        }
    }
}
