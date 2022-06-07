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

    Options options("ibdtools encode",
        "`ibdtools encode` encodes the ibd file from IBD caller, the vcf file (phased, "
        "gimputed and biallelic) and the genetic map (plink format) information into "
        "binary. NOTE: Each set of these files should contain information for a single "
        "chromosome. Files should be split if they contain more than one chromosomes. "
        "For vcf file, only the SNP positions and GT field are used. SNP positions by "
        "base-pair should be sorted numerically and unique. SNP are indexed from 0 and "
        "the indices are mapped to position in bp and cM using the plink map. For IBD, "
        "6 columns are encoded, sample1, hap1, sample2, hap2, start_position_id and "
        "end_position_id. Note the start/end poistion are encoded as SNP indices "
        "instead of bp positions to save memory ; additional trick like struct packing "
        "is used to save extra memory. The encoded ibd and meta files (as well as the "
        "*.gzi index files ) are used as input for other subcommands. The encoded IBD "
        "can be decoded into text format via `ibdtools decode` subcommand. Ex:\n\n  "
        "ibdtools encode -i 1.ibd.gz -v 1.vcf.gz -g 1.map -c 1 -o 1.eibd -m 1.meta -M "
        "1\n");
    // options.set_width(80);

    try {
        auto add = options.add_options();

        add("i,ibd_in", "raw txt ibd file (plain text or gzip compressed)",
            value<string>(ibd_in_fn));
        add("v,vcf_in", "VCF file (formats that can bed read by bcftools)",
            value<string>(vcf_in_fn));
        add("g,gmap_in",
            "recombination map (format: plink .map; separator: single whitespace",
            value<string>(map_in_fn));
        add("c,chr_name", "chromosome name or id", value<string>(chr_name));
        add("o,ibd_out", "output ibd(encoded)", value<string>(ibd_out_fn));
        add("m,meta_out", "output metafile", value<string>(meta_out_fn));
        add("M,mem", "RAM to use (Gb)", value<float>(mem)->default_value("10.0"));
        add("h,help", "print help information");

        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "vcf_in");
        check_required_option(options, result, "gmap_in");
        check_required_option(options, result, "chr_name");
        check_required_option(options, result, "ibd_out");
        check_required_option(options, result, "meta_out");

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
    Options options("ibdtools snpdens",
        "`ibdtools snpdens` calculates the number of snps in each window which can be."
        "This step does not need the enocded ibd file. Ex:\n\n  ibdtools snpdens -m "
        "1.meta -o 1_snpdens.txt\n");
    // options.set_width(80);

    try {
        auto add = options.add_options();

        add("m,meta_in", "input meta file", value<string>(meta_in_fn));
        add("W,window", "window size in cM",
            value<float>(window_cM)->default_value("2.0"));
        add("o,snp_density_out", "output snp desnity (txt file)",
            value<string>(snp_density_fn));
        add("h,help", "print help information");

        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "snp_density_out");

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

    Options options("ibdtools coverage",
        "`ibdtools coverage` caculates the numbers of IBD segments overlapping "
        "specified windows. When `-P` is used, the calculation with will be performed "
        "on the subset of IBD that are shared between a pair of samples from the "
        "specified list of samples. Ex:\n\n  ibdtools coverage -i 1.mibd -m 1.meta -M 1 "
        "-o 1_cov.txt\n");
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
        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(1);
        }
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "out");

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

    Options options("ibdtools split",
        "`ibdtools split` splits and removes IBD within specified "
        "specified regions or calculated regions with low snp density. `-W`, `-S` and "
        "`-C` options can be used to specified how low-SNP-density regions will be "
        "calculated. `-E` can be used to exclude additional regions, such as regions "
        "under strong selection. NOTE: For splitting IBD segment, all segments "
        "overapping with the boundaries of target regions are first cut into parts, "
        "those within the target regions and those outside the target regions. The "
        "parts within of the target regions are removed; the parts outside will be "
        "keeped if it is long enough (`C` option). Ex:\n\n   ibdtools split -i 1.eibd "
        "-m 1.meta -W 2 -S 10 -o 1.split -M 1\n");
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

        add("o,out_prefix",
            "output file prefix (encoded), `[out_prefix]1` file contains remaining IBD; "
            "`[out_prefix]0` file contains IBD removed",
            value<string>(out_prefix));
        add("O,out_range_label_only",
            "if nonzero, only output the range_labels for splitting but not run "
            "split",
            value<bool>(out_range_label_only)
                ->default_value("false")
                ->implicit_value("True"));
        add("h,help", "print help information");

        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(1);
        }
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "out_prefix");

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

    Options options("ibdtools sort",
        "`ibdtools sort` sorts the encoded IBD files according IBD indcies of "
        "sample1, hap1, sample2 and hap2, start. The sorted, encoded IBD files are used "
        "in the following sub commands (ibdtools merge, and ibdtools view). Ex:\n\n   "
        "ibdtools sort -i 1.split1 -o 1.sibd -M 1\n");

    try {
        auto add = options.add_options();
        add("i,ibd_in", "input IBD file (encoded, unsorted)", value<string>(ibd_in));
        add("o,ibd_out", "output ibd file (encoded", value<string>(ibd_out));
        add("M,mem", "RAM to use (Gb)", value<float>(mem)->default_value("10.0"));
        add("K,k_way", "no. way for the merging step",
            value<size_t>(k_way)->default_value("10"));
        add("h,help", "print help message");

        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "ibd_out");

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
        "`ibdtools merge` merges IBD segments if they close enough (by comparing IBD "
        "coordinates) and only separated by a few discordant sites (by checking the "
        "phased genotypes).  It can be used to reduce breaks in detected IBD due to "
        "phase errors and genotype errors. The subcommand is to a reimplementation of "
        "the Dr. Browning's tool (merge-ibd-segments.17Jan20.102.jar) for better memory "
        "control and possibly better performance. Ex:\n\n   ibdtools merge -i 1.sibd -m "
        "1.meta -o 1.mibd -M 1\n");
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
        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-2);
        }
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "ibd_out");

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
    float min_seg_cM, max_seg_cM;
    float filt_lower_cm, filt_upper_cm;
    uint16_t filt_lower_cm10x, filt_upper_cm10x;
    auto filt_max = numeric_limits<uint16_t>::max();

    Options options("ibdtools matrix",
        "`ibdtools matrix` aggregates IBD segments into sample-pair total IBD matrix. "
        "There are two modes: i) calculate chromosome-wide sample pair total IBD where "
        "input is the encoded ibd file for the corresponding chromosome (options `-i`); "
        "ii) calculate genome-wide sample pair total IBD where inputs are a "
        "comma-separated list of matrix files generated by running `ibdtools matrix` "
        "mode 1 for each chromosome (option `-x`). Additional, this subcommand allows "
        "three levels of filtering: i) filtering IBD by sample population; ii) "
        "filtering at IBD segment level (option `-T` and `-B`) and iii) filtering at "
        "total IBD level (options `-L` abd `-U`). \nEx1:\n   ibdtools matrix -i 1.mibd "
        "-m 1.meta -M 1 -B 2.5 -L 5 -o 1.mat\n\nEx2:\n   ibdtools matrix -x "
        "1.mat,2.mat,3.mat -m 1.meta -B 2.5 -L 5 -o gw.mat -M 3\n");
    try {
        auto add = options.add_options();
        add("i,ibd_in", "input ibdfile (encoded)",
            value<string>(ibd_in)->default_value(""));
        add("x,matrices_in",
            "matrices to be combined. When this is given, --meta_in can be the "
            "metafile of any chromsome as it only uses sample information from it "
            "which is shared among chromosmes. NOTE: Files should separated with "
            "',' without spaces",
            value<vector<string> >(matrices_in));
        add("m,meta_in", "input metafile", value<string>(meta_in));
        add("B,min_seg_cM",
            "IBD segment not shorter than min_seg_cM will be summed up to total "
            "IBD",
            value<float>(min_seg_cM)->default_value("2"));
        add("T,max_seg_cM",
            "IBD segment shorter than max_seg_cM will be summed up to total "
            "IBD. Any value < min_seg_cM is ignored. ",
            value<float>(max_seg_cM)->default_value("-1"));
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
        if (result.count("ibd_in") < 1 && result.count("matrices_in") < 1) {
            cerr << options.help({ "" }) << '\n';
            cerr << "At least one of --ibd_in and --matrices_in should be specified!"
                 << '\n';
            exit(-1);
        }
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "out_prefix");

        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }

        if (max_seg_cM < min_seg_cM) {
            max_seg_cM = std::numeric_limits<float>::max();
        }

        cerr << "ibdtools matrix options received: \n";
        cerr << "--ibd_in: " << ibd_in << '\n';
        cerr << "--matrices_in: ";
        for (auto s : matrices_in)
            cerr << s << ' ';
        cerr << '\n';
        cerr << "--meta_in: " << meta_in << '\n';
        cerr << "--min_seg_cM: " << min_seg_cM << '\n';
        cerr << "--max_seg_cM: " << max_seg_cM << '\n';
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
        mat.calculate_total_from_ibdfile(
            ibdfile, use_hap_pair != 0, NULL, min_seg_cM, max_seg_cM);

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

    Options options("ibdtools decode",
        "`ibdtools decode` decode binary ibd to compressed text "
        "format. It also allowed filtering IBD by samples. Ex: \n\n  ibdtools decode -i "
        "$i.split1 -m 1.meta -o 1.ibd.gz -M 1\n");

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
        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "ibd_out");

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
    Options options("ibdtools view",
        "view ibd segments shared by a pair of sample names or a pair "
        "ids. NOTE, the order of the sample names and ids are "
        "unimportant as they will bed sorted internally by ibdtools. Ex:\n\n   ibdtools "
        "view -i 1.mibd -m 1.meta -1 4 -2 6\n");

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
        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-1);
        }
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");

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
    Options options("ibdtools stat",
        "`ibdtools stat` calculate the distribution of IBD length by "
        "counting IBD segments in all length bins. The ith bin "
        "correponds to IBD with length in a range of [i, i+1) cM. Ex:\n\n   ibdtools "
        "stat -i 1.mibd -m 1.meta -o 1_stat.txt -M 1\n");

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
        if (result.count("help")) {
            cerr << options.help({ "" }) << '\n';
            exit(-2);
        }
        check_required_option(options, result, "ibd_in");
        check_required_option(options, result, "meta_in");
        check_required_option(options, result, "out");

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
