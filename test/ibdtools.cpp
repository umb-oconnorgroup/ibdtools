#include "../include/common.hpp"
#include "../include/ibdcoverage.hpp"
#include "../include/ibdfile.hpp"
#include "../include/ibdspliter.hpp"
#include "../include/metafile.hpp"
#include "htslib/bgzf.h"
#include <argp.h>
#include <boost/program_options.hpp>
#include <fmt/core.h>
#include <fstream>

using namespace std;

int
ibdtools_encode_main(int argc, char *argv[])
{
    using namespace boost::program_options;

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
    using namespace boost::program_options;
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
    using namespace boost::program_options;
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
main(int argc, char *argv[])
{
    vector<string> subcmds = { "encode", "coverage", "split" };

    auto err_msg = [&]() {
        cerr << "ibdtools <subcommand> [options]\n  Avaiable subcommands: \n";
        for (auto s : subcmds)
            cerr << "\t" << s << '\n';
        cerr << "  More help information is available by running:\n\t"
                "ibdtools <subcommand> --help\n";
    };

    if (argc < 2) {

        err_msg();
        return -1;

    } else {

        string subcmd(argv[1]);
        if (subcmd == "-h" || subcmd == "--help") {
            err_msg();
            return -1;
        } else if (subcmd == "encode") {
            return ibdtools_encode_main(argc - 1, argv + 1);
        } else if (subcmd == "coverage") {
            return ibdtools_coverage_main(argc - 1, argv + 1);
        } else if (subcmd == "split") {
            return ibdtools_split_main(argc - 1, argv + 1);
        } else {
            cerr << "Subcommand `" << subcmd << "` NOT Implemented!\n";
            exit(-1);
        }
    }
}
