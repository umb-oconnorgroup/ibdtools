
#include <charconv>
#include <fmt/core.h>

#include "chromosomes.hpp"
#include "common.hpp"
#include "ibdfile.hpp"
#include "metafile.hpp"
#include "positions.hpp"
#include "samples.hpp"
#include <filesystem>

// IbdFile object represent an packed IBD file
IbdFile::IbdFile(const char *ibd_fn, MetaFile *meta_, size_t max_rec)
    : filename{ ibd_fn }, fp{ NULL }, meta(meta_)
{
    ibd_vec.reserve(max_rec);
}

std::vector<ibd_rec1_t> &
IbdFile::get_vec()
{
    return ibd_vec;
}

const std::string &
IbdFile::get_filename()
{
    return filename;
}

BGZF *
IbdFile::get_fp()
{
    return fp;
}

MetaFile *
IbdFile::get_meta()
{
    return meta;
}

void
IbdFile::open(const char *mode)
{
    fp = bgzf_open(filename.c_str(), mode);

    my_assert(fp != NULL, "bgzf_open failed ");
    if (!(fp->is_write)) {
        // std::cout << "ibd_file open file for reading \n";
        my_assert(fp != NULL, "bgzf_open failed to open ibd out file");
        bgzf_mt(fp, 10, 256);
        std::string gzi(filename);
        gzi += ".gzi";
        if (std::filesystem::exists(gzi) && fp->is_compressed) {
            my_assert(bgzf_index_load(fp, filename.c_str(), ".gzi") == 0,
                "bgzf_index_load error");
        }
    } else {
        my_assert(fp != NULL, "bgzf_open failed to open ibd out file");
        bgzf_mt(fp, 10, 256);
        my_assert(bgzf_index_build_init(fp) == 0, " ");
    }
}

void
IbdFile::close()
{
    my_assert(fp != NULL, "bgzf_open failed ");
    if (!(fp->is_write)) {
        bgzf_close(fp);
        // std::cout << "ibdfile closed for reading\n";
    } else {
        // IbdFile object represent an packed IBD file
        if (fp->is_compressed && filename != "/dev/null")
            my_assert(bgzf_index_dump(fp, filename.c_str(), ".gzi") == 0,
                "bgzf_index_dump error");
        bgzf_close(fp);
        // std::cout << "ibdfile closed for reading\n";
    }
    fp = NULL;
}

void
IbdFile::encode_raw_ibd(const char *raw_ibd_in, int col_sample1, int col_sample2,
    int col_start, int col_end, int col_hap1, int col_hap2)
{
    my_assert(fp != NULL, "bgzf_open failed ");
    my_assert(meta != NULL, "encode need meta");
    // check args
    my_assert(col_sample1 < col_sample2 && col_start < col_end, "");

    // references
    auto &samples = meta->get_samples();
    auto &positions = meta->get_positions();

    // buffer string for reading
    kstring_t kstr = { 0 };

    // open file
    BGZF *fp_in = bgzf_open(raw_ibd_in, "r");
    my_assert(fp_in != NULL, "bgzf_open: failed to open ibd in file");

    // threading
    bgzf_mt(fp_in, 10, 256);

    // read and process and raw
    int res;
    StringViewSplitter svs("\t");
    std::string strval;
    uint32_t u32val;
    ibd_rec2_t rec2;

    while ((res = bgzf_getline(fp_in, '\n', &kstr)) > 0) {

        // flush out encoded IBD to disk and empty out memory (vector) for new raw ibd
        if (ibd_vec.size() == ibd_vec.capacity()) {
            write_to_file();
            ibd_vec.resize(0);
        }

        std::string_view sv(kstr.s, kstr.l);

        svs.split(sv);

        svs.get(col_sample1, strval);
        rec2.sid1 = samples.get_id(strval);

        svs.get(col_sample2, strval);
        rec2.sid2 = samples.get_id(strval);

        svs.get(col_start, u32val);

        rec2.pid1 = positions.get_id(u32val);
        svs.get(col_end, u32val);

        rec2.pid2 = positions.get_id(u32val);
        svs.get(col_hap1, u32val);

        // add bit operator to avoid error when reading brownning's merged file
        rec2.hid1 = ((u32val - 1) & 0x1); // raw is 1, 2 for hap; encoded is 0, 1
        svs.get(col_hap2, u32val);
        rec2.hid2 = ((u32val - 1) & 0x1); // raw is 1, 2 for hap; encoded is 0, 1

        if (rec2.sid1 < rec2.sid2) {
            std::swap(rec2.sid1, rec2.sid2);
            std::swap(rec2.hid1, rec2.hid2);
        }

        // push_back an encoded record to the vector
        ibd_vec.push_back(rec2);
    }

    // flush out ibd_vec when all input are parsed but the last chunk still
    // stored in memory
    if (ibd_vec.size() > 0) {
        write_to_file();
    }

    // close input file
    // out file will be closed and index dumped in the IbdFile::close() method.
    bgzf_close(fp_in);
    // free string buffer
    ks_free(&kstr);
}

void
IbdFile::to_raw_ibd(
    const char *raw_ibd_fn, size_t line_buffer_capcity, const char *subpop_fn)
{
    my_assert(fp != NULL, "");
    my_assert(meta != NULL, "decode need meta");
    // references
    auto samples = meta->get_samples();
    auto positions = meta->get_positions();
    auto chromosomes = meta->get_chromosomes();

    // buffer string for writing
    std::string line_buffer;
    line_buffer.reserve(line_buffer_capcity);

    // open file
    BGZF *fp_out = bgzf_open(raw_ibd_fn, "w");
    my_assert(fp_out != NULL, "bgzf_open failed to open ibd out file");

    my_assert(bgzf_index_build_init(fp_out) == 0, "");

    // threading
    bgzf_mt(fp_out, 10, 256);

    // subpoplation samples_ids

    std::vector<uint8_t> subpop_v;
    if (subpop_fn != NULL)
        samples.get_subpop_vector(subpop_fn, subpop_v);

    // loop read in from encoded file and write to decoded file
    bool did_vec_read_full;
    // size_t sid1, sid2, pid1, pid2;
    ibd_rec1_t tmp_rec;
    std::string chrom_name = positions.get_chrom_id() == -1
                                 ? "0"
                                 : chromosomes.get_name(positions.get_chrom_id());
    // debug
    // std::cerr << "chrom_name: " << chrom_name
    //           << "chrom_id: " << positions.get_chrom_id()
    //           << " gen_name: " << chromosomes.get_name(positions.get_chrom_id())
    //           << '\n';

    // The outter loop corresponds to each bulk read
    do {
        did_vec_read_full = read_from_file();

        // The inner loop corresponds to write out each record
        for (auto &rec : ibd_vec) {
            //  write to out file when line_buffer is almost full
            if (line_buffer.size() > line_buffer.capacity() - 1000) {
                auto ret = bgzf_write(fp_out, line_buffer.c_str(), line_buffer.size());
                auto cond = ret >= 0 && ((size_t) ret) == line_buffer.size();
                my_assert(
                    cond, "bgzf_write error during writting ibd_rec during decoding");

                // debug
                // std::cout << line_buffer << '\n';
                //
                line_buffer.clear();
            } else {

                uint32_t sid1 = rec.get_sid1();
                uint32_t sid2 = rec.get_sid2();

                // if either of the two samples is not in the subpop then skip
                if (subpop_v.size() > 0
                    && (subpop_v[sid1] == 0 || subpop_v[sid2] == 0)) {
                    continue;
                }

                // add the string version of each rec to line buffer
                fmt::format_to(std::back_inserter(line_buffer),
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", samples.get_name(sid1),
                    rec.get_hid1() + 1, samples.get_name(sid2), rec.get_hid2() + 1,
                    chrom_name, positions.get_bp(rec.get_pid1()),
                    positions.get_bp(rec.get_pid2()),
                    positions.get_cm(rec.get_pid2()) - positions.get_cm(rec.get_pid1()));

                //	std::cout << "+" << line_buffer.size() << '\n';
            }
        }

    } while (did_vec_read_full);

    // flush the line_buffer
    if (line_buffer.size() > 0) {
        auto ret = bgzf_write(fp_out, line_buffer.c_str(), line_buffer.size());
        auto cond = ret >= 0 && ((size_t) ret) == line_buffer.size();
        const char *msg = "bgzf_write error during writting ibd_rec during decoding";
        my_assert(cond, msg);
        line_buffer.clear();
    }

    // close file and dump index
    my_assert(bgzf_index_dump(fp_out, raw_ibd_fn, ".gzi") == 0, "");
    bgzf_close(fp_out);
}

void
IbdFile::write_to_file(size_t max)
{
    size_t total_bytes = 0;
    ssize_t ret = 0;
    bool cond = false;

    my_assert(fp != NULL, "");

    if (max == 0 || max > ibd_vec.size())
        total_bytes = ibd_vec.size() * sizeof(decltype(ibd_vec)::value_type);
    else
        total_bytes = max * sizeof(decltype(ibd_vec)::value_type);

    // std::cerr << ibd_vec.size() << filename
    //           << "  <---- write_to_file, Total bytes: " << total_bytes << '\n'
    //           << " First vector write: \n";
    // ibd_vec[0].print();
    // ibd_vec[1].print();
    // ibd_vec[2].print();

    ret = bgzf_write(fp, &ibd_vec[0], total_bytes);
    cond = ret > 0 && ((size_t) ret) == total_bytes;
    my_assert(cond, "");
}

void
IbdFile::write_to_file(IbdFile &other, size_t max)
{
    size_t total_bytes;
    ssize_t ret = 0;
    bool cond = false;

    my_assert(other.get_fp() != NULL, "");
    if (max == 0 || max > ibd_vec.size())
        total_bytes = ibd_vec.size() * sizeof(decltype(ibd_vec)::value_type);
    else
        total_bytes = max * sizeof(decltype(ibd_vec)::value_type);

    // std::cout << "write to other file: total bytes " << total_bytes << '\n';
    ret = bgzf_write(other.fp, &ibd_vec[0], total_bytes);
    cond = ret >= 0 && ((size_t) ret) == total_bytes;
    my_assert(cond, "");
}

bool
IbdFile::read_from_file(bool append, size_t new_capacity)
{
    my_assert(fp != NULL, "");

    if (!append)
        ibd_vec.clear();

    if (new_capacity > ibd_vec.capacity())
        ibd_vec.reserve(new_capacity);

    constexpr size_t element_size = sizeof(decltype(ibd_vec)::value_type);

    // exisitng number of elements, and remaining capcity
    size_t n_element = ibd_vec.size();

    // maximize the size to capcity
    ibd_vec.resize(ibd_vec.capacity());

    // calcute cyte capacity
    size_t bytes_capacity = (ibd_vec.capacity() - n_element) * element_size;

    ssize_t bytes_read = bgzf_read(fp, &ibd_vec[n_element], bytes_capacity);
    // verify(bytes_read >= 0 && "bgzf_read error");
    my_assert(bytes_read >= 0, "bgzf_read error");

    // shrink if fewer bytes are read
    if ((size_t) bytes_read < bytes_capacity) {
        ibd_vec.resize(n_element + bytes_read / element_size);

        // std::cout << ibd_vec.size() << "<---  read ibd_vec.size() \n";
        // ibd_vec[0].print();
        // ibd_vec[1].print();
        // ibd_vec[2].print();

        // make sure this happens only at end of the file
        char c;
        my_assert(bgzf_read(fp, &c, 1) == 0, "bgzf_read read error");
        return false;
    }

    return true;
}

void
IbdFile::delete_file_from_disk()
{
    std::filesystem::remove(filename);
}

bool
IbdFile::has_equal_vector(const IbdFile &other)
{
    /*
    std::cout << (ibd_vec == other.ibd_vec) << (meta.is_equal(other.meta))
              << "<- ibdfile is equal \n";
    std::cout << ibd_vec.size() << " size.vs " << other.ibd_vec.size() << '\n';
    */

    return ibd_vec == other.ibd_vec;
}

void
IbdFile::head(size_t n)
{
    for (size_t i = 0; i < n && i < ibd_vec.size(); i++) {
        auto &rec = ibd_vec[i];
        std::cout << rec.get_sid1() << '\t' << rec.get_sid2() << '\t' << rec.get_pid1()
                  << '\t' << rec.get_pid2() << '\n';
    }
}

void
IbdFile::tail(size_t n)
{
    if (n > ibd_vec.size())
        n = ibd_vec.size();
    for (size_t i = ibd_vec.size() - n; i < ibd_vec.size(); i++) {
        auto &rec = ibd_vec[i];
        std::cout << rec.get_sid1() << '\t' << rec.get_sid2() << '\t' << rec.get_pid1()
                  << '\t' << rec.get_pid2() << '\n';
    }
}

void
IbdFile::summary()
{
    std::cout << "File: " << filename << " is_open: " << (fp != NULL)
              << " is_write: " << ((fp == NULL) ? " " : std::to_string(fp->is_write))
              << " ibd_vec size: " << ibd_vec.size() << '\n';
    std::cout << "Head and tail 5 rows: \n";
    head(5);
    std::cout << "... \n";
    tail(5);
}

void
IbdFile::print_to_file(std::ofstream &ofs)
{
    for (auto rec : ibd_vec) {
        ofs << rec.get_sid1() << '\t' << rec.get_sid2() << '\t' << rec.get_pid1() << '\t'
            << rec.get_pid2() << '\n';
    }
}
