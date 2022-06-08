#ifndef __ibdfile_hpp__
#define __ibdfile_hpp__

#include <charconv>

#include "common.hpp"
#include "metafile.hpp"

// IbdFile object represent an packed IBD file
class IbdFile
{
    std::string filename;
    BGZF *fp;

    // buffer
    std::vector<ibd_rec1_t> ibd_vec;

    // pointer
    MetaFile *meta;

  public:
    IbdFile(){};
    IbdFile(const char *ibd_fn, MetaFile *meta_ = NULL, size_t max_rec = 1 * 1024 * 1024)
        : filename{ ibd_fn }, fp{ NULL }, meta(meta_)
    {
        ibd_vec.reserve(max_rec);
    }

    std::vector<ibd_rec1_t> &
    get_vec()
    {
        return ibd_vec;
    }

    const std::string &
    get_filename()
    {
        return filename;
    }

    BGZF *
    get_fp()
    {
        return fp;
    }

    MetaFile *
    get_meta()
    {
        return meta;
    }

    void
    open(const char *mode = "r")
    {
        fp = bgzf_open(filename.c_str(), mode);

        verify(fp != NULL && "bgzf_open failed ");
        if (!(fp->is_write)) {
            // std::cout << "ibd_file open file for reading \n";
            verify(fp != NULL && "bgzf_open failed to open ibd out file");
            bgzf_mt(fp, 10, 256);
            std::string gzi(filename);
            gzi += ".gzi";
            if (std::filesystem::exists(gzi) && fp->is_compressed) {
                verify(bgzf_index_load(fp, filename.c_str(), ".gzi") == 0
                       && "bgzf_index_load error");
            }
        } else {
            verify(fp != NULL && "bgzf_open failed to open ibd out file");
            bgzf_mt(fp, 10, 256);
            verify(bgzf_index_build_init(fp) == 0);
        }
    }

    void
    close()
    {
        verify(fp != NULL);
        if (!(fp->is_write)) {
            bgzf_close(fp);
            // std::cout << "ibdfile closed for reading\n";
        } else {
            if (fp->is_compressed && filename != "/dev/null")
                verify(bgzf_index_dump(fp, filename.c_str(), ".gzi") == 0);
            bgzf_close(fp);
            // std::cout << "ibdfile closed for reading\n";
        }
        fp = NULL;
    }

    void
    from_raw_ibd(const char *raw_ibd_in, int col_sample1 = 0, int col_sample2 = 2,
        int col_start = 5, int col_end = 6, int col_hap1 = 1, int col_hap2 = 3)
    {
        verify(fp != NULL);
        verify(meta != NULL && "encode need meta");
        // check args
        verify(col_sample1 < col_sample2 < col_start < col_end);

        // references
        auto samples = meta->get_samples();
        auto positions = meta->get_positions();

        // buffer string for reading
        kstring_t kstr = { 0 };

        // open file
        BGZF *fp_in = bgzf_open(raw_ibd_in, "r");
        verify(fp_in != NULL && "bgzf_open failed to open ibd in file");

        // threading
        bgzf_mt(fp_in, 10, 256);

        // read and process and raw
        int res;
        StringViewSplitter svs("\t");
        std::string strval;
        uint32_t u32val;
        ibd_rec2_t rec2;

        while ((res = bgzf_getline(fp_in, '\n', &kstr)) > 0) {
            if (ibd_vec.size() == ibd_vec.capacity()) {
                //    std::cerr << "  Write encoded file " << ibd_vec.size() << '\n';
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

            // push_back to vector
            ibd_vec.push_back(rec2);

            // std::cout << "------------\nrec2: ";
            // rec2.print();
            // std::cout << "rec1: ";
            // ibd_vec.back().print();
        }

        // flush out ibd_vec
        if (ibd_vec.size() > 0) {
            // std::cerr << "  Write encoded file " << ibd_vec.size() << '\n';
            write_to_file();

            // For debug purpose last flush not resize the vector
            // ibd_vec.resize(0);
        }

        // close file and dump index
        bgzf_close(fp_in);
        // free string buffer
        ks_free(&kstr);
    }

    void
    to_raw_ibd(const char *raw_ibd_fn, size_t line_buffer_capcity = 10 * 1024 * 1024,
        const char *subpop_fn = NULL)
    {
        verify(fp != NULL);
        verify(meta != NULL && "decode need meta");
        // references
        auto samples = meta->get_samples();
        auto positions = meta->get_positions();
        auto chromosomes = meta->get_chromosomes();

        // buffer string for writing
        std::string line_buffer;
        line_buffer.reserve(line_buffer_capcity);

        // open file
        BGZF *fp_out = bgzf_open(raw_ibd_fn, "w");
        verify(fp_out != NULL && "bgzf_open failed to open ibd out file");

        verify(bgzf_index_build_init(fp_out) == 0);

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
                    verify(
                        bgzf_write(fp_out, line_buffer.c_str(), line_buffer.size())
                            == line_buffer.size()
                        && "bgzf_write error during writting ibd_rec during decoding");

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
                        positions.get_cm(rec.get_pid2())
                            - positions.get_cm(rec.get_pid1()));

                    //	std::cout << "+" << line_buffer.size() << '\n';
                }
            }

        } while (did_vec_read_full);

        // flush the line_buffer
        if (line_buffer.size() > 0) {
            verify(bgzf_write(fp_out, line_buffer.c_str(), line_buffer.size())
                       == line_buffer.size()
                   && "bgzf_write error during writting ibd_rec during decoding");
            line_buffer.clear();
        }

        // close file and dump index
        verify(bgzf_index_dump(fp_out, raw_ibd_fn, ".gzi") == 0);
        bgzf_close(fp_out);
    }

    void
    write_to_file(size_t max = 0)
    {
        verify(fp != NULL);
        size_t total_bytes = 0;
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

        verify(bgzf_write(fp, &ibd_vec[0], total_bytes) == total_bytes);
        __used(total_bytes);
    }

    void
    write_to_file(IbdFile &other, size_t max = 0)
    {
        verify(other.get_fp() != NULL);
        size_t total_bytes;
        if (max == 0 || max > ibd_vec.size())
            total_bytes = ibd_vec.size() * sizeof(decltype(ibd_vec)::value_type);
        else
            total_bytes = max * sizeof(decltype(ibd_vec)::value_type);

        // std::cout << "write to other file: total bytes " << total_bytes << '\n';

        verify(bgzf_write(other.fp, &ibd_vec[0], total_bytes) == total_bytes);
        __used(total_bytes);
    }

    // if new_capcity >= current capcity then enlarge the capacity
    // if append is false
    // @return  true if vector is full; false if vector not full, which
    // indicates end of file
    bool
    read_from_file(bool append = false, size_t new_capacity = 0)
    {
        verify(fp != NULL);
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
        exit_on_false(bytes_read >= 0, "bgzf_read error", __FILE__, __LINE__);

        // shrink if fewer bytes are read
        if ((size_t) bytes_read < bytes_capacity) {
            ibd_vec.resize(n_element + bytes_read / element_size);

            // std::cout << ibd_vec.size() << "<---  read ibd_vec.size() \n";
            // ibd_vec[0].print();
            // ibd_vec[1].print();
            // ibd_vec[2].print();

            // make sure this happens only at end of the file
            char c;
            verify(bgzf_read(fp, &c, 1) == 0 && "bgzf_read read error");
            __used(c);
            return false;
        }

        return true;
    }

    void
    delete_file_from_disk()
    {
        std::filesystem::remove(filename);
    }

    bool
    has_equal_vector(const IbdFile &other)
    {
        /*
        std::cout << (ibd_vec == other.ibd_vec) << (meta.is_equal(other.meta))
                  << "<- ibdfile is equal \n";
        std::cout << ibd_vec.size() << " size.vs " << other.ibd_vec.size() << '\n';
        */

        return ibd_vec == other.ibd_vec;
    }

    void
    head(size_t n)
    {
        for (size_t i = 0; i < n && i < ibd_vec.size(); i++) {
            auto &rec = ibd_vec[i];
            std::cout << rec.get_sid1() << '\t' << rec.get_sid2() << '\t'
                      << rec.get_pid1() << '\t' << rec.get_pid2() << '\n';
        }
    }

    void
    tail(size_t n)
    {
        if (n > ibd_vec.size())
            n = ibd_vec.size();
        for (size_t i = ibd_vec.size() - n; i < ibd_vec.size(); i++) {
            auto &rec = ibd_vec[i];
            std::cout << rec.get_sid1() << '\t' << rec.get_sid2() << '\t'
                      << rec.get_pid1() << '\t' << rec.get_pid2() << '\n';
        }
    }

    void
    summary()
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
    print_to_file(std::ofstream &ofs)
    {
        for (auto rec : ibd_vec) {
            ofs << rec.get_sid1() << '\t' << rec.get_sid2() << '\t' << rec.get_pid1()
                << '\t' << rec.get_pid2() << '\n';
        }
    }
};

#endif
