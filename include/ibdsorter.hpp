#ifndef __ibdsorter_hpp__
#define __ibdsorter_hpp__
#include "common.hpp"
#include "ibdfile.hpp"
// #include <execution>

// Note:
// Be sure to call constructor
// Explitly use auto & instead of auto for reference
// File id should not use vector size as it needs to erase already merged items.
class IbdSorter
{
    IbdFile in;
    IbdFile out;
    std::vector<IbdFile> chunks;

    size_t max_rec_allowed_by_ram;
    std::string chunk_prefix;
    std::string out_mode;
    size_t counter;

    ibd_rec_cmp_t Cmp;

  private:
    void
    merge_first_k_chunk(IbdFile &out_put, size_t kways)
    {
        // std::cerr << "Enter merge_first_k_chunk\n";
        // The actual K
        size_t K = (chunks.size() > kways) ? kways : chunks.size();

        // Buffer size
        // Open file for reading
        // std::cerr << "set buffer and open files \n";
        size_t buffer_size = max_rec_allowed_by_ram / (K + 3);
        for (size_t i = 0; i < K; i++) {
            chunks[i].get_vec().reserve(buffer_size);
            chunks[i].open("r");
            // std::cout << "\t From file: " << chunks[i].get_filename() << '\n';
        }
        // std::cout << "\t To file: " << out_put.get_filename() << '\n';

        /*
        chunks[0].read_from_file();
        std::cerr << "\t chunk_size " <<  chunks.size() << '\n';
        std::cerr << "\t fp " <<  chunks[0].get_fp()	<< '\n';


        chunks[0].get_vec()[0].print();
        */

        out_put.get_vec().reserve(3 * buffer_size);

        // std::cerr << "set indicator \n";
        std::vector<size_t> indicator_vec;
        indicator_vec.resize(K, 0);

        // max_val
        ibd_rec1_t max_val;
        max_val.maximize();

        auto get_value
            = [&indicator_vec, this, &max_val](size_t id) mutable -> ibd_rec1_t & {
            // std::cerr << "enter get_value \n";
            if (indicator_vec[id] != chunks[id].get_vec().size()) {
                // std::cout << "get value: " << id << '\n';
                return chunks[id].get_vec()[indicator_vec[id]++];
            } else {
                auto is_full = chunks[id].read_from_file(false);
                // std::cout << "is full: " << eof << " chunk size "
                //           << chunks[id].get_vec().size() << '\n';
                // insert senitel
                if (!is_full)
                    chunks[id].get_vec().push_back(max_val);
                indicator_vec[id] = 0;
                return chunks[id].get_vec()[indicator_vec[id]++];
            }
            indicator_vec[id]++;
        };

        // The Tournament tree
        size_t winner;
        TournamentTree ttree(K, max_val);
        std::vector<ibd_rec1_t> init_vals;
        init_vals.reserve(K);
        for (size_t i = 0; i < K; i++)
            init_vals.push_back(get_value(i));

        // initial run
        auto &new_chunk_out_vec = out_put.get_vec();
        new_chunk_out_vec.push_back(ttree.init_run(init_vals, winner));

        while (!(new_chunk_out_vec.back() == max_val)) {
            // debug
            // std::cerr << "new_chunk_out_vec.size(): " << new_chunk_out_vec.size() <<
            // '\n' ; std::cerr << "new_chunk_out_vec[0]" <<
            // (new_chunk_out_vec[0].get_sid1()) << '\n'; std::cerr << "indicator_vec[0]:
            // " <<  indicator_vec[0]  << '\n'; std::cerr << "indicator_vec[1]: " <<
            // indicator_vec[1]  << '\n'; std::cerr << "indicator_vec[2]: " <<
            // indicator_vec[2]  << '\n'; flush out buffer if full
            if (new_chunk_out_vec.size() == new_chunk_out_vec.capacity()) {
                out_put.write_to_file();
                new_chunk_out_vec.clear();
            }
            // replace run
            new_chunk_out_vec.push_back(ttree.replace_run(get_value(winner), winner));

            // new_chunk_out_vec.back().print();
            // if (new_chunk_out_vec.size() > 5)
            //     verify(new_chunk_out_vec[new_chunk_out_vec.size() - 4]
            //                < new_chunk_out_vec[new_chunk_out_vec.size() - 3]
            //            || new_chunk_out_vec[new_chunk_out_vec.size() - 4]
            //                   == new_chunk_out_vec[new_chunk_out_vec.size() - 3]);
        }

        // flush out out_chunk records
        if (new_chunk_out_vec.size() > 1) {
            // remove the max_val senital
            out_put.write_to_file(new_chunk_out_vec.size() - 1);
            new_chunk_out_vec.clear();
        }

        // release memory and close file
        for (size_t i = 0; i < K; i++) {
            chunks[i].get_vec().resize(0);
            chunks[i].get_vec().shrink_to_fit();
            // std::cout << "close chunks[i]: " << chunks[i].get_filename() << '\n';
            chunks[i].close();
            // std::cout << "close chunks[i]: " << chunks[i].get_filename() << '\n';
        }
        out_put.get_vec().resize(0);
        out_put.get_vec().shrink_to_fit();
    }

  public:
    IbdSorter(const char *in_fn, const char *out_fn, const char *out_mode,
        const char *chunk_fn_prefix, size_t max_rec_allowed_by_ram,
        ibd_rec_cmp_t Cmp = std::less<ibd_rec1_t>())
        : max_rec_allowed_by_ram(max_rec_allowed_by_ram), out_mode(out_mode), Cmp(Cmp)
    {
        in = IbdFile(in_fn, NULL, max_rec_allowed_by_ram);
        out = IbdFile(out_fn, NULL, 0);
        chunk_prefix = chunk_fn_prefix;
        // std::cout << this->max_rec_allowed_by_ram << '\n';
        counter = 0;
    }

    void
    sort_into_chunks()
    {
        // std::cout << " sortting \n";

        auto &vec = in.get_vec();
        bool read_full = false;
        in.open("r");

        do {
            // fill in buffer to sort
            read_full = in.read_from_file(false);

            // sort
            // std::sort(std::execution::par_unseq, vec.begin(), vec.end());
            std::sort(vec.begin(), vec.end(), Cmp);

            // debug
            // std::cout << "In file vector size: " << in.get_vec().size() << '\n';

            // add a new file to store tempory files
            std::string chunk_file_name = chunk_prefix;
            chunk_file_name += std::to_string(counter++);
            chunks.push_back(IbdFile(chunk_file_name.c_str(), NULL, 0));
            auto &chunk = chunks.back();

            // open the file, write and close
            chunk.open("w");
            in.write_to_file(chunk);
            chunk.close();

        } while (read_full);

        // release memory
        in.get_vec().resize(0);
        in.get_vec().shrink_to_fit();

        in.close();
    }

    void
    merge_chunks(size_t kways = 10)
    {
        // std::cout << " merging\n";
        while (chunks.size() > kways) {

            // new chunk file
            std::string fn = chunk_prefix;
            fn += std::to_string(counter++);
            IbdFile new_chunk(fn.c_str(), NULL, 0);

            new_chunk.open("w");
            // std::cout << " Merge to " << fn << '\n';
            merge_first_k_chunk(new_chunk, kways);
            new_chunk.close();

            chunks.push_back(std::move(new_chunk));
            // for (size_t i = 0; i < chunks.size(); i++) {
            //      std::cout << i << " : name: " << chunks[i].get_filename() << '\n';
            // }
            // delete file from the chunk vector and from disk.
            for_each(chunks.begin(), chunks.begin() + kways,
                [](auto &chunk) { chunk.delete_file_from_disk(); });
            chunks.erase(chunks.begin(), chunks.begin() + kways);
        }

        out.open(out_mode.c_str());
        merge_first_k_chunk(out, kways);
        out.close();

        for_each(chunks.begin(), chunks.end(),
            [](auto &chunk) { chunk.delete_file_from_disk(); });
        chunks.erase(chunks.begin(), chunks.end());
    }
};

#endif
