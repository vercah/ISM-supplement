// C++ implementation skeleton for the next version.

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <omp.h>

namespace fs = std::filesystem;

static constexpr int DEFAULT_BLOCK_SIZE = 256;

static std::string trim_copy(const std::string& value) {
    size_t start = 0;
    while (start < value.size() &&
           std::isspace(static_cast<unsigned char>(value[start]))) {
        ++start;
    }

    size_t end = value.size();
    while (end > start &&
           std::isspace(static_cast<unsigned char>(value[end - 1]))) {
        --end;
    }

    return value.substr(start, end - start);
}

static std::string basename_posix_like(const std::string& value) {
    size_t end = value.size();
    while (end > 0 && value[end - 1] == '/') {
        --end;
    }
    if (end == 0) {
        return "";
    }

    const size_t slash = value.rfind('/', end - 1);
    const size_t start = (slash == std::string::npos) ? 0 : slash + 1;
    return value.substr(start, end - start);
}

static std::string python_path_stem(const std::string& name) {
    if (name.empty()) {
        return name;
    }

    bool all_dots = true;
    for (char c : name) {
        if (c != '.') {
            all_dots = false;
            break;
        }
    }
    if (all_dots) {
        return name;
    }

    if (name.back() == '.') {
        return name.substr(0, name.size() - 1);
    }

    const size_t dot = name.find_last_of('.');
    if (dot == std::string::npos || dot == 0) {
        return name;
    }

    return name.substr(0, dot);
}

static std::string canonical_sample_name(const std::string& value) {
    static const std::vector<std::string> FASTA_SUFFIXES = {
        ".fa.gz", ".fasta.gz", ".fna.gz", ".fa", ".fasta", ".fna",
    };

    const std::string stripped = trim_copy(value);
    const std::string name = basename_posix_like(stripped);

    for (const std::string& suffix : FASTA_SUFFIXES) {
        if (name.size() >= suffix.size() &&
            name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
            return name.substr(0, name.size() - suffix.size());
        }
    }

    return python_path_stem(name);
}

static int parse_int_after_equals(const std::string& token) {
    const size_t pos = token.find('=');
    if (pos == std::string::npos || pos + 1 >= token.size()) {
        throw std::runtime_error("malformed token: " + token);
    }
    return std::stoi(token.substr(pos + 1));
}

struct Row {
    std::vector<int> colors;
    int64_t unitig_weight = 0;
    int64_t kmer_weight = 0;
};

class FulgorDump {
public:
    FulgorDump(const std::string& dumpdir, const std::string& dataset, int k, const std::string& outdir)
        : dumpdir_(dumpdir),
          dataset_(dataset),
          k_(k),
          outdir_(outdir),
          color_sets_file_(dumpdir_ / (dataset_ + "_k" + std::to_string(k_) + ".color_sets.txt")),
          metadata_file_(dumpdir_ / (dataset_ + "_k" + std::to_string(k_) + ".metadata.txt")),
          unitigs_file_(dumpdir_ / (dataset_ + "_k" + std::to_string(k_) + ".unitigs.fa")),
          filenames_file_(fs::path("01_datasets") / (dataset_ + ".txt")) {
        validate_files();
        fs::create_directories(outdir_);
        load_names();
        load_color_sets();
        stream_unitigs_and_accumulate_weights();
        append_num_kmers_if_missing(total_kmers_);
    }

    const std::vector<std::string>& names() const { return names_; }
    const std::vector<Row>& rows() const { return rows_; }
    int num_samples() const { return static_cast<int>(names_.size()); }
    int64_t total_kmers() const { return total_kmers_; }

    std::string output_file(const std::string& kind) const {
        return (outdir_ / (dataset_ + "_k" + std::to_string(k_) + "_" + kind + ".dists.txt")).string();
    }

private:
    void validate_files() const {
        const std::vector<fs::path> paths = {
            color_sets_file_, metadata_file_, unitigs_file_, filenames_file_,
        };
        for (const fs::path& path : paths) {
            if (!fs::exists(path)) {
                throw std::runtime_error("missing required input: " + path.string());
            }
        }
    }

    void load_names() {
        std::ifstream in(filenames_file_);
        if (!in) {
            throw std::runtime_error("failed to open " + filenames_file_.string());
        }

        std::string line;
        while (std::getline(in, line)) {
            line = trim_copy(line);
            if (line.empty() || line[0] == '#') {
                continue;
            }
            names_.push_back(canonical_sample_name(line));
        }
    }

    void load_color_sets() {
        std::ifstream in(color_sets_file_);
        if (!in) {
            throw std::runtime_error("failed to open " + color_sets_file_.string());
        }

        std::string line;
        while (std::getline(in, line)) {
            line = trim_copy(line);
            if (line.empty()) {
                continue;
            }

            std::istringstream iss(line);
            std::string first;
            std::string second;
            if (!(iss >> first)) {
                continue;
            }

            const int color_set_id = parse_int_after_equals(first);
            std::vector<int> colors;
            if (iss >> second) {
                std::string rest;
                std::getline(iss, rest);
                rest = trim_copy(rest);
                if (!rest.empty()) {
                    std::istringstream colors_stream(rest);
                    int value = 0;
                    while (colors_stream >> value) {
                        colors.push_back(value);
                    }
                    std::sort(colors.begin(), colors.end());
                    colors.erase(std::unique(colors.begin(), colors.end()), colors.end());
                }
            }

            const auto it = row_of_color_set_id_.find(color_set_id);
            if (it == row_of_color_set_id_.end()) {
                row_of_color_set_id_[color_set_id] = rows_.size();
                rows_.push_back(Row{std::move(colors), 0, 0});
            } else {
                rows_[it->second].colors = std::move(colors);
            }
        }
    }

    void stream_unitigs_and_accumulate_weights() {
        std::ifstream in(unitigs_file_);
        if (!in) {
            throw std::runtime_error("failed to open " + unitigs_file_.string());
        }

        std::string line;
        int current_color_set_id = -1;
        bool have_color_set_id = false;

        while (std::getline(in, line)) {
            line = trim_copy(line);
            if (line.empty()) {
                continue;
            }

            if (line[0] == '>') {
                std::istringstream iss(line.substr(1));
                std::string token;
                if (!(iss >> token)) {
                    continue;
                }
                if (!(iss >> token)) {
                    continue;
                }
                current_color_set_id = parse_int_after_equals(token);
                have_color_set_id = true;
                continue;
            }

            const int n_kmers = std::max(0, static_cast<int>(line.size()) - k_ + 1);
            if (n_kmers <= 0) {
                continue;
            }

            total_kmers_ += n_kmers;

            if (!have_color_set_id) {
                continue;
            }

            const auto it = row_of_color_set_id_.find(current_color_set_id);
            if (it == row_of_color_set_id_.end()) {
                continue;
            }

            Row& row = rows_[it->second];
            row.unitig_weight += 1;
            row.kmer_weight += n_kmers;
        }
    }

    void append_num_kmers_if_missing(int64_t total_kmers) const {
        std::ifstream in(metadata_file_);
        if (!in) {
            throw std::runtime_error("failed to open " + metadata_file_.string());
        }

        std::string line;
        bool present = false;
        while (std::getline(in, line)) {
            if (line.rfind("num_kmers=", 0) == 0) {
                present = true;
                break;
            }
        }

        if (present) {
            return;
        }

        std::ofstream out(metadata_file_, std::ios::app);
        if (!out) {
            throw std::runtime_error("failed to append to " + metadata_file_.string());
        }
        out << "num_kmers=" << total_kmers << '\n';
    }

    fs::path dumpdir_;
    std::string dataset_;
    int k_;
    fs::path outdir_;
    fs::path color_sets_file_;
    fs::path metadata_file_;
    fs::path unitigs_file_;
    fs::path filenames_file_;

    std::vector<std::string> names_;
    std::vector<Row> rows_;
    std::unordered_map<int, size_t> row_of_color_set_id_;
    int64_t total_kmers_ = 0;
};

class BlockwiseDistanceWriter {
public:
    explicit BlockwiseDistanceWriter(const FulgorDump& dump, int block_size)
        : dump_(dump), block_size_(block_size) {
        if (block_size_ <= 0) {
            throw std::runtime_error("block size must be positive");
        }
    }

    void write_all() const {
        std::vector<int64_t> unitig_marginal;
        std::vector<int64_t> kmer_marginal;
        std::vector<int64_t> uniqrow_marginal;
        compute_marginals(unitig_marginal, kmer_marginal, uniqrow_marginal);

        const int n = dump_.num_samples();
        const std::vector<std::string>& names = dump_.names();
        const int num_blocks = (n + block_size_ - 1) / block_size_;
        std::vector<std::string> name_prefixes;
        name_prefixes.reserve(static_cast<size_t>(n));
        size_t total_name_chars = 0;
        for (const std::string& name : names) {
            total_name_chars += name.size();
            name_prefixes.push_back(name + '\t');
        }
        const size_t avg_name_chars =
            n > 0 ? total_name_chars / static_cast<size_t>(n) : 0;

        const fs::path unitig_path = dump_.output_file("unitig");
        const fs::path kmer_path = dump_.output_file("kmer");
        const fs::path uniqrow_path = dump_.output_file("uniqrow");
        std::vector<BlockTextResult> block_texts(static_cast<size_t>(num_blocks));

        std::atomic<bool> failed{false};
        std::string error_message;

#pragma omp parallel for schedule(dynamic)
        for (int block_index = 0; block_index < num_blocks; ++block_index) {
            if (failed.load(std::memory_order_relaxed)) {
                continue;
            }

            try {
                const int i0 = block_index * block_size_;
                const int i1 = std::min(i0 + block_size_, n);
                const BlockResult block = compute_block(i0, i1);
                if (failed.load(std::memory_order_relaxed)) {
                    continue;
                }
                block_texts[static_cast<size_t>(block_index)] =
                    render_block_text(block, names, name_prefixes,
                                      unitig_marginal, kmer_marginal, uniqrow_marginal,
                                      i0, i1, n, avg_name_chars);
            } catch (const std::exception& e) {
#pragma omp critical
                {
                    if (!failed.load(std::memory_order_relaxed)) {
                        failed.store(true, std::memory_order_relaxed);
                        error_message = e.what();
                    }
                }
            } catch (...) {
#pragma omp critical
                {
                    if (!failed.load(std::memory_order_relaxed)) {
                        failed.store(true, std::memory_order_relaxed);
                        error_message = "unknown error while computing block output";
                    }
                }
            }
        }

        if (failed.load(std::memory_order_relaxed)) {
            throw std::runtime_error(error_message);
        }

        std::ofstream unitig_out(unitig_path, std::ios::binary);
        std::ofstream kmer_out(kmer_path, std::ios::binary);
        std::ofstream uniqrow_out(uniqrow_path, std::ios::binary);
        if (!unitig_out || !kmer_out || !uniqrow_out) {
            throw std::runtime_error("failed to open output files");
        }

        for (const BlockTextResult& block_text : block_texts) {
            unitig_out.write(block_text.unitig_text.data(),
                             static_cast<std::streamsize>(block_text.unitig_text.size()));
            kmer_out.write(block_text.kmer_text.data(),
                           static_cast<std::streamsize>(block_text.kmer_text.size()));
            uniqrow_out.write(block_text.uniqrow_text.data(),
                              static_cast<std::streamsize>(block_text.uniqrow_text.size()));
        }

        if (!unitig_out || !kmer_out || !uniqrow_out) {
            throw std::runtime_error("failed while writing output files");
        }
    }

private:
    struct BlockResult {
        std::vector<int64_t> unitig;
        std::vector<int64_t> kmer;
        std::vector<int64_t> uniqrow;
    };

    struct BlockTextResult {
        std::string unitig_text;
        std::string kmer_text;
        std::string uniqrow_text;
    };

    void compute_marginals(std::vector<int64_t>& unitig,
                           std::vector<int64_t>& kmer,
                           std::vector<int64_t>& uniqrow) const {
        const int n = dump_.num_samples();
        unitig.assign(n, 0);
        kmer.assign(n, 0);
        uniqrow.assign(n, 0);

        for (const Row& row : dump_.rows()) {
            for (int sample : row.colors) {
                if (sample < 0 || sample >= n) {
                    throw std::runtime_error("sample index out of range");
                }
                unitig[sample] += row.unitig_weight;
                kmer[sample] += row.kmer_weight;
                uniqrow[sample] += 1;
            }
        }
    }

    BlockResult compute_block(int i0, int i1) const {
        const int n = dump_.num_samples();
        const size_t block_len = static_cast<size_t>(i1 - i0);
        const size_t block_area = block_len * static_cast<size_t>(n);

        BlockResult result;
        result.unitig.assign(block_area, 0);
        result.kmer.assign(block_area, 0);
        result.uniqrow.assign(block_area, 0);

        const std::vector<Row>& rows = dump_.rows();
        for (const Row& row : rows) {
            if (row.colors.empty()) {
                continue;
            }

            if (row.colors.front() >= i1 || row.colors.back() < i0) {
                continue;
            }

            const auto lo = std::lower_bound(row.colors.begin(), row.colors.end(), i0);
            if (lo == row.colors.end()) {
                continue;
            }

            const auto hi = std::lower_bound(lo, row.colors.end(), i1);
            if (lo == hi) {
                continue;
            }

            const int64_t unitig_weight = row.unitig_weight;
            const int64_t kmer_weight = row.kmer_weight;

            for (auto it = lo; it != hi; ++it) {
                const size_t offset = static_cast<size_t>(*it - i0) * static_cast<size_t>(n);
                int64_t* unitig_row = result.unitig.data() + offset;
                int64_t* kmer_row = result.kmer.data() + offset;
                int64_t* uniqrow_row = result.uniqrow.data() + offset;

                // row.colors is sorted, so the suffix after i is exactly the j > i frontier.
                for (auto jt = std::next(it); jt != row.colors.end(); ++jt) {
                    const int sample = *jt;
                    unitig_row[sample] += unitig_weight;
                    kmer_row[sample] += kmer_weight;
                    uniqrow_row[sample] += 1;
                }
            }
        }

        return result;
    }

    static void append_distance_line(std::string& out,
                                     const std::string& name_prefix,
                                     const std::string& other_name,
                                     int64_t dist) {
        out.append(name_prefix);
        out.append(other_name);
        out.push_back('\t');
        out.append(std::to_string(dist));
        out.push_back('\n');
    }

    static size_t estimate_block_text_bytes(int i0,
                                            int i1,
                                            int n,
                                            size_t avg_name_chars) {
        const size_t block_len = static_cast<size_t>(i1 - i0);
        const size_t pair_count =
            block_len * static_cast<size_t>(2 * n - i0 - i1 - 1) / 2;
        const size_t approx_line_bytes = 2 * avg_name_chars + 24;
        return pair_count * approx_line_bytes;
    }

    static BlockTextResult render_block_text(const BlockResult& block,
                                             const std::vector<std::string>& names,
                                             const std::vector<std::string>& name_prefixes,
                                             const std::vector<int64_t>& unitig_marginal,
                                             const std::vector<int64_t>& kmer_marginal,
                                             const std::vector<int64_t>& uniqrow_marginal,
                                             int i0,
                                             int i1,
                                             int n,
                                             size_t avg_name_chars) {
        BlockTextResult text;
        const size_t reserve_bytes = estimate_block_text_bytes(i0, i1, n, avg_name_chars);
        text.unitig_text.reserve(reserve_bytes);
        text.kmer_text.reserve(reserve_bytes);
        text.uniqrow_text.reserve(reserve_bytes);

        for (int i = i0; i < i1; ++i) {
            const std::string& name_prefix = name_prefixes[static_cast<size_t>(i)];
            const int local_i = i - i0;
            const int64_t* unitig_row =
                block.unitig.data() + static_cast<size_t>(local_i) * static_cast<size_t>(n);
            const int64_t* kmer_row =
                block.kmer.data() + static_cast<size_t>(local_i) * static_cast<size_t>(n);
            const int64_t* uniqrow_row =
                block.uniqrow.data() + static_cast<size_t>(local_i) * static_cast<size_t>(n);
            for (int j = i + 1; j < n; ++j) {
                const int64_t unitig_dist =
                    unitig_marginal[i] + unitig_marginal[j] - 2 * unitig_row[j];
                const int64_t kmer_dist =
                    kmer_marginal[i] + kmer_marginal[j] - 2 * kmer_row[j];
                const int64_t uniqrow_dist =
                    uniqrow_marginal[i] + uniqrow_marginal[j] - 2 * uniqrow_row[j];
                const std::string& other_name = names[static_cast<size_t>(j)];
                append_distance_line(text.unitig_text, name_prefix, other_name, unitig_dist);
                append_distance_line(text.kmer_text, name_prefix, other_name, kmer_dist);
                append_distance_line(text.uniqrow_text, name_prefix, other_name, uniqrow_dist);
            }
        }

        return text;
    }

    const FulgorDump& dump_;
    int block_size_;
};

static void print_usage(const char* argv0) {
    std::cerr << "usage: " << argv0
              << " dumpdir dataset k out_prefix n_cores [--block-size N]\n";
}

int main(int argc, char** argv) {
    try {
        if (argc != 6 && argc != 8) {
            print_usage(argv[0]);
            return 1;
        }

        std::string dumpdir = argv[1];
        std::string dataset = argv[2];
        int k = std::stoi(argv[3]);
        std::string out_prefix = argv[4];
        int n_cores = std::stoi(argv[5]);
        int block_size = DEFAULT_BLOCK_SIZE;

        if (argc == 8) {
            if (std::string(argv[6]) != "--block-size") {
                print_usage(argv[0]);
                return 1;
            }
            block_size = std::stoi(argv[7]);
        }

        if (n_cores < 1) {
            n_cores = 1;
        }

        omp_set_dynamic(0);
        omp_set_num_threads(n_cores);

        FulgorDump dump(dumpdir, dataset, k, out_prefix);
        BlockwiseDistanceWriter writer(dump, block_size);
        writer.write_all();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
