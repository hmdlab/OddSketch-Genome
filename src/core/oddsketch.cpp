#include "oddsketch/oddsketch_cli.hpp"
#include "third_party/xxhash.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include <zlib.h>

namespace {

namespace fs = std::filesystem;

constexpr uint64_t kHashSeed = 0x9E3779B97F4A7C15ULL;
constexpr uint64_t kBucketHashSeed = 0xD1B54A32D192ED03ULL;
constexpr uint64_t kRankHashSeed = 0xABC98388FB8FAC03ULL;
constexpr uint64_t kEmptyBucket = std::numeric_limits<uint64_t>::max();
constexpr uint64_t kBaseSeeds[4] = {
    0x3C8BFBB395C60474ULL, // A
    0x3193C18562A02B4CULL, // C
    0x20323ED082572324ULL, // G
    0x295549F54BE24456ULL, // T
};

enum class PosMode { Value = 0, Mix = 1, Stripe = 2 };

// Runtime parameters supplied by the CLI or experiment scripts.
struct OddsketchOptions {
    size_t sketch_size = 8192; // Multiple of 64.
    size_t kmer = 64;
    bool canonical = true;
    double j0 = 0.75;
    PosMode pos_mode = PosMode::Value;
    size_t threads = 1;
};

struct CliArgs {
    std::string mode;
    std::string input_paths_path;
    std::string qlist_path;
    std::string dblist_path;
    std::string pairlist_path;
    std::string out_dir;
    std::string sketch_paths_out;
    bool explicit_all_to_all = false;
    bool explicit_bipartite = false;
    bool skip_existing = false;
    OddsketchOptions options;
};

struct SketchInput {
    std::string sequence_path;
};

// Densified MinHash values from OPH and the actual bucket count used.
struct MinhashResult {
    std::vector<uint64_t> values;
    size_t num_buckets = 0;
};

// Sketch words plus metadata used by downstream distance estimation.
struct SketchBuildResult {
    std::vector<uint64_t> words;
    size_t num_buckets = 0;
};

// Compact on-disk header. The version field keeps compatibility explicit.
struct SketchHeader {
    char magic[4];
    uint32_t version;
    uint64_t nbits;
    uint64_t kbuckets;
    double j0;
};

struct Sketch {
    std::vector<uint64_t> words;
    uint64_t bit_size = 0;
    uint64_t k_buckets = 0;
    double j0 = 0.75;
};

struct SketchPair {
    std::string left;
    std::string right;
};

struct PairlistGroup {
    std::string left;
    std::vector<size_t> pair_indices;
};

struct UnivHash {
    uint64_t a;
    uint64_t b;
    size_t k;
    unsigned shift;

    explicit UnivHash(size_t k_, uint64_t seed = kHashSeed) : a(0), b(0), k(k_), shift(64) {
        a = 0x2545F4914F6CDD1DULL ^ seed;
        if ((a & 1ULL) == 0) {
            a ^= 1ULL;
        }
        b = 0x9E3779B185EBCA87ULL + (seed << 1);
        unsigned m = 0;
        while ((1ULL << m) < k) {
            ++m;
        }
        shift = 64 - m;
    }

    size_t operator()(uint64_t x) const {
        return static_cast<size_t>((a * x + b) >> shift);
    }
};

uint64_t rotl64(uint64_t value, size_t count) {
    const unsigned r = static_cast<unsigned>(count & 63U);
    if (r == 0) {
        return value;
    }
    return (value << r) | (value >> (64U - r));
}

uint64_t rotr64(uint64_t value, size_t count) {
    const unsigned r = static_cast<unsigned>(count & 63U);
    if (r == 0) {
        return value;
    }
    return (value >> r) | (value << (64U - r));
}

uint64_t mix64(uint64_t value) {
    value += 0x9E3779B97F4A7C15ULL;
    value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31);
}

int base_index(char c) {
    switch (c) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return -1;
    }
}

char normalize_base(char c) {
    switch (c) {
        case 'A':
        case 'a':
            return 'A';
        case 'C':
        case 'c':
            return 'C';
        case 'G':
        case 'g':
            return 'G';
        case 'T':
        case 't':
            return 'T';
        default:
            return '\0';
    }
}

size_t floor_pow2(size_t x) {
    if (x == 0) {
        return 1;
    }

    size_t p = 1;
    while ((p << 1) <= x) {
        p <<= 1;
    }
    return p;
}

bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

size_t compute_num_buckets(const OddsketchOptions& options) {
    // Choose the OPH bucket count from L = n / (4 * (1 - J0)) in the OddSketch paper.
    // Round down to a power of two for uniformity and simpler indexing.
    const double denom = 4.0 * (1.0 - options.j0);
    const size_t desired = (denom > 0.0)
        ? static_cast<size_t>(std::max(1.0, std::round(static_cast<double>(options.sketch_size) / denom)))
        : options.sketch_size;
    return std::max<size_t>(1, floor_pow2(desired));
}

size_t normalize_thread_count(size_t requested_threads, size_t task_count) {
    const size_t safe_threads = std::max<size_t>(1, requested_threads);
    return task_count == 0 ? 1 : std::min(safe_threads, task_count);
}

template <typename Func>
void run_in_parallel(size_t task_count, size_t requested_threads, Func&& func) {
    const size_t worker_count = normalize_thread_count(requested_threads, task_count);
    if (task_count == 0) {
        return;
    }
    if (worker_count == 1) {
        for (size_t i = 0; i < task_count; ++i) {
            func(i);
        }
        return;
    }

    std::atomic<size_t> next_index{0};
    std::atomic<bool> failed{false};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    auto worker = [&]() {
        while (!failed.load(std::memory_order_relaxed)) {
            const size_t idx = next_index.fetch_add(1, std::memory_order_relaxed);
            if (idx >= task_count) {
                break;
            }
            try {
                func(idx);
            } catch (...) {
                {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!first_error) {
                        first_error = std::current_exception();
                    }
                }
                failed.store(true, std::memory_order_relaxed);
                break;
            }
        }
    };

    for (size_t i = 0; i < worker_count; ++i) {
        workers.emplace_back(worker);
    }
    for (auto& thread : workers) {
        thread.join();
    }

    if (first_error) {
        std::rethrow_exception(first_error);
    }
}

uint64_t pack_pair(uint32_t i, uint32_t attempt) {
    return (uint64_t(i) << 32) ^ uint64_t(attempt);
}

std::vector<uint64_t> densify_optimal(const std::vector<uint64_t>& buckets) {
    const size_t k = buckets.size();
    std::vector<uint64_t> out(k);
    if (k == 0) {
        return out;
    }

    const bool any = std::any_of(buckets.begin(), buckets.end(), [](uint64_t v) {
        return v != kEmptyBucket;
    });
    if (!any) {
        std::fill(out.begin(), out.end(), 0ULL);
        return out;
    }

    // Fill empty buckets using optimal densification from Shrivastava 2017.
    const UnivHash huniv(k);
    for (uint32_t i = 0; i < static_cast<uint32_t>(k); ++i) {
        if (buckets[i] != kEmptyBucket) {
            out[i] = buckets[i];
            continue;
        }

        uint32_t attempt = 1;
        while (true) {
            const size_t next = huniv(pack_pair(i, attempt));
            if (buckets[next] != kEmptyBucket) {
                out[i] = buckets[next];
                break;
            }
            ++attempt;
        }
    }
    return out;
}

void update_minhash_buckets_from_sequence(
    std::string_view seq,
    const OddsketchOptions& options,
    unsigned shift,
    std::vector<uint64_t>& bucket_min_hash
) {
    const size_t k = options.kmer;
    if (k == 0 || seq.size() < k) {
        return;
    }

    auto update_bucket = [&](uint64_t raw_hash) {
        const uint64_t h_bucket = mix64(raw_hash ^ kBucketHashSeed);
        uint64_t h_rank = mix64(raw_hash ^ kRankHashSeed);
        if (h_rank == kEmptyBucket) {
            --h_rank;
        }
        const size_t bucket_id = static_cast<size_t>(h_bucket >> shift);
        if (h_rank < bucket_min_hash[bucket_id]) {
            bucket_min_hash[bucket_id] = h_rank;
        }
    };

    uint64_t fwd = 0;
    uint64_t rev = 0;
    for (size_t i = 0; i < k; ++i) {
        const int idx = base_index(seq[i]);
        if (idx < 0) {
            return;
        }
        fwd ^= rotl64(kBaseSeeds[idx], k - 1 - i);
        rev ^= rotl64(kBaseSeeds[3 - idx], i);
    }
    update_bucket(options.canonical ? std::min(fwd, rev) : fwd);

    for (size_t i = k; i < seq.size(); ++i) {
        const int out_idx = base_index(seq[i - k]);
        const int in_idx = base_index(seq[i]);
        if (out_idx < 0 || in_idx < 0) {
            return;
        }

        fwd = rotl64(fwd, 1) ^ rotl64(kBaseSeeds[out_idx], k) ^ kBaseSeeds[in_idx];
        rev = rotr64(rev, 1) ^ rotr64(kBaseSeeds[3 - out_idx], 1) ^ rotl64(kBaseSeeds[3 - in_idx], k - 1);
        update_bucket(options.canonical ? std::min(fwd, rev) : fwd);
    }
}

void flip_bit(std::vector<uint64_t>& words, size_t idx) {
    const size_t w = idx / 64;
    const size_t b = idx % 64;
    words[w] ^= (uint64_t(1) << b);
}

size_t map_position(uint32_t idx, uint64_t hv, size_t nbits, size_t kbuckets, PosMode pos_mode) {
    if (pos_mode == PosMode::Value) {
        return static_cast<size_t>(hv % nbits);
    }

    if (pos_mode == PosMode::Stripe) {
        // Assign each bucket a local stripe of bit positions.
        const size_t stride = (kbuckets > 0) ? (nbits / kbuckets) : 0;
        if (stride >= 2) {
            const size_t base = static_cast<size_t>(idx) * stride;
            const size_t offset = static_cast<size_t>(hv % stride);
            const size_t pos = base + offset;
            return (pos < nbits) ? pos : (pos % nbits);
        }
    }

    // Mix recombines the bucket index and hash value to reduce positional correlation.
    uint64_t buf[2];
    buf[0] = hv;
    buf[1] = static_cast<uint64_t>(idx);
    const uint64_t mixed = XXH64(reinterpret_cast<const void*>(buf), sizeof(buf), kHashSeed);
    return static_cast<size_t>(mixed % nbits);
}

void update_from_sequence_line(
    const std::string& line,
    const OddsketchOptions& options,
    unsigned shift,
    std::vector<uint64_t>& bucket_min_hash,
    std::string& segment
) {
    auto flush_segment = [&]() {
        if (!segment.empty()) {
            update_minhash_buckets_from_sequence(segment, options, shift, bucket_min_hash);
            segment.clear();
        }
    };

    for (char c : line) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            continue;
        }
        const char base = normalize_base(c);
        if (base == '\0') {
            flush_segment();
            continue;
        }
        segment.push_back(base);
    }
}

template <typename LineFunc>
void process_sequence_line(
    std::string line,
    const OddsketchOptions& options,
    unsigned shift,
    std::vector<uint64_t>& bucket_min_hash,
    std::string& segment,
    bool& skip_quality_line,
    LineFunc&& update_line
) {
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    const auto first = std::find_if_not(line.begin(), line.end(), [](char c) {
        return std::isspace(static_cast<unsigned char>(c));
    });
    if (first == line.end()) {
        return;
    }

    auto flush_segment = [&]() {
        if (!segment.empty()) {
            update_minhash_buckets_from_sequence(segment, options, shift, bucket_min_hash);
            segment.clear();
        }
    };

    if (skip_quality_line) {
        skip_quality_line = false;
        flush_segment();
        return;
    }

    const char lead = *first;
    if (lead == '>' || lead == '@') {
        flush_segment();
        return;
    }
    if (lead == '+') {
        flush_segment();
        skip_quality_line = true;
        return;
    }

    update_line(line, options, shift, bucket_min_hash, segment);
}

void update_minhash_from_plain_sequence_file(
    const std::string& fname,
    const OddsketchOptions& options,
    unsigned shift,
    std::vector<uint64_t>& bucket_min_hash
) {
    std::ifstream ifs(fname);
    if (!ifs) {
        throw std::runtime_error("Cannot open " + fname);
    }

    std::string line;
    std::string segment;
    bool skip_quality_line = false;
    while (std::getline(ifs, line)) {
        process_sequence_line(
            line,
            options,
            shift,
            bucket_min_hash,
            segment,
            skip_quality_line,
            update_from_sequence_line
        );
    }
    if (!segment.empty()) {
        update_minhash_buckets_from_sequence(segment, options, shift, bucket_min_hash);
    }
}

void update_minhash_from_gzip_sequence_file(
    const std::string& fname,
    const OddsketchOptions& options,
    unsigned shift,
    std::vector<uint64_t>& bucket_min_hash
) {
    gzFile file = gzopen(fname.c_str(), "rb");
    if (file == nullptr) {
        throw std::runtime_error("Cannot open gzip FASTA: " + fname);
    }

    std::vector<char> buffer(1024 * 1024);
    std::string line;
    std::string segment;
    bool skip_quality_line = false;
    while (gzgets(file, buffer.data(), static_cast<int>(buffer.size())) != nullptr) {
        const std::string_view chunk(buffer.data());
        for (char c : chunk) {
            if (c == '\n') {
                process_sequence_line(
                    line,
                    options,
                    shift,
                    bucket_min_hash,
                    segment,
                    skip_quality_line,
                    update_from_sequence_line
                );
                line.clear();
            } else {
                line.push_back(c);
            }
        }
    }
    if (!line.empty()) {
        process_sequence_line(
            line,
            options,
            shift,
            bucket_min_hash,
            segment,
            skip_quality_line,
            update_from_sequence_line
        );
    }
    if (!segment.empty()) {
        update_minhash_buckets_from_sequence(segment, options, shift, bucket_min_hash);
    }

    int gzerrnum = Z_OK;
    const char* gzerr = gzerror(file, &gzerrnum);
    const int close_status = gzclose(file);
    if (gzerrnum != Z_OK && gzerrnum != Z_STREAM_END) {
        throw std::runtime_error("Error reading gzip FASTA " + fname + ": " + std::string(gzerr ? gzerr : "unknown error"));
    }
    if (close_status != Z_OK) {
        throw std::runtime_error("Error closing gzip FASTA: " + fname);
    }
}

MinhashResult get_minhash_one_permutation_from_sequence_file(
    const std::string& fname,
    const OddsketchOptions& options
) {
    const size_t num_buckets = compute_num_buckets(options);
    std::vector<uint64_t> bucket_min_hash(num_buckets, kEmptyBucket);

    unsigned m = 0;
    while ((1ULL << m) < num_buckets) {
        ++m;
    }
    const unsigned shift = 64 - m;

    if (ends_with(fname, ".gz")) {
        update_minhash_from_gzip_sequence_file(fname, options, shift, bucket_min_hash);
    } else {
        update_minhash_from_plain_sequence_file(fname, options, shift, bucket_min_hash);
    }
    return {densify_optimal(bucket_min_hash), num_buckets};
}

SketchBuildResult make_odd_sketch_from_fasta(const std::string& fname, const OddsketchOptions& options) {
    if (options.sketch_size == 0 || (options.sketch_size % 64) != 0) {
        throw std::runtime_error("Invalid sketch size (must be a positive multiple of 64)");
    }

    const MinhashResult minhash = get_minhash_one_permutation_from_sequence_file(fname, options);
    std::vector<uint64_t> words(options.sketch_size / 64, 0);
    for (size_t i = 0; i < minhash.values.size(); ++i) {
        const uint64_t hv = minhash.values[i];
        const size_t pos = map_position(
            static_cast<uint32_t>(i),
            hv,
            options.sketch_size,
            minhash.values.size(),
            options.pos_mode
        );
        // OddSketch updates the sketch by flipping the selected bit.
        flip_bit(words, pos);
    }

    return {std::move(words), minhash.num_buckets};
}

bool has_valid_sketch_header(const SketchHeader& h) {
    return h.magic[0] == 'O' &&
           h.magic[1] == 'D' &&
           h.magic[2] == 'S' &&
           h.magic[3] == 'K' &&
           h.version >= 1;
}

Sketch load_legacy_sketch(std::ifstream& ifs, const std::string& fname) {
    // The legacy format is just the word array without a header.
    ifs.clear();
    ifs.seekg(0, std::ios::end);
    const std::streampos end = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    if (end <= 0 || (end % static_cast<std::streampos>(sizeof(uint64_t))) != 0) {
        throw std::runtime_error("Invalid sketch file size: " + fname);
    }

    const size_t num_words = static_cast<size_t>(end / static_cast<std::streampos>(sizeof(uint64_t)));
    Sketch s;
    s.words.resize(num_words);
    ifs.read(reinterpret_cast<char*>(s.words.data()), num_words * sizeof(uint64_t));
    if (!ifs) {
        throw std::runtime_error("Error reading sketch file: " + fname);
    }
    s.bit_size = static_cast<uint64_t>(num_words * 64);
    return s;
}

void write_sketch_binary(const std::vector<uint64_t>& words,
                         const std::string& outfname,
                         uint64_t nbits,
                         uint64_t kbuckets,
                         double j0val) {
    std::ofstream ofs(outfname, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Cannot open output file: " + outfname);
    }

    SketchHeader h{};
    h.magic[0] = 'O';
    h.magic[1] = 'D';
    h.magic[2] = 'S';
    h.magic[3] = 'K';
    h.version = 1;
    h.nbits = nbits;
    h.kbuckets = kbuckets;
    h.j0 = j0val;

    ofs.write(reinterpret_cast<const char*>(&h), sizeof(h));
    if (!ofs) {
        throw std::runtime_error("Error writing header: " + outfname);
    }

    ofs.write(reinterpret_cast<const char*>(words.data()), words.size() * sizeof(uint64_t));
    if (!ofs) {
        throw std::runtime_error("Error writing body: " + outfname);
    }
}

Sketch load_sketch(const std::string& fname) {
    std::ifstream ifs(fname, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Cannot open sketch file: " + fname);
    }

    SketchHeader h{};
    ifs.read(reinterpret_cast<char*>(&h), sizeof(h));
    // Fall back to the legacy reader if the file does not start with a valid header.
    if (!ifs || !has_valid_sketch_header(h)) {
        return load_legacy_sketch(ifs, fname);
    }

    ifs.seekg(0, std::ios::end);
    const std::streampos end = ifs.tellg();
    const std::streampos body = end - static_cast<std::streampos>(sizeof(SketchHeader));
    if (body < 0 || (body % static_cast<std::streampos>(sizeof(uint64_t))) != 0) {
        throw std::runtime_error("Invalid sketch body size: " + fname);
    }

    const size_t num_words = static_cast<size_t>(body / static_cast<std::streampos>(sizeof(uint64_t)));
    Sketch s;
    s.words.resize(num_words);
    ifs.seekg(sizeof(SketchHeader), std::ios::beg);
    ifs.read(reinterpret_cast<char*>(s.words.data()), num_words * sizeof(uint64_t));
    if (!ifs) {
        throw std::runtime_error("Error reading sketch body: " + fname);
    }

    s.bit_size = h.nbits ? h.nbits : static_cast<uint64_t>(s.words.size() * 64);
    s.k_buckets = h.kbuckets;
    s.j0 = (h.j0 > 0.0 && h.j0 < 1.0) ? h.j0 : 0.75;
    return s;
}

std::vector<Sketch> load_all_sketches(const std::vector<std::string>& paths) {
    std::vector<Sketch> sketches;
    sketches.reserve(paths.size());
    for (const auto& path : paths) {
        sketches.push_back(load_sketch(path));
    }
    return sketches;
}

double jaccard_distance(const Sketch& a, const Sketch& b, const OddsketchOptions& options) {
    if (a.words.size() != b.words.size()) {
        throw std::runtime_error("Sketch size mismatch");
    }

    // Evaluate the paper's estimator from the number of differing XOR bits.
    uint64_t popcnt = 0;
    for (size_t w = 0; w < a.words.size(); ++w) {
        popcnt += __builtin_popcountll(a.words[w] ^ b.words[w]);
    }

    const double n = static_cast<double>(a.words.size() * 64);
    const uint64_t k_header = (a.k_buckets > 0 && b.k_buckets > 0)
        ? std::min(a.k_buckets, b.k_buckets)
        : (a.k_buckets > 0 ? a.k_buckets : b.k_buckets);

    double k = 0.0;
    // Prefer the generation-time k from the header; otherwise re-estimate it from j0.
    if (k_header > 0) {
        k = static_cast<double>(k_header);
    } else {
        const double denom = 4.0 * (1.0 - options.j0);
        const size_t k_est = (denom > 0.0)
            ? floor_pow2(static_cast<size_t>(std::max(1.0, std::round(n / denom))))
            : static_cast<size_t>(n);
        k = static_cast<double>(k_est);
    }

    const double x = 1.0 - (2.0 * static_cast<double>(popcnt)) / n;
    const double term = std::log(std::max(1e-12, x));
    double jacc = 1.0 + (n / (4.0 * k)) * term;
    jacc = std::clamp(jacc, 0.0, 1.0);
    return jacc;
}

std::vector<std::string> read_nonempty_lines(std::istream& input) {
    std::vector<std::string> paths;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty()) {
            paths.push_back(line);
        }
    }
    return paths;
}

std::vector<std::string> read_paths_from_stdin() {
    return read_nonempty_lines(std::cin);
}

std::vector<std::string> read_paths_from_list_file(const std::string& list_path) {
    std::ifstream ifs(list_path);
    if (!ifs) {
        throw std::runtime_error("Cannot open list file: " + list_path);
    }
    return read_nonempty_lines(ifs);
}

std::vector<SketchInput> read_sketch_inputs_from_stdin() {
    const auto paths = read_paths_from_stdin();
    std::vector<SketchInput> inputs;
    inputs.reserve(paths.size());
    for (const auto& path : paths) {
        inputs.push_back({path});
    }
    return inputs;
}

std::vector<SketchInput> make_sketch_inputs_from_paths(const std::vector<std::string>& paths) {
    std::vector<SketchInput> inputs;
    inputs.reserve(paths.size());
    for (const auto& path : paths) {
        inputs.push_back({path});
    }
    return inputs;
}

std::vector<SketchInput> read_sketch_inputs_from_input_paths_file(const std::string& input_paths_path) {
    return make_sketch_inputs_from_paths(read_paths_from_list_file(input_paths_path));
}

std::vector<SketchPair> read_pairs_from_pairlist_file(const std::string& pairlist_path) {
    std::ifstream ifs(pairlist_path);
    if (!ifs) {
        throw std::runtime_error("Cannot open pairlist file: " + pairlist_path);
    }

    std::vector<SketchPair> pairs;
    std::string line;
    size_t lineno = 0;
    while (std::getline(ifs, line)) {
        ++lineno;
        if (line.empty()) {
            continue;
        }

        const auto tab = line.find('\t');
        if (tab == std::string::npos) {
            throw std::runtime_error("Invalid pairlist line " + std::to_string(lineno) + ": expected two tab-separated paths");
        }

        const std::string left = line.substr(0, tab);
        const std::string right = line.substr(tab + 1);
        if (left.empty() || right.empty() || right.find('\t') != std::string::npos) {
            throw std::runtime_error("Invalid pairlist line " + std::to_string(lineno) + ": expected exactly two tab-separated paths");
        }
        pairs.push_back({left, right});
    }
    return pairs;
}

std::string sketch_output_path(const SketchInput& input, const CliArgs& args) {
    if (args.out_dir.empty()) {
        return input.sequence_path + ".sketch";
    }

    std::string filename = fs::path(input.sequence_path).filename().string();
    if (ends_with(filename, ".gz")) {
        filename.resize(filename.size() - 3);
    }
    return (fs::path(args.out_dir) / (filename + ".sketch")).string();
}

void ensure_parent_dir(const std::string& path) {
    const fs::path parent = fs::path(path).parent_path();
    if (!parent.empty()) {
        fs::create_directories(parent);
    }
}

void write_sketch_paths_file(const std::string& sketch_paths_out, const std::vector<std::string>& output_paths) {
    ensure_parent_dir(sketch_paths_out);
    std::ofstream ofs(sketch_paths_out);
    if (!ofs) {
        throw std::runtime_error("Cannot write sketch path list: " + sketch_paths_out);
    }
    for (const auto& output_path : output_paths) {
        ofs << output_path << "\n";
    }
}

void reject_duplicate_outputs(const std::vector<std::string>& output_paths) {
    std::unordered_map<std::string, size_t> seen;
    for (size_t i = 0; i < output_paths.size(); ++i) {
        const auto inserted = seen.emplace(output_paths[i], i);
        if (!inserted.second) {
            throw std::runtime_error(
                "Duplicate sketch output path: " + output_paths[i] +
                " (inputs " + std::to_string(inserted.first->second + 1) +
                " and " + std::to_string(i + 1) + ")"
            );
        }
    }
}

std::vector<PairlistGroup> group_pairs_by_left(const std::vector<SketchPair>& pairs) {
    std::unordered_map<std::string, size_t> group_index_by_left;
    std::vector<PairlistGroup> groups;
    groups.reserve(pairs.size());

    for (size_t i = 0; i < pairs.size(); ++i) {
        const auto& left = pairs[i].left;
        const auto found = group_index_by_left.find(left);
        if (found == group_index_by_left.end()) {
            const size_t new_group_index = groups.size();
            group_index_by_left.emplace(left, new_group_index);
            groups.push_back({left, {i}});
        } else {
            groups[found->second].pair_indices.push_back(i);
        }
    }
    return groups;
}

const char* pos_mode_name(PosMode pos_mode) {
    switch (pos_mode) {
        case PosMode::Value:
            return "value";
        case PosMode::Mix:
            return "mix";
        case PosMode::Stripe:
            return "stripe";
    }
    return "value";
}

void run_dist_all_vs_all(const std::vector<std::string>& paths, const OddsketchOptions& options) {
    const auto sketches = load_all_sketches(paths);
    if (!sketches.empty()) {
        const auto& s0 = sketches.front();
        const uint64_t nbits = s0.bit_size ? s0.bit_size : static_cast<uint64_t>(s0.words.size() * 64);
        std::cerr << "[oddsketch] dist(all-to-all): nbits=" << nbits
                  << ", threads=" << normalize_thread_count(options.threads, sketches.size())
                  << ", kbuckets(header)=" << s0.k_buckets
                  << ", j0(current)=" << std::fixed << std::setprecision(6) << options.j0 << "\n";
    }

    std::mutex output_mutex;
    run_in_parallel(sketches.size(), options.threads, [&](size_t i) {
        std::ostringstream chunk;
        for (size_t j = i + 1; j < sketches.size(); ++j) {
            const double d = jaccard_distance(sketches[i], sketches[j], options);
            chunk << paths[i] << '\t' << paths[j] << '\t' << d << "\n";
        }
        const std::string text = chunk.str();
        if (!text.empty()) {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << text;
        }
    });
}

void run_dist_bipartite(const std::vector<std::string>& qpaths,
                        const std::vector<std::string>& dbpaths,
                        const OddsketchOptions& options) {
    const auto qsketches = load_all_sketches(qpaths);
    const auto dbsketches = load_all_sketches(dbpaths);
    if (!qsketches.empty()) {
        const auto& s0 = qsketches.front();
        const uint64_t nbits = s0.bit_size ? s0.bit_size : static_cast<uint64_t>(s0.words.size() * 64);
        std::cerr << "[oddsketch] dist(bipartite): nbits=" << nbits
                  << ", queries=" << qsketches.size()
                  << ", db=" << dbsketches.size()
                  << ", threads=" << normalize_thread_count(options.threads, qsketches.size())
                  << ", kbuckets(header)=" << s0.k_buckets
                  << ", j0(current)=" << std::fixed << std::setprecision(6) << options.j0 << "\n";
    }

    std::mutex output_mutex;
    run_in_parallel(qsketches.size(), options.threads, [&](size_t qi) {
        std::ostringstream chunk;
        for (size_t di = 0; di < dbsketches.size(); ++di) {
            const double d = jaccard_distance(qsketches[qi], dbsketches[di], options);
            chunk << qpaths[qi] << '\t' << dbpaths[di] << '\t' << d << "\n";
        }
        const std::string text = chunk.str();
        if (!text.empty()) {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << text;
        }
    });
}

void run_dist_pairlist(const std::vector<SketchPair>& pairs, const OddsketchOptions& options) {
    std::unordered_map<std::string, Sketch> cache;
    cache.reserve(pairs.size() * 2);

    for (const auto& pair : pairs) {
        if (cache.find(pair.left) == cache.end()) {
            cache.emplace(pair.left, load_sketch(pair.left));
        }
        if (cache.find(pair.right) == cache.end()) {
            cache.emplace(pair.right, load_sketch(pair.right));
        }
    }

    const auto groups = group_pairs_by_left(pairs);

    if (!pairs.empty()) {
        const auto& s0 = cache.at(pairs.front().left);
        const uint64_t nbits = s0.bit_size ? s0.bit_size : static_cast<uint64_t>(s0.words.size() * 64);
        std::cerr << "[oddsketch] dist(pairlist): nbits=" << nbits
                  << ", pairs=" << pairs.size()
                  << ", query_groups=" << groups.size()
                  << ", threads=" << normalize_thread_count(options.threads, groups.size())
                  << ", unique_sketches=" << cache.size()
                  << ", kbuckets(header)=" << s0.k_buckets
                  << ", j0(current)=" << std::fixed << std::setprecision(6) << options.j0 << "\n";
    }

    std::vector<double> results(pairs.size(), 0.0);
    run_in_parallel(groups.size(), options.threads, [&](size_t group_index) {
        const auto& group = groups[group_index];
        const auto& left = cache.at(group.left);
        for (const size_t pair_index : group.pair_indices) {
            results[pair_index] = jaccard_distance(left, cache.at(pairs[pair_index].right), options);
        }
    });

    for (size_t i = 0; i < pairs.size(); ++i) {
        std::cout << pairs[i].left << '\t' << pairs[i].right << '\t' << results[i] << "\n";
    }
}

void set_kmer(OddsketchOptions& options, const std::string& value) {
    if (value.empty()) {
        return;
    }
    options.kmer = static_cast<size_t>(std::stoul(value));
}

void set_sketch_size(OddsketchOptions& options, const std::string& value) {
    if (value.empty()) {
        return;
    }

    const size_t sketch_size = static_cast<size_t>(std::stoul(value));
    if (sketch_size == 0 || (sketch_size % 64) != 0) {
        throw std::runtime_error("sketch-size must be a positive multiple of 64");
    }
    options.sketch_size = sketch_size;
}

void set_j0(OddsketchOptions& options, const std::string& value) {
    if (value.empty()) {
        return;
    }

    const double j0 = std::stod(value);
    if (j0 <= 0.0 || j0 >= 1.0) {
        throw std::runtime_error("j0 must be in (0,1)");
    }
    options.j0 = j0;
}

void set_pos_mode(OddsketchOptions& options, const std::string& value) {
    if (value == "value") {
        options.pos_mode = PosMode::Value;
    } else if (value == "mix") {
        options.pos_mode = PosMode::Mix;
    } else if (value == "stripe") {
        options.pos_mode = PosMode::Stripe;
    } else {
        throw std::runtime_error("pos-mode must be one of: value, mix, stripe");
    }
}

void set_canonical(OddsketchOptions& options, const std::string& value) {
    if (value.empty()) {
        options.canonical = true;
        return;
    }

    std::string normalized = value;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    if (normalized == "1" || normalized == "true" || normalized == "yes" || normalized == "on") {
        options.canonical = true;
    } else if (normalized == "0" || normalized == "false" || normalized == "no" || normalized == "off") {
        options.canonical = false;
    } else {
        throw std::runtime_error("canonical must be true/false/1/0");
    }
}

void set_threads(OddsketchOptions& options, const std::string& value) {
    if (value.empty()) {
        return;
    }

    const size_t threads = static_cast<size_t>(std::stoul(value));
    if (threads == 0) {
        throw std::runtime_error("threads must be positive");
    }
    options.threads = threads;
}

bool option_requires_value(const std::string& key) {
    return key == "input-paths" ||
           key == "qlist" ||
           key == "dblist" ||
           key == "pairlist" ||
           key == "out-dir" ||
           key == "sketch-paths-out" ||
           key == "manifest" ||
           key == "kmer" ||
           key == "kmerlen" ||
           key == "sketch-size" ||
           key == "j0" ||
           key == "j-threshold" ||
           key == "pos-mode" ||
           key == "canonical" ||
           key == "threads";
}

void require_value(const std::string& key, const std::string& value) {
    if (value.empty()) {
        throw std::runtime_error("--" + key + " requires a value");
    }
}

CliArgs parse_options(int argc, char** argv) {
    CliArgs args;
    if (argc < 2) {
        return args;
    }

    args.mode = argv[1];
    if (args.mode == "--help" || args.mode == "-h") {
        return args;
    }

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg.rfind("--", 0) != 0) {
            continue;
        }

        const auto eq = arg.find('=');
        const std::string key = (eq == std::string::npos) ? arg.substr(2) : arg.substr(2, eq - 2);
        std::string value = (eq == std::string::npos) ? std::string() : arg.substr(eq + 1);

        if (eq == std::string::npos && option_requires_value(key)) {
            if (i + 1 < argc) {
                const std::string next = argv[i + 1];
                if (next.rfind("--", 0) != 0) {
                    value = next;
                    ++i;
                }
            }
            require_value(key, value);
        }

        if (key == "input-paths") {
            require_value(key, value);
            args.input_paths_path = value;
        } else if (key == "qlist") {
            require_value(key, value);
            args.qlist_path = value;
        } else if (key == "dblist") {
            require_value(key, value);
            args.dblist_path = value;
        } else if (key == "pairlist") {
            require_value(key, value);
            args.pairlist_path = value;
        } else if (key == "out-dir") {
            require_value(key, value);
            args.out_dir = value;
        } else if (key == "sketch-paths-out" || key == "manifest") {
            require_value(key, value);
            args.sketch_paths_out = value;
        } else if (key == "skip-existing") {
            args.skip_existing = true;
        } else if (key == "all-to-all" || key == "alltoall") {
            args.explicit_all_to_all = true;
        } else if (key == "bipartite") {
            args.explicit_bipartite = true;
        } else if (key == "kmer" || key == "kmerlen") {
            set_kmer(args.options, value);
        } else if (key == "sketch-size") {
            set_sketch_size(args.options, value);
        } else if (key == "j0" || key == "j-threshold") {
            set_j0(args.options, value);
        } else if (key == "pos-mode") {
            set_pos_mode(args.options, value);
        } else if (key == "canonical") {
            set_canonical(args.options, value);
        } else if (key == "threads") {
            set_threads(args.options, value);
        } else {
            throw std::runtime_error("unknown option: --" + key);
        }
    }

    return args;
}

}  // namespace

void print_usage() {
    std::cerr
        << "Usage:\n"
        << "  oddsketch sketch --input-paths genomes.list [options]\n"
        << "  oddsketch dist --all-to-all [options] < sketches.list\n"
        << "  oddsketch dist --bipartite --qlist queries.list --dblist db.list [options]\n"
        << "  oddsketch dist --pairlist pairs.tsv [options]\n"
        << "\n"
        << "Options:\n"
        << "  --input-paths=genomes.list\n"
        << "  --out-dir=sketches\n"
        << "  --sketch-paths-out=sketch_paths.txt\n"
        << "  --skip-existing\n"
        << "  --kmer=N, --kmerlen=N        k-mer length (default: 64)\n"
        << "  --sketch-size=M              sketch size in bits (default: 8192)\n"
        << "  --canonical=0|1              use canonical k-mers (default: 1)\n"
        << "  --pos-mode=value|mix|stripe  bit-position mapping mode (default: value)\n"
        << "  --j0=F                       design/reference Jaccard value (default: 0.75)\n"
        << "  --threads=N                  worker threads (default: 1)\n"
        << "\n"
        << "Compatibility:\n"
        << "  oddsketch sketch < paths.list still accepts one FASTA or FASTA.gz path per line.\n"
        << "  oddsketch dist < sketches.list still runs all-to-all.\n"
        << "  --name=value and --name value are both accepted for value options.\n";
}

int oddsketch_cli_main(int argc, char** argv) {
    CliArgs args;
    try {
        args = parse_options(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Usage error: " << e.what() << "\n";
        return 1;
    }

    if (args.mode == "--help" || args.mode == "-h") {
        print_usage();
        return 0;
    }

    if (args.mode != "sketch" && args.mode != "dist") {
        print_usage();
        return 1;
    }
    if ((args.qlist_path.empty() && !args.dblist_path.empty()) ||
        (!args.qlist_path.empty() && args.dblist_path.empty())) {
        std::cerr << "Usage error: --qlist and --dblist must be specified together\n";
        return 1;
    }

    const bool pairlist_mode = !args.pairlist_path.empty();
    const bool bipartite_mode = args.explicit_bipartite || !args.qlist_path.empty();
    const bool all_to_all_mode = args.explicit_all_to_all;
    const int dist_mode_count = (pairlist_mode ? 1 : 0) + (bipartite_mode ? 1 : 0) + (all_to_all_mode ? 1 : 0);

    if (args.mode == "dist" && dist_mode_count > 1) {
        std::cerr << "Usage error: choose only one dist mode: --all-to-all, --bipartite, or --pairlist\n";
        return 1;
    }
    if (args.mode == "dist" && !args.input_paths_path.empty()) {
        std::cerr << "Usage error: --input-paths is only valid for sketch\n";
        return 1;
    }
    if (args.mode == "dist" && (!args.out_dir.empty() || !args.sketch_paths_out.empty() || args.skip_existing)) {
        std::cerr << "Usage error: --out-dir, --sketch-paths-out, and --skip-existing are only valid for sketch\n";
        return 1;
    }
    if (args.explicit_bipartite && (args.qlist_path.empty() || args.dblist_path.empty())) {
        std::cerr << "Usage error: --bipartite requires --qlist and --dblist\n";
        return 1;
    }

    try {
        if (args.mode == "sketch") {
            const auto inputs = !args.input_paths_path.empty()
                ? read_sketch_inputs_from_input_paths_file(args.input_paths_path)
                : read_sketch_inputs_from_stdin();
            std::cerr << "[oddsketch] sketch: nbits=" << args.options.sketch_size
                      << ", kbuckets=" << compute_num_buckets(args.options)
                      << ", inputs=" << inputs.size()
                      << ", threads=" << normalize_thread_count(args.options.threads, inputs.size())
                      << ", j0=" << std::fixed << std::setprecision(6) << args.options.j0
                      << ", canonical=" << (args.options.canonical ? "true" : "false")
                      << ", pos-mode=" << pos_mode_name(args.options.pos_mode)
                      << ", out-dir=" << (args.out_dir.empty() ? "<input-dir>" : args.out_dir)
                      << ", skip-existing=" << (args.skip_existing ? "true" : "false")
                      << "\n";

            // By default sketch writes <input>.sketch; --out-dir writes one file per input basename.
            std::vector<std::string> output_paths(inputs.size());
            for (size_t i = 0; i < inputs.size(); ++i) {
                output_paths[i] = sketch_output_path(inputs[i], args);
            }
            reject_duplicate_outputs(output_paths);
            if (!args.out_dir.empty()) {
                fs::create_directories(args.out_dir);
            }

            run_in_parallel(inputs.size(), args.options.threads, [&](size_t i) {
                const auto& input = inputs[i];
                const std::string& output_path = output_paths[i];
                if (args.skip_existing && fs::exists(output_path) && fs::file_size(output_path) > 0) {
                    return;
                }
                const SketchBuildResult sketch = make_odd_sketch_from_fasta(input.sequence_path, args.options);
                write_sketch_binary(
                    sketch.words,
                    output_path,
                    static_cast<uint64_t>(args.options.sketch_size),
                    static_cast<uint64_t>(sketch.num_buckets),
                    args.options.j0
                );
            });

            if (!args.sketch_paths_out.empty()) {
                write_sketch_paths_file(args.sketch_paths_out, output_paths);
            }

            for (const auto& output_path : output_paths) {
                std::cout << output_path << "\n";
            }
        } else if (bipartite_mode) {
            run_dist_bipartite(
                read_paths_from_list_file(args.qlist_path),
                read_paths_from_list_file(args.dblist_path),
                args.options
            );
        } else if (pairlist_mode) {
            run_dist_pairlist(read_pairs_from_pairlist_file(args.pairlist_path), args.options);
        } else {
            run_dist_all_vs_all(read_paths_from_stdin(), args.options);
        }
    } catch (const std::exception& e) {
        std::cerr << "[oddsketch] error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
