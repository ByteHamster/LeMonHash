#include <unordered_set>
#include <tlx/cmdline_parser.hpp>
#include <DirectRankStoring.hpp>
#include "simpleMmphfBenchmark.hpp"
#include <XorShift64.h>
#include "ChunkedStringMmphf.hpp"
#include "simpleMmphfBenchmark.hpp"
#include "RecursiveDirectRankStoring.hpp"

std::string stringOfLength(size_t length, util::XorShift64 &prng) {
    std::string key(length, 'x');
    size_t i = 0;
    while (i < (length & ~0xf)) {
        *((uint64_t*) (key.data() + i)) = prng();
        i += 8;
    }
    while (i < length) {
        key[i] = (char) prng();
        i += 1;
    }
    return key;
}

std::vector<std::string> randomUniformStrings(size_t n, size_t length, size_t commonPrefixLength) {
    std::unordered_set<std::string> keys;
    keys.reserve(n);
    size_t seed = time(nullptr);
    std::cout<<"Seed: "<<seed<<std::endl;
    util::XorShift64 prng(seed);
    std::string prefix = stringOfLength(commonPrefixLength, prng);
    while (keys.size() < n) {
        std::string key = prefix + stringOfLength(length, prng);
        keys.insert(key);
    }
    std::vector<std::string> dataset;
    dataset.insert(dataset.end(), keys.begin(), keys.end());
    std::sort(dataset.begin(), dataset.end());
    return dataset;
}

std::vector<std::string> loadFile(std::string &filename, size_t maxStrings) {
    std::vector<std::string> inputData;
    inputData.reserve(maxStrings);
    std::ifstream stream(filename);
    const int MAX_LENGTH = 524288;
    char* line = new char[MAX_LENGTH];
    while (stream.getline(line, MAX_LENGTH)) {
        inputData.emplace_back(line);
        if (inputData.size() >= maxStrings) {
            break;
        }
    }
    delete[] line;
    std::cout<<"Sorting input data"<<std::endl;
    std::sort(inputData.begin(), inputData.end());
    std::cout<<"Loaded "<<inputData.size()<<" strings"<<std::endl;
    return inputData;
}

int main(int argc, char** argv) {
    size_t N = 1e6;
    std::string filename;
    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "num_keys", N, "Number of keys to generate");
    cmd.add_string('f', "filename", filename, "File with input data");
    if (!cmd.process(argc, argv)) {
        return 1;
    }

    std::vector<std::string> inputData;
    if (filename.empty()) {
        std::cout<<"Generating input data"<<std::endl;
        inputData = randomUniformStrings(N, 10, 4);
        /*inputData.emplace_back("aaaaaaaaaaaaaaaa");
        inputData.emplace_back("aaaaaaaabbbbbbbbcccccccc");
        inputData.emplace_back("bbbbbbbbaaaaaaaa");*/
    } else {
        std::cout<<"Loading file "<<filename<<std::endl;
        inputData = loadFile(filename, N);
    }

    struct FullChunking {
        static ChunkingStrategy *createLayer(size_t maxLCP, size_t layer) {
            return new FullChunkingStrategy(maxLCP, 8);
        }
    };
    simpleMmphfBenchmark<ChunkedStringMmphf<FullChunking>>(inputData);

    //struct GreedyChunking {
    //    static ChunkingStrategy *createLayer(size_t maxLCP, size_t layer) {
    //        return new GreedyChunkingStrategy(maxLCP, 8);
    //    }
    //};
    //simpleMmphfBenchmark<ChunkedStringMmphf<GreedyChunking>>(inputData); // Very slow

    struct SeparateChunking {
        static ChunkingStrategy *createLayer(size_t maxLCP, size_t layer) {
            return new SeparateChunkingStrategy(maxLCP, 8);
        }
    };
    simpleMmphfBenchmark<ChunkedStringMmphf<SeparateChunking>>(inputData);

    struct LayeredChunking {
        static ChunkingStrategy *createLayer(size_t maxLCP, size_t layer) {
            if (layer == 0) {
                return new FullChunkingStrategy(maxLCP, 3); // Alphabet reduction
            } else {
                return new FullChunkingStrategy(maxLCP, 8);
            }
        }
    };
    simpleMmphfBenchmark<ChunkedStringMmphf<LayeredChunking>>(inputData);

    struct LayeredChunkingForUrls {
        static ChunkingStrategy *createLayer(size_t maxLCP, size_t layer) {
            if (layer == 0) {
                return new FullChunkingStrategy(maxLCP, 3); // Alphabet reduction
            } else if (layer == 1) {
                return new BackChunkingStrategy(maxLCP, 8); // Make long strings shorter
            } else {
                return new FullChunkingStrategy(maxLCP, 8);
            }
        }
    };
    simpleMmphfBenchmark<ChunkedStringMmphf<LayeredChunkingForUrls>>(inputData);

    simpleMmphfBenchmark<RecursiveDirectRankStoringMmphf>(inputData);

    return 0;
}
