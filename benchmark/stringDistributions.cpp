#include <unordered_set>
#include <tlx/cmdline_parser.hpp>
#include <DirectRankStoring.hpp>
#include "simpleMmphfBenchmark.hpp"
#include <XorShift64.h>
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
    if (!stream) throw std::system_error(errno, std::system_category(), "failed to open " + filename);
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
    size_t N = std::numeric_limits<size_t>::max();
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
    } else {
        std::cout<<"Loading file "<<filename<<std::endl;
        inputData = loadFile(filename, N);
    }
    if (inputData.empty()) {
        std::cout<<"Empty input data. Used invalid file?"<<std::endl;
        return 1;
    }

    size_t positionOfSlash = filename.find_last_of('/');
    std::string baseFilename = positionOfSlash == std::string::npos ? filename : filename.substr(positionOfSlash + 1);
    simpleMmphfBenchmark<RecursiveDirectRankStoringMmphf>(inputData, baseFilename);
    simpleMmphfBenchmark<RecursiveDirectRankStoringV2Mmphf>(inputData, baseFilename);

    return 0;
}
