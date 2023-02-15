#include <unordered_set>
#include <tlx/cmdline_parser.hpp>
#include <DirectRankStoring.hpp>
#include <RecursiveDirectRankStoring.hpp>
#include "simpleMmphfBenchmark.hpp"

std::vector<uint64_t> loadIntegerFile(std::string &filename, size_t maxN) {
    std::cout<<"Loading input file"<<std::endl;
    std::ifstream fileIn(filename, std::ios::in | std::ios::binary);
    if (!fileIn) throw std::system_error(errno, std::system_category(), "failed to open " + filename);
    size_t size = 0;
    fileIn.read(reinterpret_cast<char *>(&size), sizeof(size_t));
    size = std::min(size, maxN);
    std::vector<uint64_t> inputData(size);
    fileIn.read(reinterpret_cast<char *>(inputData.data()), size * sizeof(uint64_t));
    fileIn.close();
    std::cout<<"Sorting input data"<<std::endl;
    std::sort(inputData.begin(), inputData.end());
    std::cout<<"Loaded "<<inputData.size()<<" integers"<<std::endl;
    return inputData;
}

std::vector<uint64_t> randomUniform(size_t n, uint64_t u) {
    std::vector<uint64_t> dataset;
    dataset.reserve(n);
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint64_t> dist(0, u);
    {
        std::unordered_set<uint64_t> keys;
        keys.reserve(n);
        while (dataset.size() < n) {
            auto key = dist(rng);
            if (keys.insert(key).second)
                dataset.push_back(key);
        }
    }
    std::sort(dataset.begin(), dataset.end());
    if (dataset.back() != u)
        dataset.back() = u;
    return dataset;
}

std::vector<uint64_t> randomPareto(size_t n, double shape = 1.1) {
    std::vector<uint64_t> dataset;
    dataset.reserve(n);
    std::mt19937_64 rng(42);
    std::exponential_distribution<double> dist(shape);
    {
        std::unordered_set<uint64_t> keys;
        keys.reserve(n);
        while (dataset.size() < n) {
            auto key = uint64_t(std::exp(dist(rng)));
            if (keys.insert(key).second)
                dataset.push_back(key);
        }
    }
    std::sort(dataset.begin(), dataset.end());
    return dataset;
}

int main(int argc, char** argv) {
    size_t N = std::numeric_limits<size_t>::max();
    std::string filename;
    std::string type = "uniform";
    std::string datasetName = "";

    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "num_keys", N, "Number of keys to generate");
    cmd.add_string('t', "type", type, "Type of data to generate (uniform or pareto)");
    cmd.add_string('f', "filename", filename, "Input data set to load. First 64 bits must be length, then all following words are integers");
    if (!cmd.process(argc, argv)) {
        return 1;
    }

    std::cout<<"Generating input data"<<std::endl;
    std::vector<uint64_t> inputData;
    if (!filename.empty()) {
        inputData = loadIntegerFile(filename, N);
        size_t positionOfSlash = filename.find_last_of('/');
        datasetName = positionOfSlash == std::string::npos ? filename : filename.substr(positionOfSlash + 1);
    } else if (type == "pareto") {
        inputData = randomPareto(N, 1);
        datasetName = "pareto";
    } else if (type == "uniform") {
        inputData = randomUniform(N, UINT64_MAX - 1);
        datasetName = "uniform";
    } else {
        cmd.print_usage();
        return 1;
    }

    //doTest<BinarySearchMmphf>(inputData);
    //doTest<DirectRetrievalMmphf>(inputData);
    //simpleMmphfBenchmark<DirectRankStoringMmphf<LinearBucketMapper<1.0f>>>(inputData);
    //simpleMmphfBenchmark<DirectRankStoringMmphf<LinearBucketMapper<1.125f>>>(inputData);
    //simpleMmphfBenchmark<DirectRankStoringMmphf<PgmBucketMapper>>(inputData);
    simpleMmphfBenchmark<DirectRankStoringMmphf<SuccinctPgmBucketMapper>>(inputData, datasetName);

    std::vector<std::string> inputAsString;
    inputAsString.reserve(inputData.size());
    for (uint64_t x : inputData) {
        uint64_t swapped = __builtin_bswap64(x);
        inputAsString.emplace_back((char*) &swapped, sizeof(uint64_t));
    }
    simpleMmphfBenchmark<RecursiveDirectRankStoringMmphf>(inputAsString, datasetName);

    return 0;
}
