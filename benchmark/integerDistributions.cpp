#include <unordered_set>
#include <tlx/cmdline_parser.hpp>
#include <DirectRankStoring.hpp>
#include <ExampleCompetitorMmphfs.hpp>
#include "simpleMmphfBenchmark.hpp"

template <float rankMin, float rankDelta, float rankMax>
void dispatchDirectRankStoring(std::vector<uint64_t> &inputData) {
    simpleMmphfBenchmark<DirectRankStoringMmphf<LinearBucketMapper<rankMin>>>(inputData);
    if constexpr (rankMin < rankMax) {
        dispatchDirectRankStoring<rankMin+rankDelta, rankDelta, rankMax>(inputData);
    }
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
    size_t N = 1e4;
    std::string type = "uniform";

    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "num_keys", N, "Number of keys to generate");
    cmd.add_string('t', "type", type, "Type of data to generate (uniform or pareto)");
    if (!cmd.process(argc, argv)) {
        return 1;
    }

    std::cout<<"Generating input data"<<std::endl;
    std::vector<uint64_t> inputData;
    if (type == "pareto") {
        inputData = randomPareto(1e4, 1);
    } else if (type == "uniform") {
        inputData = randomUniform(1e6, UINT64_MAX - 1);
    } else {
        cmd.print_usage();
        return 1;
    }

    //doTest<BinarySearchMmphf>(inputData);
    //doTest<DirectRetrievalMmphf>(inputData);
    //dispatchDirectRankStoring<1.0f, 0.01f, 1.2f>(inputData);
    simpleMmphfBenchmark<DirectRankStoringMmphf<LinearBucketMapper<1.0f>>>(inputData);
    simpleMmphfBenchmark<DirectRankStoringMmphf<LinearBucketMapper<1.125f>>>(inputData);
    simpleMmphfBenchmark<DirectRankStoringMmphf<PgmBucketMapper>>(inputData);
    simpleMmphfBenchmark<DirectRankStoringMmphf<SuccinctPgmBucketMapper>>(inputData);

    return 0;
}
