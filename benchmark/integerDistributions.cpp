#include <ctime>
#include <chrono>
#include <unordered_set>
#include <XorShift64.h>
#include <tlx/cmdline_parser.hpp>
#include <DirectRankStoring.hpp>
#include <ExampleCompetitorMmphfs.hpp>

template <typename MMphf>
void doTest(std::vector<uint64_t> &inputData) {
    std::cout<<"\r\033[K"<<"Generating "<<MMphf::name()<<std::flush;
    std::chrono::steady_clock::time_point beginConstr = std::chrono::steady_clock::now();
    MMphf mmphf(inputData);
    std::chrono::steady_clock::time_point endConstr = std::chrono::steady_clock::now();

    std::cout<<"\r\033[K"<<"Verifying "<<MMphf::name()<<std::flush;
    for (size_t i = 0; i < inputData.size(); i++) {
        if (mmphf(inputData.at(i)) != i) {
            std::cerr<<"Error verifying"<<std::endl;
        }
    }

    std::cout<<"\r\033[K"<<"Benchmarking "<<MMphf::name()<<std::flush;
    size_t numQueries = 5e6;
    size_t N = inputData.size();
    util::XorShift64 prng(time(nullptr));
    uint64_t h = prng();
    std::chrono::steady_clock::time_point beginQuery = std::chrono::steady_clock::now();
    for (size_t i = 0; i < numQueries; i++) {
        h ^= mmphf(inputData.at(h % N)) ^ prng();
    }
    std::chrono::steady_clock::time_point endQuery = std::chrono::steady_clock::now();
    std::cout<<h<<"\r\033[K"<<"Done."<<std::endl;

    double spacePerElement = (double) mmphf.spaceBits() / N;

    // This format can be converted to TikZ/PGF plots using https://github.com/bingmann/sqlplot-tools
    std::cout << "RESULT"
              << " competitor=" << MMphf::name()
              << " N=" << N
              << " numQueries=" << numQueries
              << " spacePerObject=" << spacePerElement
              << " constructionMs=" << std::chrono::duration_cast<std::chrono::milliseconds>(endConstr - beginConstr).count()
              << " queryMs=" << std::chrono::duration_cast<std::chrono::milliseconds>(endQuery - beginQuery).count()
              << std::endl;
    std::cout << std::endl;
}

template <float rankMin, float rankDelta, float rankMax>
void dispatchDirectRankStoring(std::vector<uint64_t> &inputData) {
    doTest<DirectRankStoringMmphf<LinearBucketMapper<rankMin>>>(inputData);
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
    doTest<DirectRankStoringMmphf<LinearBucketMapper<1.0f>>>(inputData);
    doTest<DirectRankStoringMmphf<LinearBucketMapper<1.125f>>>(inputData);
    doTest<DirectRankStoringMmphf<PgmBucketMapper<1.0f>>>(inputData);

    return 0;
}
