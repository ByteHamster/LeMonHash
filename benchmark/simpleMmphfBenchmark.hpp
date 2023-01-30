#pragma once

#include <ctime>
#include <chrono>
#include <XorShift64.h>
#include <vector>
#include <iostream>
#include <malloc.h>

template <typename MMphf, typename key_t>
void simpleMmphfBenchmark(std::vector<key_t> &inputData) {
    std::cout<<"\r\033[K"<<"Generating "<<MMphf::name()<<std::flush;

    size_t spaceBefore = mallinfo2().uordblks;
    std::chrono::steady_clock::time_point beginConstr = std::chrono::steady_clock::now();
    MMphf mmphf(inputData);
    std::chrono::steady_clock::time_point endConstr = std::chrono::steady_clock::now();
    size_t spaceAfter = mallinfo2().uordblks + sizeof(MMphf);

    std::cout<<"\r\033[K"<<"Verifying "<<MMphf::name()<<std::flush;
    for (size_t i = 0; i < inputData.size(); i++) {
        if (mmphf(inputData.at(i)) != i) {
            std::cerr<<std::endl<<std::endl<<"Error verifying key "<<i<<std::endl;
            return;
        }
    }

    std::cout<<"\r\033[K"<<"Benchmarking "<<MMphf::name()<<std::flush;
    size_t numQueries = 1e6;
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
    // Note that this is a bit larger than what is calculated by the data structure itself.
    // The discrepancy decreases the larger we chose N.
    // My guess would be that this is because containers allocate more space than necessary to store their content.
    std::cout<<"Space sanity check: mallinfo2 returns "<<8.0*(spaceAfter-spaceBefore)/N<<"/object"<<std::endl;
    std::cout << std::endl;
}
