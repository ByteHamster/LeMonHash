#pragma once

#include <ctime>
#include <chrono>
#include <XorShift64.h>
#include <vector>
#include <iostream>

template <typename MMphf, typename key_t>
void simpleMmphfBenchmark(std::vector<key_t> &inputData) {
    std::cout<<"\r\033[K"<<"Generating "<<MMphf::name()<<std::flush;
    std::chrono::steady_clock::time_point beginConstr = std::chrono::steady_clock::now();
    MMphf mmphf(inputData);
    std::chrono::steady_clock::time_point endConstr = std::chrono::steady_clock::now();

    std::cout<<"\r\033[K"<<"Verifying "<<MMphf::name()<<std::flush;
    for (size_t i = 0; i < inputData.size(); i++) {
        if (mmphf(inputData.at(i)) != i) {
            std::cerr<<std::endl<<std::endl<<"Error verifying key "<<i<<std::endl;
            return;
        }
    }

    std::cout<<"\r\033[K"<<"Benchmarking "<<MMphf::name()<<std::flush;
    size_t numQueries = 1e5;
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
