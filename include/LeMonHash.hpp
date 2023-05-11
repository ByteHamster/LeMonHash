#pragma once

#include <vector>
#include <cstdint>
#include <algorithm>
#include <iostream>
#include <Function.h>
#include "support/sequence/BucketOffsets.hpp"
#include "support/MultiRetrievalDataStructure.hpp"
#include "bucket_mapping/LinearBucketMapper.hpp"
#include "bucket_mapping/PGMBucketMapper.hpp"
#include "bucket_mapping/SuccinctPGMBucketMapper.hpp"

namespace lemonhash {
/**
 * Learned Monotone Minimal Perfect Hash Function (MMPHF).
 * Each object is mapped to a bucket using the PGM index.
 * Within the buckets, a retrieval data structure explicitly stores the ranks of all objects.
 * The prefix sum of bucket sizes is stored with Elias-Fano.
 */
template <typename BucketMapper = SuccinctPGMBucketMapper<>, size_t retrievalCoeffBits = 64>
class LeMonHash {
    public:
        size_t N;
    private:
        BucketMapper bucketMapper;
        MultiRetrievalDataStructure<retrievalCoeffBits> retrievalDataStructure;
        BucketOffsets bucketOffsets;
    public:
        static std::string name() {
            return std::string("LeMonHash")
                   + " bucketMapper=" + BucketMapper::name()
                   + " retrievalCoeffBits=" + std::to_string(retrievalCoeffBits);
        }

        explicit LeMonHash(const std::vector<uint64_t> &data)
                : N(data.size()),
                  bucketMapper(data.begin(), data.end()),
                  bucketOffsets(bucketMapper.numBuckets(), data.size()) {
            size_t prevBucket = 0;
            auto currentBucketBegin = data.begin();
            auto currentBucketEnd = data.begin();

            auto constructBucket = [&] {
                size_t bucketSize = currentBucketEnd - currentBucketBegin;
                if (bucketSize > 1) {
                    for (size_t j = 0; j < bucketSize; j++) {
                        retrievalDataStructure.addInput(bucketSize, currentBucketBegin[j], j);
                    }
                }
                bucketOffsets.push(bucketSize);
            };

            bucketMapper.bucketOf(data.begin(), data.end(), [&] (auto it, size_t bucket) {
                if (it != data.begin()) { [[likely]]
                    if (*it <= *std::prev(it)) [[unlikely]]
                        throw std::invalid_argument("Data not sorted or duplicates found");
                    if (bucket < prevBucket) [[unlikely]]
                        throw std::runtime_error("Non-monotonic bucket mapper");
                }
                while (prevBucket < bucket) {
                    currentBucketEnd = it;
                    constructBucket();
                    currentBucketBegin = currentBucketEnd;
                    prevBucket++;
                }
            });

            currentBucketEnd = data.end();
            constructBucket();
            bucketOffsets.done();

            retrievalDataStructure.build();
        }

        size_t operator()(const uint64_t key) {
            auto [bucketOffset, bucketSize] = bucketOffsets.at(bucketMapper.bucketOf(key));
            if (bucketSize <= 1) {
                return bucketOffset;
            } else {
                return bucketOffset + retrievalDataStructure.query(bucketSize, key);
            }
        }

        size_t spaceBits(bool print = true) {
            if (print) {
                std::cout << "Retrieval:      " << (double(retrievalDataStructure.spaceBits(0)) / N) << std::endl;
                std::cout << "Bucket mapper:  " << (8.0 * bucketMapper.space() / N) << " (" << bucketMapper.info() << ")" << std::endl;
                std::cout << "Bucket offsets: " << (8.0 * bucketOffsets.space() / N) << std::endl;
            }
            size_t bytes = bucketMapper.space() + sizeof(*this) + bucketOffsets.space();
            return 8 * bytes + retrievalDataStructure.spaceBits(0);
        }
};
} // namespace lemonhash
