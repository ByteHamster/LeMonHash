#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include <EliasFano.h>
#include <Function.h>
#include "support/MultiRetrievalDataStructure.hpp"
#include "bucket_mapping/LinearBucketMapper.hpp"
#include "bucket_mapping/PGMBucketMapper.hpp"
#include "bucket_mapping/SuccinctPGMBucketMapper.hpp"

/**
 * Learned Monotone Minimal Perfect Hash Function (MMPHF).
 * Each object is mapped to a bucket using the PGM index.
 * Within the buckets, a retrieval data structure explicitly stores the ranks of all objects.
 * The prefix sum of bucket sizes is stored with Elias-Fano.
 */
template <typename BucketMapper = SuccinctPGMBucketMapper, size_t retrievalCoeffBits = 64>
class LeMonHash {
    public:
        size_t N;
    private:
        BucketMapper bucketMapper;
        MultiRetrievalDataStructure<retrievalCoeffBits> retrievalDataStructure;
        util::EliasFano<util::floorlog2(std::max(1.0f, BucketMapper::elementsPerBucket()))> bucketSizePrefix;
    public:
        static std::string name() {
            return std::string("LeMonHash")
                   + " bucketMapper=" + BucketMapper::name()
                   + " retrievalCoeffBits=" + std::to_string(retrievalCoeffBits);
        }

        explicit LeMonHash(const std::vector<uint64_t> &data)
                : N(data.size()),
                  bucketMapper(data.begin(), data.end()),
                  bucketSizePrefix(bucketMapper.numBuckets() + 1, data.size() + 1) {
            size_t prevBucket = 0;
            size_t bucketSizePrefixTemp = 0;
            auto currentBucketBegin = data.begin();
            auto currentBucketEnd = data.begin();

            auto constructBucket = [&] {
                size_t bucketSize = currentBucketEnd - currentBucketBegin;
                if (bucketSize > 1) {
                    for (size_t j = 0; j < bucketSize; j++) {
                        retrievalDataStructure.addInput(bucketSize, currentBucketBegin[j], j);
                    }
                }
                bucketSizePrefix.push_back(bucketSizePrefixTemp);
                bucketSizePrefixTemp += bucketSize;
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
            bucketSizePrefix.push_back(bucketSizePrefixTemp);

            bucketSizePrefix.buildRankSelect();

            retrievalDataStructure.build();
        }

        size_t operator()(const uint64_t key) {
            auto ptr = bucketSizePrefix.at(bucketMapper.bucketOf(key));
            size_t bucketOffset = *ptr;
            ++ptr;
            size_t nextBucketOffset = *ptr;
            size_t bucketSize = nextBucketOffset - bucketOffset;
            if (bucketSize <= 1) {
                return bucketOffset;
            } else {
                return bucketOffset + retrievalDataStructure.query(bucketSize, key);
            }
        }

        size_t spaceBits(bool print = true) {
            if (print) {
                std::cout << "EliasFano:    " << (8.0 * bucketSizePrefix.space() / N) << std::endl;
                std::cout << "BucketMapper: " << (8.0 * bucketMapper.size() / N) << " (" << bucketMapper.info() << ")" << std::endl;
                std::cout << "Retrieval:    " << (double(retrievalDataStructure.spaceBits(0)) / N) << std::endl;
            }
            size_t bytes = bucketMapper.size() + sizeof(*this) + bucketSizePrefix.space();
            return 8 * bytes + retrievalDataStructure.spaceBits(0);
        }
};