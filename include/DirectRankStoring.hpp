#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include <EliasFano.h>
#include <Function.h>
#include "MultiRetrievalDataStructure.hpp"
#include "bucket_mapping/LinearBucketMapper.hpp"
#include "bucket_mapping/PgmBucketMapper.hpp"
#include "bucket_mapping/SuccinctPgmBucketMapper.hpp"

/**
 * Monotone Minimal Perfect Hash Function (MMPHF) using the Direct Rank Storing (DRS) technique.
 * Each object is mapped to a bucket (using different techniques).
 * Within the buckets, a retrieval data structure explicitly stores the ranks of all objects.
 * Given that we expect only one object per bucket, we can usually use a 1-bit retrieval data structure.
 * The bucket sizes are stored with Elias-Fano.
 */
template <typename BucketMapper>
class DirectRankStoringMmphf {
    public:
        size_t N;
    private:
        BucketMapper bucketMapper;
        MultiRetrievalDataStructure retrievalDataStructure;
        util::EliasFano<util::floorlog2(std::max(1.0f, BucketMapper::elementsPerBucket()))> bucketSizePrefix;
    public:
        static std::string name() {
            return "DirectRankStoringMmphf bucketMapper=" + BucketMapper::name();
        }

        explicit DirectRankStoringMmphf(std::vector<uint64_t> &data)
                : N(data.size()),
                  bucketMapper(data.begin(), data.end()),
                  bucketSizePrefix(bucketMapper.numBuckets() + 1, data.size() + 1) {

            std::vector<uint64_t> currentBucket;
            size_t prev_bucket = 0;
            size_t bucketSizePrefixTemp = 0;

            auto constructBucket = [&] {
                size_t bucketSize = currentBucket.size();
                if (bucketSize > 1) {
                    for (size_t j = 0; j < bucketSize; j++) {
                        retrievalDataStructure.addInput(bucketSize, currentBucket.at(j), j);
                    }
                }
                bucketSizePrefix.push_back(bucketSizePrefixTemp);
                bucketSizePrefixTemp += bucketSize;
                currentBucket.clear();
            };

            bucketMapper.bucketOf(data.begin(), data.end(), [&] (auto it, auto bucket) {
                if (it != data.begin()) {
                    if (*it <= *std::prev(it))
                        throw std::invalid_argument("Data not sorted or duplicates found");
                    if (bucket < prev_bucket)
                        throw std::runtime_error("Non-monotonic bucket mapper");
                }
                while (bucket != prev_bucket) {
                    constructBucket();
                    prev_bucket++;
                }
                currentBucket.push_back(*it);
            });

            while (prev_bucket < bucketMapper.numBuckets() + 1) {
                constructBucket();
                prev_bucket++;
            }
            bucketSizePrefix.buildRankSelect();

            retrievalDataStructure.build();
        }

        size_t operator()(uint64_t key) {
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
                std::cout << "BucketMapper: " << (8.0 * bucketMapper.size() / N)
                          << " (" << bucketMapper.info() << ")" << std::endl;
            }
            size_t bytes = bucketMapper.size() + sizeof(*this) + bucketSizePrefix.space();
            return 8 * bytes + retrievalDataStructure.spaceBits(0);
        }
};
