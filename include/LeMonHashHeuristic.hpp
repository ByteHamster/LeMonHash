#pragma once

#include <vector>
#include <cstdint>
#include <algorithm>
#include <iostream>
#include <bytehamster/util/Function.h>
#include <MultiRetrievalDataStructure.hpp>
#include "support/sequence/BucketOffsets.hpp"
#include "bucket_mapping/LinearBucketMapper.hpp"
#include "bucket_mapping/PGMBucketMapper.hpp"
#include "bucket_mapping/SuccinctPGMBucketMapper.hpp"

namespace lemonhash {
/**
 * Learned Monotone Minimal Perfect Hash Function (MMPHF).
 * Each object is mapped to a bucket using a heuristic segmented linear mapper.
 * Within the buckets, a retrieval data structure explicitly stores the ranks of all objects.
 * The prefix sum of bucket sizes is stored with Elias-Fano.
 */
template <size_t keysPerSegment, size_t retrievalCoeffBits = 64>
class LeMonHashHeuristic {
    public:
        const size_t N;
    private:
        const size_t segments;
        EliasFanoM eliasFano;
        BucketOffsets bucketOffsets;
        MultiRetrievalDataStructure<retrievalCoeffBits> retrievalDataStructure;
    public:
        static std::string name() {
            return std::string("LeMonHash")
                   + " bucketMapper=HeuristicSegmented"
                   + " retrievalCoeffBits=" + std::to_string(retrievalCoeffBits);
        }

        explicit LeMonHashHeuristic(const std::vector<uint64_t> &data)
                : N(data.size()),
                  segments((N + keysPerSegment - 1) / keysPerSegment),
                  eliasFano(segments + 2, data.back()),
                  bucketOffsets(N, data.size()) {
            size_t prevBucket = 0;
            size_t currentBucketBegin = 0;
            size_t currentBucketEnd = 0;

            auto constructBucket = [&] {
                size_t bucketSize = currentBucketEnd - currentBucketBegin;
                if (bucketSize > 1) {
                    for (size_t j = 0; j < bucketSize; j++) {
                        retrievalDataStructure.addInput(bucketSize, data.at(currentBucketBegin + j), j);
                    }
                }
                bucketOffsets.push(bucketSize);
            };

            for (size_t segment = 0; segment < segments; segment++) {
                size_t firstIndexOfSegment = segment * keysPerSegment;
                uint64_t segmentStart = data.at(firstIndexOfSegment);
                uint64_t segmentEnd = data.at(std::min(firstIndexOfSegment + keysPerSegment, N - 1));
                eliasFano.push_back(segmentStart);
                for (size_t i = 0; i < keysPerSegment && firstIndexOfSegment + i < N; i++) {
                    size_t keyIndex = firstIndexOfSegment + i;
                    size_t bucket = bucketOf(data.at(keyIndex), segment, segmentStart, segmentEnd);
                    if (keyIndex != 0) { [[likely]]
                        if (data.at(keyIndex) <= data.at(keyIndex - 1)) [[unlikely]]
                            throw std::invalid_argument("Data not sorted or duplicates found");
                        if (bucket < prevBucket) [[unlikely]]
                            throw std::runtime_error("Non-monotonic bucket mapper");
                    }
                    while (prevBucket < bucket) {
                        currentBucketEnd = keyIndex;
                        constructBucket();
                        currentBucketBegin = currentBucketEnd;
                        prevBucket++;
                    }
                }
            }
            eliasFano.push_back(data.back());
            currentBucketEnd = N;
            constructBucket();

            bucketOffsets.done();
            eliasFano.buildRankSelect();
            retrievalDataStructure.build();
        }

        [[nodiscard]] size_t bucketOf(uint64_t key, size_t segment, uint64_t offset, uint64_t nextOffset) const {
            double slope = (float) (nextOffset - offset) / keysPerSegment;
            size_t rankInSegment = std::min<size_t>((double) (key - offset) / slope, keysPerSegment);
            size_t estimatedRank = segment * keysPerSegment + rankInSegment;
            return std::min(N - 1, estimatedRank);
        }

        size_t operator()(const uint64_t key) {
            auto segmentPtr = eliasFano.predecessorPosition(key);
            size_t segment = segmentPtr.index();
            uint64_t offset = *segmentPtr;
            uint64_t nextOffset = UINT64_MAX;
            if (segment < eliasFano.size() - 1) {
                ++segmentPtr;
                nextOffset = *segmentPtr;
            }
            size_t bucket = bucketOf(key, segment, offset, nextOffset);

            auto [bucketOffset, bucketSize] = bucketOffsets.at(bucket);
            if (bucketSize <= 1) {
                return bucketOffset;
            } else {
                return bucketOffset + retrievalDataStructure.query(bucketSize, key);
            }
        }

        size_t spaceBits(bool print = true) {
            if (print) {
                std::cout << "Retrieval:      " << (double(retrievalDataStructure.spaceBits(0)) / N) << std::endl;
                std::cout << "Segments:       " << (8.0 * eliasFano.space() / N) << std::endl;
                std::cout << "Bucket offsets: " << (8.0 * bucketOffsets.space() / N) << std::endl;
            }
            size_t bytes = eliasFano.space() + sizeof(*this) + bucketOffsets.space();
            return 8 * bytes + retrievalDataStructure.spaceBits(0);
        }
};
} // namespace lemonhash
