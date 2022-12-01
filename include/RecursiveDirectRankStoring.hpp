#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include "chunking_strategy/FullChunkingStrategy.h"
#include "chunking_strategy/GreedyChunkingStrategy.h"
#include "chunking_strategy/SeparateChunkingStrategy.h"
#include "bucket_mapping/SuccinctPgmBucketMapper.hpp"
#include "DirectRankStoring.hpp"
#include "bucket_mapping/support/EliasFanoModified.hpp"
#include <MurmurHash64.h>

/**
 * Same idea as DirectRankStoring, but for strings.
 * Because we cannot use PGM to directly map full strings to buckets,
 * this only maps the first 8 bytes to buckets.
 * Some buckets might therefore be very full.
 * If a bucket is larger than DIRECT_RANK_STORING_THRESHOLD, we recurse
 * the data structure inside that bucket. For recursion, we cut off the
 * characters that all strings of that bucket have in common.
 */
class RecursiveDirectRankStoringMmphf {
    private:
        static constexpr size_t DIRECT_RANK_STORING_THRESHOLD = 4096;

        // As in normal DirectRankStoring:
        SuccinctPgmBucketMapper *bucketMapper = nullptr;
        util::EliasFanoM *bucketSizePrefix = nullptr;
        MultiRetrievalDataStructure retrieval;

        // If a bucket is very full, recurse
        pasta::BitVector recurseBucket;
        pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES> *recurseBucketRank = nullptr;
        std::vector<RecursiveDirectRankStoringMmphf *> children;

        // Cut off characters that are the same in all input strings. Useful especially for recursion.
        size_t minLCP = 9999999;
    public:
        explicit RecursiveDirectRankStoringMmphf(std::vector<std::string> strings) {
            std::unordered_set<uint64_t> chunks;
            for (size_t i = 0; i < strings.size(); i++) {
                if (i > 0) {
                    minLCP = std::min(minLCP, LCP(strings.at(i), strings.at(i-1)));
                }
            }

            // Trim strings and extract chunks
            for (auto & string : strings) {
                string = string.substr(minLCP);
                uint64_t chunk = readChunk(string.c_str(), string.length(), 8);
                chunks.insert(chunk);
            }

            std::vector<uint64_t> sortedChunks;
            sortedChunks.insert(sortedChunks.end(), chunks.begin(), chunks.end());
            std::sort(sortedChunks.begin(), sortedChunks.end());
            assert(!sortedChunks.empty());

            bucketMapper = new SuccinctPgmBucketMapper(sortedChunks.begin(), sortedChunks.end());

            std::vector<std::vector<std::string>> buckets;
            buckets.resize(bucketMapper->numBuckets);
            // This can later be done as a linear scan, just like in DirectRankStoring
            for (std::string &string : strings) {
                uint64_t chunk = readChunk(string.c_str(), string.length(), 8);
                size_t bucket = bucketMapper->bucketOf(chunk);
                buckets.at(bucket).push_back(string);
            }

            bucketSizePrefix = new util::EliasFanoM(buckets.size() + 1, strings.size() + 1);
            recurseBucket.resize(buckets.size());

            size_t bucketSizePrefixTemp = 0;
            for (size_t i = 0; i < buckets.size(); i++) {
                bucketSizePrefix->push_back(bucketSizePrefixTemp);
                bucketSizePrefixTemp += buckets.at(i).size();

                if (buckets.at(i).size() < DIRECT_RANK_STORING_THRESHOLD) {
                    recurseBucket[i] = false;

                    // Perform direct rank storing
                    for (size_t k = 0; k < buckets.at(i).size(); k++) {
                        retrieval.addInput(buckets.at(i).size(), util::MurmurHash64(buckets.at(i).at(k)), k);
                    }
                } else {
                    recurseBucket[i] = true;
                    children.push_back(new RecursiveDirectRankStoringMmphf(buckets.at(i)));
                }
            }
            bucketSizePrefix->push_back(bucketSizePrefixTemp);
            bucketSizePrefix->buildRankSelect();
            recurseBucketRank = new pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES>(recurseBucket);
            retrieval.build();
        }

        ~RecursiveDirectRankStoringMmphf() {
            delete bucketMapper;
            delete bucketSizePrefix;
            delete recurseBucketRank;
            for (auto *child : children) {
                delete child;
            }
        }

        static std::string name() {
            return "RecursiveDirectRankStoringMmphf";
        }

        size_t spaceBits() {
            size_t bits = 8 * bucketMapper->size()
                    + retrieval.spaceTheory
                    + 8 * bucketSizePrefix->space()
                    + 8 * sizeof(*this)
                    + recurseBucket.size() + 8 * recurseBucketRank->space_usage()
                    + children.size() * sizeof(uint64_t);
            for (auto *child : children) {
                bits += child->spaceBits();
            }
            return bits;
        }

        uint64_t operator ()(std::string string) {
            string = string.substr(minLCP);
            uint64_t chunk = readChunk(string.c_str(), string.length(), 8);
            size_t bucket = bucketMapper->bucketOf(chunk);
            auto bucketOffsetPtr = bucketSizePrefix->at(bucket);
            size_t bucketOffset = *bucketOffsetPtr;

            if (recurseBucket[bucket]) {
                size_t childPos = recurseBucketRank->rank1(bucket);
                return bucketOffset + children.at(childPos)->operator()(string);
            } else {
                // Perform direct rank storing
                ++bucketOffsetPtr;
                size_t nextBucketOffset = *bucketOffsetPtr;
                size_t bucketSize = nextBucketOffset - bucketOffset;
                return bucketOffset + retrieval.query(bucketSize, util::MurmurHash64(string));
            }
        }
};
