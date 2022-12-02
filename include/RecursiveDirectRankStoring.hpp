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
        // We use the same retrieval data structure for all children,
        // so we need to prefix strings with something unique in each layer.
        static uint32_t retrievalPrefixCounter;

        // As in normal DirectRankStoring:
        SuccinctPgmBucketMapper *bucketMapper = nullptr;
        util::EliasFanoM *bucketSizePrefix = nullptr;
        MultiRetrievalDataStructure *retrieval;
        uint32_t retrievalPrefix;
        bool isRibbonOwner = false;

        // If a bucket is very full, recurse
        util::EliasFanoM *recurseBucket = nullptr;
        std::vector<RecursiveDirectRankStoringMmphf *> children;

        // Cut off characters that are the same in all input strings. Useful especially for recursion.
        size_t minLCP = 9999999;
    public:
        explicit RecursiveDirectRankStoringMmphf(std::vector<std::string> strings,
                                                 MultiRetrievalDataStructure *retrieval_ = nullptr) {
            retrieval = retrieval_;
            if (retrieval == nullptr) {
                retrieval = new MultiRetrievalDataStructure;
                isRibbonOwner = true;
            }
            retrievalPrefix = retrievalPrefixCounter++;
            assert(retrievalPrefixCounter < 99999); // All counters need to have the same number of digits

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
            std::vector<size_t> recurseBucketInput;

            size_t bucketSizePrefixTemp = 0;
            for (size_t i = 0; i < buckets.size(); i++) {
                bucketSizePrefix->push_back(bucketSizePrefixTemp);
                bucketSizePrefixTemp += buckets.at(i).size();

                if (buckets.at(i).size() < DIRECT_RANK_STORING_THRESHOLD) {
                    // Perform direct rank storing
                    if (buckets.at(i).size() > 1) {
                        for (size_t k = 0; k < buckets.at(i).size(); k++) {
                            std::string prefixedKey = std::to_string(retrievalPrefix) + buckets.at(i).at(k);
                            retrieval->addInput(buckets.at(i).size(), util::MurmurHash64(prefixedKey), k);
                        }
                    }
                } else {
                    recurseBucketInput.push_back(i);
                    children.push_back(new RecursiveDirectRankStoringMmphf(buckets.at(i), retrieval));
                }
            }
            bucketSizePrefix->push_back(bucketSizePrefixTemp);
            bucketSizePrefix->buildRankSelect();
            if (!recurseBucketInput.empty()) {
                recurseBucket = new util::EliasFanoM(buckets.size(), buckets.size() + 1);
                for (size_t b : recurseBucketInput) {
                    recurseBucket->push_back(b);
                }
                recurseBucket->buildRankSelect();
            }
            if (isRibbonOwner) {
                retrieval->build();
            }
        }

        ~RecursiveDirectRankStoringMmphf() {
            delete bucketMapper;
            delete bucketSizePrefix;
            delete recurseBucket;
            for (auto *child : children) {
                delete child;
            }
        }

        static std::string name() {
            return "RecursiveDirectRankStoringMmphf";
        }

        size_t spaceBits() {
            size_t bits = 8 * bucketMapper->size()
                    + (isRibbonOwner ? retrieval->spaceBits() : 0)
                    + 8 * bucketSizePrefix->space()
                    + 8 * sizeof(*this)
                    + 8 * (recurseBucket == nullptr ? 0 : recurseBucket->space())
                    + 8 * children.size() * sizeof(uint64_t);
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
            ++bucketOffsetPtr;
            size_t nextBucketOffset = *bucketOffsetPtr;
            size_t bucketSize = nextBucketOffset - bucketOffset;

            if (bucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                // Perform direct rank storing
                if (bucketSize == 1) {
                    return bucketOffset;
                } else {
                    std::string prefixedKey = std::to_string(retrievalPrefix) + string;
                    return bucketOffset + retrieval->query(bucketSize, util::MurmurHash64(prefixedKey));
                }
            } else {
                size_t childPos = recurseBucket->predecessorPosition(bucket).index(); // Rank query
                return bucketOffset + children.at(childPos)->operator()(string);
            }
        }
};
uint32_t RecursiveDirectRankStoringMmphf::retrievalPrefixCounter = 10000;
