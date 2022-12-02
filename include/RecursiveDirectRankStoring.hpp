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
        MultiRetrievalDataStructure * const retrieval;
        uint32_t retrievalPrefix;
        const bool isTopLevel;

        // If a bucket is very full, recurse
        util::EliasFanoM *recursingBuckets = nullptr;
        std::vector<RecursiveDirectRankStoringMmphf *> children;

        // Cut off characters that are the same in all input strings. Useful especially for recursion.
        size_t minLCP = 9999999;
    public:
        explicit RecursiveDirectRankStoringMmphf(std::vector<std::string> &strings)
            : RecursiveDirectRankStoringMmphf(strings, new MultiRetrievalDataStructure, true) {
        }
    private:
        RecursiveDirectRankStoringMmphf(std::vector<std::string> strings,
                                        MultiRetrievalDataStructure *retrieval_, bool isTopLevel_ = false)
                 : retrieval(retrieval_), isTopLevel(isTopLevel_) {
            retrievalPrefix = retrievalPrefixCounter++;
            assert(retrievalPrefixCounter < 99999); // All counters need to have the same number of digits

            for (size_t i = 0; i < strings.size(); i++) {
                if (i > 0) {
                    minLCP = std::min(minLCP, LCP(strings.at(i), strings.at(i-1)));
                }
            }

            // Trim strings and extract chunks
            std::vector<uint64_t> chunks;
            uint64_t previousChunk = 0;
            for (auto & string : strings) {
                string = string.substr(minLCP);
                uint64_t chunk = readChunk(string.c_str(), string.length(), 8);
                assert(chunk >= previousChunk);
                if (chunk != previousChunk || chunks.empty()) {
                    chunks.push_back(chunk);
                    previousChunk = chunk;
                }
            }
            assert(chunks.size() >= 2); // If all were the same, we would have not cut off enough
            bucketMapper = new SuccinctPgmBucketMapper(chunks.begin(), chunks.end());

            std::vector<std::string> currentBucket;
            size_t prev_bucket = 0;
            size_t bucketSizePrefixTemp = 0;
            bucketSizePrefix = new util::EliasFanoM(bucketMapper->numBuckets + 1, strings.size() + 1);
            std::vector<size_t> recursingBucketsInput;

            auto constructBucket = [&] {
                if (currentBucket.size() <= 1) {
                    // Nothing to do
                } else if (currentBucket.size() < DIRECT_RANK_STORING_THRESHOLD) {
                    // Perform direct rank storing
                    for (size_t k = 0; k < currentBucket.size(); k++) {
                        std::string prefixedKey = std::to_string(retrievalPrefix) + currentBucket.at(k);
                        retrieval->addInput(currentBucket.size(), util::MurmurHash64(prefixedKey), k);
                    }
                } else {
                    // Recurse
                    recursingBucketsInput.push_back(prev_bucket);
                    children.push_back(new RecursiveDirectRankStoringMmphf(currentBucket, retrieval));
                }

                bucketSizePrefix->push_back(bucketSizePrefixTemp);
                bucketSizePrefixTemp += currentBucket.size();
                currentBucket.clear();
            };
            // This loop should be replaced with a variant of SuccinctPGM::for_each
            for (std::string &string : strings) {
                uint64_t chunk = readChunk(string.c_str(), string.length(), 8);
                size_t bucket = bucketMapper->bucketOf(chunk);
                while (bucket != prev_bucket) {
                    constructBucket();
                    prev_bucket++;
                }
                currentBucket.push_back(string);
            }
            while (prev_bucket < bucketMapper->numBuckets + 1) {
                constructBucket();
                prev_bucket++;
            }
            bucketSizePrefix->buildRankSelect();
            if (!recursingBucketsInput.empty()) {
                recursingBuckets = new util::EliasFanoM(bucketMapper->numBuckets, bucketMapper->numBuckets + 1);
                for (size_t b : recursingBucketsInput) {
                    recursingBuckets->push_back(b);
                }
                recursingBuckets->buildRankSelect();
            }
            if (isTopLevel) {
                retrieval->build();
            }
        }

    public:
        ~RecursiveDirectRankStoringMmphf() {
            delete bucketMapper;
            delete bucketSizePrefix;
            delete recursingBuckets;
            for (auto *child : children) {
                delete child;
            }
        }

        static std::string name() {
            return "RecursiveDirectRankStoringMmphf";
        }

        size_t spaceBits(size_t N = 0) {
            assert(isTopLevel);
            if (N != 0) {
                std::cout<<"Retrieval:          "<<1.0*retrieval->spaceBits()/N<<std::endl;
                std::cout<<"Bucket size prefix: "<<8.0*visit([] (auto &obj) { return obj.bucketSizePrefix->size(); })/N<<std::endl;
                std::cout<<"PGM:                "<<8.0*visit([] (auto &obj) { return obj.bucketMapper->size(); })/N<<std::endl;
                std::cout<<"Recursion pointers: "<<8.0*visit([] (auto &obj) { return (obj.recursingBuckets == nullptr ? 0 : obj.recursingBuckets->space())
                                                                            + obj.children.size() * sizeof(uint64_t); })/N<<std::endl;
                std::cout<<"sizeof(*this):      "<<8.0*visit([] (auto &obj) { return sizeof(obj); })/N<<std::endl;
            }

            return retrieval->spaceBits() + 8 * visit([] (RecursiveDirectRankStoringMmphf &obj) {
                return obj.bucketMapper->size()
                     + obj.bucketSizePrefix->space()
                     + sizeof(obj)
                     + (obj.recursingBuckets == nullptr ? 0 : obj.recursingBuckets->space())
                     + obj.children.size() * sizeof(uint64_t);
                });
        }

        template <typename T>
        size_t visit(T visitor) {
            size_t sum = visitor(*this);
            for (auto *child : children) {
                sum += child->visit(visitor);
            }
            return sum;
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
                size_t childPos = recursingBuckets->predecessorPosition(bucket).index(); // Rank query
                return bucketOffset + children.at(childPos)->operator()(string);
            }
        }
};
uint32_t RecursiveDirectRankStoringMmphf::retrievalPrefixCounter = 10000;
