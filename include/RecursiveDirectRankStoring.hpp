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
        MultiRetrievalDataStructure * const retrieval;
        uint32_t retrievalSeed; // Common retrieval data structure for all layers ==> each layer needs a unique seed
        const bool isTopLevel;

        // If a bucket is very full, recurse
        util::EliasFanoM *recursingBuckets = nullptr;
        std::vector<RecursiveDirectRankStoringMmphf *> children;

        // Cut off characters that are the same in all input strings. Useful especially for recursion.
        size_t minLCP = 9999999;
    public:
        explicit RecursiveDirectRankStoringMmphf(std::vector<std::string> &strings)
            : RecursiveDirectRankStoringMmphf(strings.begin(), strings.end(),
                                              new MultiRetrievalDataStructure, true, 0) {
        }
    private:
        RecursiveDirectRankStoringMmphf(const auto begin, const auto end, MultiRetrievalDataStructure *retrieval_,
                                        const bool isTopLevel_, const size_t knownCommonPrefixLength)
                 : retrieval(retrieval_), isTopLevel(isTopLevel_) {
            static uint32_t nextRetrievalSeed = 0;
            retrievalSeed = nextRetrievalSeed++;

            assert(std::distance(begin, end) >= 2);
            auto it = begin + 1;
            minLCP = (*begin).length();
            while (it != end) {
                size_t lengthToCheck = std::min(minLCP, (*it).length());
                char* s1ptr = (*it).data();
                char* s2ptr = (*std::prev(it)).data();
                size_t lcp = knownCommonPrefixLength;
                while (lcp < lengthToCheck && s1ptr[lcp] == s2ptr[lcp]) {
                    lcp++;
                }
                minLCP = std::min(minLCP, lcp);
                it++;
            }

            // Trim strings and extract chunks
            std::vector<uint64_t> chunks;
            uint64_t previousChunk = 0;
            it = begin;
            while (it != end) {
                uint64_t chunk = readChunk((*it).c_str() + minLCP, (*it).length() - minLCP, 8);
                assert(chunk >= previousChunk);
                if (chunk != previousChunk || chunks.empty()) {
                    chunks.push_back(chunk);
                    previousChunk = chunk;
                }
                it++;
            }
            assert(chunks.size() >= 2); // If all were the same, we would have not cut off enough
            bucketMapper = new SuccinctPgmBucketMapper(chunks.begin(), chunks.end());

            it = begin;
            auto currentBucketBegin = begin;
            size_t prev_bucket = 0;
            size_t bucketSizePrefixTemp = 0;
            bucketSizePrefix = new util::EliasFanoM(bucketMapper->numBuckets + 1, std::distance(begin, end) + 1);
            std::vector<size_t> recursingBucketsInput;

            auto constructBucket = [&] {
                size_t currentBucketSize = std::distance(currentBucketBegin, it);
                if (currentBucketSize <= 1) {
                    // Nothing to do
                } else if (currentBucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                    // Perform direct rank storing
                    uint32_t indexInBucket = 0;
                    while (currentBucketBegin != it) {
                        retrieval->addInput(currentBucketSize,
                                    util::remix(util::MurmurHash64(*currentBucketBegin) + retrievalSeed), indexInBucket);
                        currentBucketBegin++;
                        indexInBucket++;
                    }
                } else {
                    // Recurse
                    recursingBucketsInput.push_back(prev_bucket);
                    children.push_back(new RecursiveDirectRankStoringMmphf(
                                currentBucketBegin, it, retrieval, false, minLCP));
                }

                bucketSizePrefix->push_back(bucketSizePrefixTemp);
                currentBucketBegin = it;
                bucketSizePrefixTemp += currentBucketSize;
            };
            bucketMapper->bucketOf(chunks.begin(), chunks.end(), [&](auto chunks_it, size_t bucket) {
                while (it != end && readChunk(it->c_str() + minLCP, it->length() - minLCP, 8) < *chunks_it) {
                    ++it;
                }
                while (bucket != prev_bucket) {
                    constructBucket();
                    prev_bucket++;
                }
            });
            it = end;
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

        size_t spaceBits() {
            assert(isTopLevel);
            size_t N = *(bucketSizePrefix->at(bucketSizePrefix->size() - 1));
            std::cout<<"Retrieval:          "<<1.0*retrieval->spaceBits()/N<<std::endl;
            std::cout<<"Bucket size prefix: "<<8.0*visit([] (auto &obj) { return obj.bucketSizePrefix->size(); })/N<<std::endl;
            std::cout<<"PGM:                "<<8.0*visit([] (auto &obj) { return obj.bucketMapper->size(); })/N<<std::endl;
            std::cout<<"Recursion pointers: "<<8.0*visit([] (auto &obj) { return (obj.recursingBuckets == nullptr ? 0 : obj.recursingBuckets->space())
                                                                        + obj.children.size() * sizeof(uint64_t); })/N<<std::endl;
            std::cout<<"sizeof(*this):      "<<8.0*visit([] (auto &obj) { return sizeof(obj); })/N<<std::endl;

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

        uint64_t operator ()(std::string &string) {
            return operator()(string, util::MurmurHash64(string));
        }

    private:
        uint64_t operator ()(std::string &string, uint64_t stringHash) {
            uint64_t chunk = readChunk(string.c_str() + minLCP, string.length() - minLCP, 8);
            size_t bucket = bucketMapper->bucketOf(chunk);
            size_t bucketOffset = *bucketSizePrefix->at(bucket);
            size_t nextBucketOffset = *bucketSizePrefix->at(bucket + 1);
            size_t bucketSize = nextBucketOffset - bucketOffset;

            if (bucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                // Perform direct rank storing
                if (bucketSize == 1) {
                    return bucketOffset;
                } else {
                    return bucketOffset + retrieval->query(bucketSize, util::remix(stringHash + retrievalSeed));
                }
            } else {
                size_t childPos = recursingBuckets->predecessorPosition(bucket).index(); // Rank query
                return bucketOffset + children.at(childPos)->operator()(string, stringHash);
            }
        }
};
