#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
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

        struct TreeNode {
            // As in normal DirectRankStoring:
            SuccinctPgmBucketMapper *bucketMapper = nullptr;
            util::EliasFanoM *bucketSizePrefix = nullptr;

            // If a bucket is very full, recurse. This is treated as pointers to the next layer data structure.
            std::vector<uint16_t> children; // TODO: Succinct tree representation

            // Cut off characters that are the same in all input strings. Useful especially for recursion.
            size_t minLCP = 9999999;
        };

        MultiRetrievalDataStructure retrieval;
        size_t N = 0;
        std::vector<TreeNode> treeNodes;

    public:
        explicit RecursiveDirectRankStoringMmphf(const std::vector<std::string> &strings) {
            N = strings.size();
            constructNode(strings.begin(), strings.end(), 0);
            retrieval.build();
        }
    private:
        size_t constructNode(const auto begin, const auto end, const size_t knownCommonPrefixLength) {
            size_t nodeId = treeNodes.size();
            treeNodes.emplace_back();
            TreeNode treeNode;

            assert(std::distance(begin, end) >= 2);
            auto it = begin + 1;
            treeNode.minLCP = (*begin).length();
            while (it != end) {
                size_t lengthToCheck = std::min(treeNode.minLCP, (*it).length());
                const char* s1ptr = (*it).data();
                const char* s2ptr = (*std::prev(it)).data();
                size_t lcp = knownCommonPrefixLength;
                while (lcp < lengthToCheck && s1ptr[lcp] == s2ptr[lcp]) {
                    lcp++;
                }
                treeNode.minLCP = std::min(treeNode.minLCP, lcp);
                it++;
            }

            // Trim strings and extract chunks
            std::vector<uint64_t> chunks;
            uint64_t previousChunk = 0;
            it = begin;
            while (it != end) {
                uint64_t chunk = readChunk((*it).c_str() + treeNode.minLCP, (*it).length() - treeNode.minLCP, 8);
                assert(chunk >= previousChunk);
                if (chunk != previousChunk || chunks.empty()) {
                    chunks.push_back(chunk);
                    previousChunk = chunk;
                }
                it++;
            }
            assert(chunks.size() >= 2); // If all were the same, we would have not cut off enough
            treeNode.bucketMapper = new SuccinctPgmBucketMapper(chunks.begin(), chunks.end());
            assert(treeNode.bucketMapper->numBuckets >= 2);

            it = begin;
            auto currentBucketBegin = begin;
            size_t prev_bucket = 0;
            size_t bucketSizePrefixTemp = 0;
            treeNode.bucketSizePrefix = new util::EliasFanoM(treeNode.bucketMapper->numBuckets + 1, std::distance(begin, end) + 1);

            auto constructBucket = [&] {
                size_t currentBucketSize = std::distance(currentBucketBegin, it);
                if (currentBucketSize <= 1) {
                    // Nothing to do
                    treeNode.children.push_back(-1);
                } else if (currentBucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                    // Perform direct rank storing
                    uint32_t indexInBucket = 0;
                    while (currentBucketBegin != it) {
                        retrieval.addInput(currentBucketSize,
                                    util::remix(util::MurmurHash64(*currentBucketBegin) + nodeId), indexInBucket);
                        currentBucketBegin++;
                        indexInBucket++;
                    }
                    treeNode.children.push_back(-1);
                } else {
                    // Recurse
                    size_t nodeId = constructNode(currentBucketBegin, it, treeNode.minLCP);
                    treeNode.children.push_back(nodeId);
                }

                treeNode.bucketSizePrefix->push_back(bucketSizePrefixTemp);
                currentBucketBegin = it;
                bucketSizePrefixTemp += currentBucketSize;
            };
            treeNode.bucketMapper->bucketOf(chunks.begin(), chunks.end(), [&](auto chunks_it, size_t bucket) {
                while (it != end && readChunk(it->c_str() + treeNode.minLCP, it->length() - treeNode.minLCP, 8) < *chunks_it) {
                    ++it;
                }
                while (bucket != prev_bucket) {
                    constructBucket();
                    prev_bucket++;
                }
            });
            it = end;
            while (prev_bucket < treeNode.bucketMapper->numBuckets + 1) {
                constructBucket();
                prev_bucket++;
            }
            treeNode.bucketSizePrefix->buildRankSelect();

            treeNodes.at(nodeId) = treeNode;
            return nodeId;
        }

    public:
        ~RecursiveDirectRankStoringMmphf() {
            for (TreeNode &node : treeNodes) {
                delete node.bucketMapper;
                delete node.bucketSizePrefix;
            }
        }

        static std::string name() {
            return "RecursiveDirectRankStoringMmphf";
        }

        size_t spaceBits() {
            std::cout<<"Retrieval:             "<<1.0*retrieval.spaceBits()/N<<std::endl;
            std::cout<<"Bucket size prefix:    "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, TreeNode &node) { return size + node.bucketSizePrefix->space(); }) / N << std::endl;
            std::cout<<"PGM:                   "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, TreeNode &node) { return size + node.bucketMapper->size(); }) / N<<std::endl;
            std::cout<<"Child pointers:        "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, TreeNode &node) { return size + node.children.size()*sizeof(uint16_t); }) / N<<std::endl;
            std::cout<<"sizeof(*this + nodes): "<<8.0*(sizeof(*this) + treeNodes.size() * sizeof(TreeNode))/N<<std::endl;
            std::cout<<"Total tree nodes:      "<<treeNodes.size()<<std::endl;

            return retrieval.spaceBits()
                        + 8 * (sizeof(*this) + treeNodes.size() * sizeof(TreeNode))
                        + 8 * std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                          [] (size_t size, TreeNode &node) {
                            return size + node.children.size() * sizeof(uint16_t)
                                        + node.bucketMapper->size()
                                        + node.bucketSizePrefix->space();
                        });
        }

        uint64_t operator ()(const std::string &string) {
            uint64_t stringHash = util::MurmurHash64(string);

            TreeNode node = treeNodes.at(0);
            size_t nodeId = 0;
            size_t keysBefore = 0;
            while (true) {
                uint64_t chunk = readChunk(string.c_str() + node.minLCP, string.length() - node.minLCP, 8);
                size_t bucket = node.bucketMapper->bucketOf(chunk);
                size_t bucketOffset = *node.bucketSizePrefix->at(bucket);
                size_t nextBucketOffset = *node.bucketSizePrefix->at(bucket + 1);
                size_t bucketSize = nextBucketOffset - bucketOffset;
                keysBefore += bucketOffset;

                if (bucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                    // Perform direct rank storing
                    if (bucketSize == 1) {
                        return keysBefore;
                    } else {
                        return keysBefore + retrieval.query(bucketSize, util::remix(stringHash + nodeId));
                    }
                } else {
                    nodeId = node.children.at(bucket);
                    node = treeNodes.at(nodeId);
                }
            }
        }
};
