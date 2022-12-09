#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include "bucket_mapping/SuccinctPgmBucketMapper.hpp"
#include "DirectRankStoring.hpp"
#include "bucket_mapping/support/EliasFanoModified.hpp"
#include <MurmurHash64.h>
#include <support/BitVectorBuilder.hpp>
#include <support/SuccinctTree.hpp>

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
            SuccinctPgmBucketMapper *bucketMapper = nullptr;
            size_t minLCP = 9999999;
        };

        MultiRetrievalDataStructure retrieval;
        size_t N = 0;
        std::vector<TreeNode> treeNodes;
        util::EliasFanoM *bucketOffsets;
        SuccinctTree tree;
    public:
        explicit RecursiveDirectRankStoringMmphf(const std::vector<std::string> &strings) {
            std::vector<std::size_t> bucketOffsetsInput;
            N = strings.size();
            constructNode(strings.begin(), strings.end(), 0, 0, bucketOffsetsInput);
            retrieval.build();
            tree.build();

            bucketOffsetsInput.push_back(N);
            bucketOffsets = new util::EliasFanoM(bucketOffsetsInput.size(), N + 1);
            for (size_t i : bucketOffsetsInput) {
                bucketOffsets->push_back(i);
            }
            bucketOffsets->buildRankSelect();
        }
    private:
        void constructNode(const auto begin, const auto end, const size_t knownCommonPrefixLength, size_t offset,
                           std::vector<std::size_t> &bucketOffsetsInput) {
            size_t nodeId = treeNodes.size();
            treeNodes.emplace_back();
            TreeNode treeNode;
            tree.openInnerNode();

            bucketOffsetsInput.push_back(offset);

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

            auto constructBucket = [&] {
                size_t currentBucketSize = std::distance(currentBucketBegin, it);
                if (currentBucketSize <= 1) {
                    // Nothing to do
                    tree.appendLeaf();
                    bucketOffsetsInput.push_back(offset + bucketSizePrefixTemp);
                } else if (currentBucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                    // Perform direct rank storing
                    uint32_t indexInBucket = 0;
                    while (currentBucketBegin != it) {
                        retrieval.addInput(currentBucketSize, util::MurmurHash64(*currentBucketBegin), indexInBucket);
                        currentBucketBegin++;
                        indexInBucket++;
                    }
                    tree.appendLeaf();
                    bucketOffsetsInput.push_back(offset + bucketSizePrefixTemp);
                } else {
                    // Recurse
                    constructNode(currentBucketBegin, it, treeNode.minLCP, offset + bucketSizePrefixTemp, bucketOffsetsInput);
                }

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
            while (prev_bucket < treeNode.bucketMapper->numBuckets) {
                constructBucket();
                prev_bucket++;
            }

            treeNodes.at(nodeId) = treeNode;
            tree.closeInnerNode();
        }

    public:
        ~RecursiveDirectRankStoringMmphf() {
            for (TreeNode &node : treeNodes) {
                delete node.bucketMapper;
            }
        }

        static std::string name() {
            return "RecursiveDirectRankStoringMmphf";
        }

        size_t spaceBits() {
            std::cout<<"Retrieval:             "<<1.0*retrieval.spaceBits()/N<<std::endl;
            std::cout<<"PGM in tree nodes:     "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, TreeNode &node) { return size + node.bucketMapper->size(); }) / N<<std::endl;
            std::cout<<"Bucket size prefix:    "<<8.0*bucketOffsets->space()/N << std::endl;
            std::cout<<"Tree representation:   "<<1.0*tree.spaceBits()/N<<std::endl;
            std::cout<<"Tree node data:        "<<8.0*(treeNodes.size() * sizeof(TreeNode))/N<<std::endl;

            return retrieval.spaceBits()
                        + 8 * (sizeof(*this) + treeNodes.size() * sizeof(TreeNode))
                        + tree.spaceBits()
                        + 8 * bucketOffsets->space()
                        + 8 * std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                  [] (size_t size, TreeNode &node) { return size + node.bucketMapper->size(); });
        }

        uint64_t operator ()(const std::string &string) {
            auto treeReader = SuccinctTree::Reader(tree);
            TreeNode node = treeNodes.at(treeReader.nodeId);
            while (true) {
                uint64_t chunk = readChunk(string.c_str() + node.minLCP, string.length() - node.minLCP, 8);
                size_t bucket = node.bucketMapper->bucketOf(chunk);

                treeReader.skipToNthChild(bucket);
                SuccinctTree::Reader nextBucket = treeReader;
                nextBucket.nextSibling();
                size_t bucketOffset = *bucketOffsets->at(treeReader.nodeId);
                size_t nextBucketOffset = *bucketOffsets->at(nextBucket.nodeId);
                size_t bucketSize = nextBucketOffset - bucketOffset;

                assert(bucketSize >= 1 && "Key not in original data set, bucket size is 0");
                if (bucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                    // Perform direct rank storing
                    if (bucketSize == 1) {
                        return bucketOffset;
                    } else {
                        return bucketOffset + retrieval.query(bucketSize, util::MurmurHash64(string));
                    }
                } else {
                    node = treeNodes.at(treeReader.innerNodeRank());
                }
            }
        }
};
