#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include "bucket_mapping/SuccinctPgmBucketMapper.hpp"
#include "DirectRankStoring.hpp"
#include "bucket_mapping/support/EliasFanoModified.hpp"
#include <MurmurHash64.h>
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
        static constexpr size_t DIRECT_RANK_STORING_THRESHOLD = 1024;
        static constexpr size_t LINEAR_MAPPING_THRESHOLD = 500;

        struct TreeNode {
            BucketMapper *bucketMapper = nullptr;
            size_t minLCP = 9999999;
        };

        MultiRetrievalDataStructure retrieval;
        size_t N = 0;
        std::vector<TreeNode> treeNodes;
        util::EliasFanoM *bucketOffsets;
        SuccinctTree tree;
        std::vector<std::size_t> bucketOffsetsInput;
    public:
        explicit RecursiveDirectRankStoringMmphf(const std::vector<std::string> &strings) {
            N = strings.size();
            constructNode(strings.begin(), strings.end(), 0, 0);
            retrieval.build();
            tree.build();

            bucketOffsetsInput.push_back(N);
            bucketOffsets = new util::EliasFanoM(bucketOffsetsInput.size(), N + 1);
            for (size_t i : bucketOffsetsInput) {
                bucketOffsets->push_back(i);
            }
            bucketOffsets->buildRankSelect();
            bucketOffsetsInput.clear();
            bucketOffsetsInput.shrink_to_fit();
            treeNodes.shrink_to_fit();
        }
    private:
        void constructNode(const auto begin, const auto end, const size_t knownCommonPrefixLength, size_t offset) {
            size_t nodeId = treeNodes.size();
            treeNodes.emplace_back();
            TreeNode treeNode;
            tree.openInnerNode();

            bucketOffsetsInput.push_back(offset);
            treeNode.minLCP = findMinLCP(begin, end, knownCommonPrefixLength);

            std::vector<uint64_t> chunks;
            extractChunks(begin, end, chunks, treeNode.minLCP);
            assert(chunks.size() >= 2); // If all were the same, we would have not cut off enough
            if (chunks.size() <= LINEAR_MAPPING_THRESHOLD) {
                treeNode.bucketMapper = new LinearBucketMapper(chunks.begin(), chunks.end());
            } else {
                treeNode.bucketMapper = new SuccinctPgmBucketMapper(chunks.begin(), chunks.end());
            }
            assert(treeNode.bucketMapper->numBuckets() >= 2);

            auto it = begin;
            auto currentBucketBegin = begin;
            size_t prev_bucket = 0;
            size_t bucketSizePrefixTemp = 0;

            treeNode.bucketMapper->bucketOf(chunks.begin(), chunks.end(), [&](auto chunks_it, size_t bucket) {
                while (it != end && readChunk(it->c_str() + treeNode.minLCP, it->length() - treeNode.minLCP, 8) < *chunks_it) {
                    ++it;
                }
                while (prev_bucket < bucket) {
                    constructChild(currentBucketBegin, it, offset + bucketSizePrefixTemp, treeNode.minLCP);
                    bucketSizePrefixTemp += std::distance(currentBucketBegin, it);
                    currentBucketBegin = it;
                    prev_bucket++;
                }
            });
            it = end;
            while (prev_bucket < treeNode.bucketMapper->numBuckets()) {
                constructChild(currentBucketBegin, it, offset + bucketSizePrefixTemp, treeNode.minLCP);
                bucketSizePrefixTemp += std::distance(currentBucketBegin, it);
                currentBucketBegin = it;
                prev_bucket++;
            }

            treeNodes.at(nodeId) = treeNode;
            tree.closeInnerNode();
        }

        size_t findMinLCP(const auto begin, const auto end, size_t knownCommonPrefixLength) {
            assert(std::distance(begin, end) >= 2);
            auto it = begin + 1;
            size_t minLCP = (*begin).length();
            while (it != end) {
                size_t lengthToCheck = std::min(minLCP, (*it).length());
                const char* s1ptr = (*it).data();
                const char* s2ptr = (*std::prev(it)).data();
                size_t lcp = knownCommonPrefixLength;
                while (lcp < lengthToCheck && s1ptr[lcp] == s2ptr[lcp]) {
                    lcp++;
                }
                minLCP = std::min(minLCP, lcp);
                it++;
            }
            return minLCP;
        }

        void extractChunks(const auto begin, const auto end, std::vector<uint64_t> &chunks, size_t minLCP) {
            uint64_t previousChunk = 0;
            auto it = begin;
            while (it != end) {
                uint64_t chunk = readChunk((*it).c_str() + minLCP, (*it).length() - minLCP, 8);
                assert(chunk >= previousChunk);
                if (chunk != previousChunk || chunks.empty()) {
                    chunks.push_back(chunk);
                    previousChunk = chunk;
                }
                it++;
            }
        }

        void constructChild(auto begin, auto end, size_t offset, size_t minLCP) {
            size_t currentBucketSize = std::distance(begin, end);
            if (currentBucketSize <= 1) {
                // Nothing to do
                tree.appendLeaf();
                bucketOffsetsInput.push_back(offset);
            } else if (currentBucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                // Perform direct rank storing
                uint32_t indexInBucket = 0;
                auto it = begin;
                while (it != end) {
                    retrieval.addInput(currentBucketSize, util::MurmurHash64(*it), indexInBucket);
                    it++;
                    indexInBucket++;
                }
                tree.appendLeaf();
                bucketOffsetsInput.push_back(offset);
            } else {
                // Recurse
                constructNode(begin, end, minLCP, offset);
            }
        }

    public:
        ~RecursiveDirectRankStoringMmphf() {
            delete bucketOffsets;
            for (TreeNode &node : treeNodes) {
                delete node.bucketMapper;
            }
        }

        static std::string name() {
            return "RecursiveDirectRankStoringMmphf";
        }

        size_t spaceBits() {
            std::cout<<"Retrieval:           "<<1.0*retrieval.spaceBits()/N<<std::endl;
            std::cout<<"Bucket mappers:      "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, TreeNode &node) { return size + node.bucketMapper->size(); }) / N<<std::endl;
            std::cout<<"  PGM:                 "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, TreeNode &node) { return size + ((node.bucketMapper->numBuckets() > LINEAR_MAPPING_THRESHOLD) ? node.bucketMapper->size() : 0); }) / N<<std::endl;
            std::cout<<"  Linear:              "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, TreeNode &node) { return size + ((node.bucketMapper->numBuckets() <= LINEAR_MAPPING_THRESHOLD) ? node.bucketMapper->size() : 0); }) / N<<std::endl;
            std::cout<<"Bucket size prefix:  "<<8.0*bucketOffsets->space()/N << std::endl;
            std::cout<<"Tree representation: "<<1.0*tree.spaceBits()/N<<std::endl;
            std::cout<<"  Without indices:     "<<1.0*tree.spaceBitsWithoutIndices()/N<<std::endl;
            std::cout<<"Tree node data:      "<<8.0*(treeNodes.size() * sizeof(TreeNode))/N<<std::endl;

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

class RecursiveDirectRankStoringV2Mmphf {
    private:
        static constexpr size_t DIRECT_RANK_STORING_THRESHOLD = 1024;
        static constexpr size_t LINEAR_MAPPING_THRESHOLD = 500;

        struct TreeNode {
            using index_type = uint16_t;
            BucketMapper *bucketMapper = nullptr;
            index_type indexes[8] = {0, 0, 0, 0, 0, 0, 0, 0};

            void setIndexes(const auto &other) {
                if (other[0] > std::numeric_limits<index_type>::max())
                    throw std::runtime_error("Not yet implemented (1)");

                indexes[0] = other[0];
                for (size_t i = 1; i < other.size(); ++i) {
                    if (other[i] - other[i - 1] > std::numeric_limits<index_type>::max())
                        throw std::runtime_error("Not yet implemented (2)");
                    indexes[i] = other[i] - other[i - 1];
                }
            }

            auto getIndexes() const {
                std::vector<uint32_t> result;
                result.push_back(indexes[0]);
                for (size_t i = 1; i < 8 && indexes[i]; ++i)
                    result.push_back(result[i - 1] + indexes[i]);
                return result;
            }
        };

        MultiRetrievalDataStructure retrieval;
        size_t N = 0;
        std::vector<TreeNode> treeNodes;
        util::EliasFanoM *bucketOffsets;
        SuccinctTree tree;
        std::vector<std::size_t> bucketOffsetsInput;
    public:
        explicit RecursiveDirectRankStoringV2Mmphf(const std::vector<std::string> &strings) {
            N = strings.size();
            auto lcps = computeLCPs(strings.begin(), strings.end());
            constructNode(strings.begin(), strings.end(), lcps.begin(), lcps.end(), 0, 0);
            retrieval.build();
            tree.build();

            bucketOffsetsInput.push_back(N);
            bucketOffsets = new util::EliasFanoM(bucketOffsetsInput.size(), N + 1);
            for (size_t i : bucketOffsetsInput) {
                bucketOffsets->push_back(i);
            }
            bucketOffsets->buildRankSelect();
            bucketOffsetsInput.clear();
        }
    private:
        void constructNode(const auto begin, const auto end, const auto lcpsBegin, const auto lcpsEnd,
                           const size_t knownCommonPrefixLength, size_t offset) {
            size_t nodeId = treeNodes.size();
            treeNodes.emplace_back();
            TreeNode treeNode;
            tree.openInnerNode();

            bucketOffsetsInput.push_back(offset);
            auto indexes = distinctMinima(lcpsBegin + 1, lcpsEnd, 8, knownCommonPrefixLength);
            treeNode.setIndexes(indexes);

            std::vector<uint64_t> chunks;
            extractChunks(begin, end, chunks, indexes);
            assert(chunks.size() >= 2); // If all were the same, we would have not cut off enough
            if (chunks.size() <= LINEAR_MAPPING_THRESHOLD) {
                treeNode.bucketMapper = new LinearBucketMapper(chunks.begin(), chunks.end());
            } else {
                treeNode.bucketMapper = new SuccinctPgmBucketMapper(chunks.begin(), chunks.end());
            }
            assert(treeNode.bucketMapper->numBuckets() >= 2);

            auto it = begin;
            auto currentBucketBegin = begin;
            size_t prev_bucket = 0;
            size_t bucketSizePrefixTemp = 0;

            treeNode.bucketMapper->bucketOf(chunks.begin(), chunks.end(), [&](auto chunks_it, size_t bucket) {
                while (it != end && readChunk(it->c_str(), it->length(), indexes.begin(), indexes.end()) < *chunks_it) {
                    ++it;
                }
                while (prev_bucket < bucket) {
                    constructChild(currentBucketBegin, it,
                                   lcpsBegin + (currentBucketBegin - begin), lcpsBegin + (it - begin),
                                   offset + bucketSizePrefixTemp, indexes[0]);
                    bucketSizePrefixTemp += std::distance(currentBucketBegin, it);
                    currentBucketBegin = it;
                    prev_bucket++;
                }
            });
            it = end;
            while (prev_bucket < treeNode.bucketMapper->numBuckets()) {
                constructChild(currentBucketBegin, it,
                               lcpsBegin + (currentBucketBegin - begin), lcpsBegin + (it - begin),
                               offset + bucketSizePrefixTemp, indexes[0]);
                bucketSizePrefixTemp += std::distance(currentBucketBegin, it);
                currentBucketBegin = it;
                prev_bucket++;
            }

            treeNodes.at(nodeId) = treeNode;
            tree.closeInnerNode();
        }

        void extractChunks(const auto begin, const auto end, std::vector<uint64_t> &chunks, const auto &indexes) {
            uint64_t previousChunk = 0;
            auto it = begin;
            while (it != end) {
                uint64_t chunk = readChunk(it->c_str(), it->length(), indexes.begin(), indexes.end());
                assert(chunk >= previousChunk);
                if (chunk != previousChunk || chunks.empty()) {
                    chunks.push_back(chunk);
                    previousChunk = chunk;
                }
                it++;
            }
        }

        void constructChild(auto begin, auto end, auto lcpsBegin, auto lcpsEnd, size_t offset, size_t minLCP) {
            size_t currentBucketSize = std::distance(begin, end);
            if (currentBucketSize <= 1) {
                // Nothing to do
                tree.appendLeaf();
                bucketOffsetsInput.push_back(offset);
            } else if (currentBucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                // Perform direct rank storing
                uint32_t indexInBucket = 0;
                auto it = begin;
                while (it != end) {
                    retrieval.addInput(currentBucketSize, util::MurmurHash64(*it), indexInBucket);
                    it++;
                    indexInBucket++;
                }
                tree.appendLeaf();
                bucketOffsetsInput.push_back(offset);
            } else {
                // Recurse
                constructNode(begin, end, lcpsBegin, lcpsEnd, minLCP, offset);
            }
        }

    public:
        ~RecursiveDirectRankStoringV2Mmphf() {
            delete bucketOffsets;
            for (TreeNode &node : treeNodes) {
                delete node.bucketMapper;
            }
        }

        static std::string name() {
            return "RecursiveDirectRankStoringV2Mmphf";
        }

        size_t spaceBits() {
            std::cout<<"Retrieval:           "<<1.0*retrieval.spaceBits()/N<<std::endl;
            std::cout<<"Bucket mappers:      "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, TreeNode &node) { return size + node.bucketMapper->size(); }) / N<<std::endl;
            std::cout<<"  PGM:                 "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, TreeNode &node) { return size + ((node.bucketMapper->numBuckets() > LINEAR_MAPPING_THRESHOLD) ? node.bucketMapper->size() : 0); }) / N<<std::endl;
            std::cout<<"  Linear:              "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, TreeNode &node) { return size + ((node.bucketMapper->numBuckets() <= LINEAR_MAPPING_THRESHOLD) ? node.bucketMapper->size() : 0); }) / N<<std::endl;
            std::cout<<"Bucket size prefix:  "<<8.0*bucketOffsets->space()/N << std::endl;
            std::cout<<"Tree representation: "<<1.0*tree.spaceBits()/N<<std::endl;
            std::cout<<"  Without indices:     "<<1.0*tree.spaceBitsWithoutIndices()/N<<std::endl;
            std::cout<<"Tree node data:      "<<8.0*(treeNodes.size() * sizeof(TreeNode))/N<<std::endl;

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
                auto indexes = node.getIndexes();
                uint64_t chunk = readChunk(string.c_str(), string.length(), indexes.begin(), indexes.end());
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