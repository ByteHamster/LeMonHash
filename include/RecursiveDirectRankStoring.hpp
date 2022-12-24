#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include "bucket_mapping/SuccinctPgmBucketMapper.hpp"
#include "bucket_mapping/ShortPgmBucketMapper.hpp"
#include "DirectRankStoring.hpp"
#include "bucket_mapping/support/EliasFanoModified.hpp"
#include <MurmurHash64.h>
#include "support/PartitionedEliasFano.hpp"

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
        // For fewer objects, it is not worth creating a tree structure. Just store the ranks explicitly.
        static constexpr size_t DIRECT_RANK_STORING_THRESHOLD = 128;

        // For fewer than 25 chunks, it is cheaper storing the ranks explicitly than mapping them - even with the LinearBucketMapper.
        // sizeof(LinearBucketMapper) < 25×ceil(log2(25))
        // For some more it is not explicitly cheaper but given that the ranks are "perfect", it is still worth it.
        static constexpr size_t CHUNK_DIRECT_RANK_STORING_THRESHOLD = 128;

        const std::string PATH_ROOT = "root";

        struct TreeNode {
            union {
                ShortPgmBucketMapper *bucketMapper = nullptr;
                size_t numChildren;
            };
            size_t offsetsOffset : 48 = 0;
            size_t minLCP : 15 = ((1 << 15) - 1);
            bool directRankStoring : 1 = false;
        };

        MultiRetrievalDataStructure retrieval;
        size_t N = 0;
        size_t height = 0;
        std::unordered_map<std::string, TreeNode> treeNodes; // TODO: Use perfect hashing
        std::vector<PartitionedEliasFano *> bucketOffsets;
        std::vector<std::vector<size_t>> bucketOffsetsInput;
    public:
        explicit RecursiveDirectRankStoringMmphf(const std::vector<std::string> &strings) {
            N = strings.size();
            constructNode(strings.begin(), strings.end(), 0, 0, PATH_ROOT, 0);
            retrieval.build();

            for (std::vector<size_t> v : bucketOffsetsInput) {
                auto *fano = new PartitionedEliasFano(v);
                bucketOffsets.push_back(fano);
            }
            height = bucketOffsetsInput.size();
            bucketOffsetsInput.clear();
            bucketOffsetsInput.shrink_to_fit();

            //std::ofstream myfile("tree.dot");
            //exportTreeStructure(myfile);
            //myfile.close();
        }
    private:
        void constructNode(const auto begin, const auto end, const size_t knownCommonPrefixLength,
                           size_t offset, const std::string &path, size_t layer) {
            TreeNode treeNode;
            if (bucketOffsetsInput.size() <= layer) {
                bucketOffsetsInput.resize(layer + 1);
            }
            treeNode.offsetsOffset = bucketOffsetsInput.at(layer).size();
            treeNode.minLCP = findMinLCP(begin, end, knownCommonPrefixLength);
            size_t nThisNode = std::distance(begin, end);

            std::vector<uint64_t> chunks;
            extractChunks(begin, end, chunks, treeNode.minLCP);
            assert(chunks.size() >= 2); // If all were the same, we would have not cut off enough

            if (chunks.size() < CHUNK_DIRECT_RANK_STORING_THRESHOLD) {
                treeNode.directRankStoring = true;
                treeNode.numChildren = chunks.size();

                auto it = begin;
                size_t bucketSizePrefixTemp = 0;
                for (size_t chunk = 0; chunk < chunks.size(); chunk++) {
                    uint64_t chunkValue = chunks.at(chunk);
                    auto currentBucketBegin = it;
                    while (it != end &&
                           readChunk(it->c_str() + treeNode.minLCP, it->length() - treeNode.minLCP, 8) == chunkValue) {
                        ++it;
                    }
                    std::string key = "Chunk" + path + "/" + std::to_string(chunkValue);
                    retrieval.addInput(chunks.size(), key, chunk);
                    constructChild(currentBucketBegin, it, offset + bucketSizePrefixTemp, treeNode.minLCP,
                                   path + "/" + std::to_string(chunk), layer);
                    bucketSizePrefixTemp += std::distance(currentBucketBegin, it);
                }
            } else {
                treeNode.directRankStoring = false;
                treeNode.bucketMapper = new ShortPgmBucketMapper(chunks.begin(), chunks.end());
                size_t numBuckets = treeNode.bucketMapper->numBuckets();
                assert(numBuckets >= 2);

                auto it = begin;
                auto currentBucketBegin = begin;
                size_t prev_bucket = 0;
                size_t bucketSizePrefixTemp = 0;

                treeNode.bucketMapper->bucketOf(chunks.begin(), chunks.end(), [&](auto chunks_it, size_t bucket) {
                    while (it != end &&
                           readChunk(it->c_str() + treeNode.minLCP, it->length() - treeNode.minLCP, 8) < *chunks_it) {
                        ++it;
                    }
                    while (prev_bucket < bucket) {
                        constructChild(currentBucketBegin, it, offset + bucketSizePrefixTemp, treeNode.minLCP,
                                       path + "/" + std::to_string(prev_bucket), layer);
                        bucketSizePrefixTemp += std::distance(currentBucketBegin, it);
                        currentBucketBegin = it;
                        prev_bucket++;
                    }
                });
                it = end;
                while (prev_bucket < numBuckets) {
                    constructChild(currentBucketBegin, it, offset + bucketSizePrefixTemp, treeNode.minLCP,
                                   path + "/" + std::to_string(prev_bucket), layer);
                    bucketSizePrefixTemp += std::distance(currentBucketBegin, it);
                    currentBucketBegin = it;
                    prev_bucket++;
                }
            }

            bucketOffsetsInput.at(layer).push_back(offset + nThisNode);
            treeNodes.insert(std::make_pair(path, treeNode));
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

        void constructChild(auto begin, auto end, size_t offset, size_t minLCP, std::string path, size_t layer) {
            bucketOffsetsInput.at(layer).push_back(offset);
            size_t currentBucketSize = std::distance(begin, end);
            if (currentBucketSize <= 1) {
                // Nothing to do
            } else if (currentBucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                // Perform direct rank storing
                uint32_t indexInBucket = 0;
                auto it = begin;
                while (it != end) {
                    retrieval.addInput(currentBucketSize, util::MurmurHash64(*it), indexInBucket);
                    it++;
                    indexInBucket++;
                }
            } else {
                // Recurse
                constructNode(begin, end, minLCP, offset, path, layer + 1);
            }
        }

    public:
        ~RecursiveDirectRankStoringMmphf() {
            for (auto fano : bucketOffsets) {
                delete fano;
            }
            for (auto &node : treeNodes) {
                if (!node.second.directRankStoring) {
                    delete node.second.bucketMapper;
                }
            }
        }

        static std::string name() {
            return "RecursiveDirectRankStoringMmphf";
        }

        size_t spaceBits() {
            std::cout<<"Retrieval:           "<<1.0*retrieval.spaceBits()/N<<std::endl;
            std::cout<<"Bucket mapper:       "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                 [] (size_t size, auto &node) { return size + (node.second.directRankStoring ? 0 : node.second.bucketMapper->size()); }) / N<<std::endl;
            std::cout<<"Bucket offsets:      "<<1.0 * std::accumulate(bucketOffsets.begin(), bucketOffsets.end(), 0,
                                 [] (size_t size, PartitionedEliasFano *fano) { return size + fano->bit_size(); }) / N<<std::endl;
            std::cout<<"Tree node data:      "<<(8.0 * treeNodes.size() * sizeof(TreeNode) + 3.0 * treeNodes.size())/N<<std::endl;
            std::cout<<"Height:              "<<height<<std::endl;

            return retrieval.spaceBits()
                        + 8 * (sizeof(*this)
                            + treeNodes.size() * sizeof(TreeNode)
                            + bucketOffsets.size() * sizeof(void*))
                        + std::accumulate(bucketOffsets.begin(), bucketOffsets.end(), 0,
                                    [] (size_t size, PartitionedEliasFano *fano) { return size + fano->bit_size(); })
                        + 8 * std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                    [] (size_t size, auto &node) { return size + (node.second.directRankStoring ? 0 : node.second.bucketMapper->size()); })
                        + 3 * treeNodes.size(); // Use an MPHF instead of std::unordered_map
        }

        void exportTreeStructure(std::ostream &os) {
            os<<"digraph {"<<std::endl;
            exportTreeStructureInternal(os, PATH_ROOT, 0);
            os<<"}"<<std::endl;
        }

        void exportTreeStructureInternal(std::ostream &os, std::string path, size_t layer) {
            float scaleY = 200;
            float scaleX = 1.6 * (bucketOffsets.size() * scaleY) / N;

            TreeNode node = treeNodes.at(path);
            size_t beginX = bucketOffsets.at(layer)->at(node.offsetsOffset);
            size_t endX = bucketOffsets.at(layer)->at(node.offsetsOffset + node.bucketMapper->numBuckets());
            os<<"  \""<<path<<"\" [ "<<std::endl;
            os<<"    pos = \""<<+scaleX*(beginX+endX)/2<<","<<scaleY*layer<<"\""<<std::endl;
            os<<"    layer = \""<<+layer<<"\""<<std::endl;
            os<<"    label = \""<<+(endX - beginX)<<"\""<<std::endl;
            os<<"  ]"<<std::endl;

            for (size_t i = 0; i < node.bucketMapper->numBuckets(); i++) {
                std::string childPath = path + "/" + std::to_string(i);
                os<<"  \""<<path<<"\" -> \""<<childPath<<"\""<<std::endl;
                size_t childBegin = bucketOffsets.at(layer)->at(node.offsetsOffset + i);
                size_t childSize = bucketOffsets.at(layer)->at(node.offsetsOffset + i + 1) - childBegin;
                if (childSize < DIRECT_RANK_STORING_THRESHOLD) {
                    os<<"  \""<<childPath<<"\" [ "<<std::endl;
                    os<<"    pos = \""<<+scaleX*childBegin<<","<<scaleY*(layer+1)<<"\""<<std::endl;
                    os<<"    layer = \""<<+(layer+1)<<"\""<<std::endl;
                    os<<"    label = \""<<+childSize<<"\""<<std::endl;
                    os<<"  ]"<<std::endl;
                } else {
                    exportTreeStructureInternal(os, childPath, layer + 1);
                }
            }
        }

        uint64_t operator ()(const std::string &string) {
            std::string path = PATH_ROOT;
            TreeNode node = treeNodes.at(path);
            size_t layer = 0;
            while (true) {
                uint64_t chunk = readChunk(string.c_str() + node.minLCP, string.length() - node.minLCP, 8);
                size_t bucket;
                if (node.directRankStoring) {
                    std::string key = "Chunk" + path + "/" + std::to_string(chunk);
                    bucket = retrieval.query(node.numChildren, key);
                } else {
                    bucket = node.bucketMapper->bucketOf(chunk);
                }
                size_t bucketOffset = bucketOffsets.at(layer)->at(node.offsetsOffset + bucket);
                size_t nextBucketOffset = bucketOffsets.at(layer)->at(node.offsetsOffset + bucket + 1);
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
                    layer++;
                    path += "/" + std::to_string(bucket);
                    node = treeNodes.at(path);
                }
            }
        }
};

class RecursiveDirectRankStoringV2Mmphf {
    private:
        static constexpr size_t DIRECT_RANK_STORING_THRESHOLD = 128;

        static constexpr size_t CHUNK_DIRECT_RANK_STORING_THRESHOLD = 128;

        const std::string PATH_ROOT = "root";

        #pragma pack(push, 1)
        struct TreeNode {
            using index_type = uint16_t;
            union {
                ShortPgmBucketMapper *bucketMapper = nullptr;
                size_t numChildren;
            };
            bool directRankStoring : 1 = false;
            uint32_t offsetsOffset : 31 = 0;
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
        #pragma pack(pop)

        MultiRetrievalDataStructure retrieval;
        size_t N = 0;
        size_t height = 0;
        std::unordered_map<std::string, TreeNode> treeNodes; // TODO: Use perfect hashing
        std::vector<PartitionedEliasFano *> bucketOffsets;
        std::vector<std::vector<size_t>> bucketOffsetsInput;
    public:
        explicit RecursiveDirectRankStoringV2Mmphf(const std::vector<std::string> &strings) {
            N = strings.size();
            auto lcps = computeLCPs(strings.begin(), strings.end());
            constructNode(strings.begin(), strings.end(), lcps.begin(), lcps.end(), 0, 0, PATH_ROOT, 0);
            retrieval.build();

            for (std::vector<size_t> v : bucketOffsetsInput) {
                auto *fano = new PartitionedEliasFano(v);
                bucketOffsets.push_back(fano);
            }
            height = bucketOffsets.size();
            bucketOffsetsInput.clear();
            bucketOffsetsInput.shrink_to_fit();

            //std::ofstream myfile("tree2.dot");
            //exportTreeStructure(myfile);
            //myfile.close();
        }
    private:
        void constructNode(const auto begin, const auto end, const auto lcpsBegin, const auto lcpsEnd,
                           const size_t knownCommonPrefixLength, size_t offset, const std::string &path, size_t layer) {
            TreeNode treeNode;
            auto indexes = distinctMinima(lcpsBegin + 1, lcpsEnd, 8, knownCommonPrefixLength);
            treeNode.setIndexes(indexes);
            if (bucketOffsetsInput.size() <= layer) {
                bucketOffsetsInput.resize(layer + 1);
            }
            treeNode.offsetsOffset = bucketOffsetsInput.at(layer).size();
            size_t nThisNode = std::distance(begin, end);

            std::vector<uint64_t> chunks;
            extractChunks(begin, end, chunks, indexes);
            assert(chunks.size() >= 2); // If all were the same, we would have not cut off enough

            if (chunks.size() < CHUNK_DIRECT_RANK_STORING_THRESHOLD) {
                treeNode.directRankStoring = true;
                treeNode.numChildren = chunks.size();

                auto it = begin;
                size_t bucketSizePrefixTemp = 0;
                for (size_t chunk = 0; chunk < chunks.size(); chunk++) {
                    uint64_t chunkValue = chunks.at(chunk);
                    auto currentBucketBegin = it;
                    while (it != end &&
                        readChunk(it->c_str(), it->length(), indexes.begin(), indexes.end()) == chunkValue) {
                        ++it;
                    }
                    std::string key = "Chunk" + path + "/" + std::to_string(chunkValue);
                    retrieval.addInput(chunks.size(), key, chunk);
                    constructChild(currentBucketBegin, it,
                                   lcpsBegin + (currentBucketBegin - begin), lcpsBegin + (it - begin),
                                   offset + bucketSizePrefixTemp, indexes[0],
                                   path + "/" + std::to_string(chunk), layer);
                    bucketSizePrefixTemp += std::distance(currentBucketBegin, it);
                }
            } else {
                treeNode.directRankStoring = false;
                treeNode.bucketMapper = new ShortPgmBucketMapper(chunks.begin(), chunks.end());
                assert(treeNode.bucketMapper->numBuckets() >= 2);

                auto it = begin;
                auto currentBucketBegin = begin;
                size_t prev_bucket = 0;
                size_t bucketSizePrefixTemp = 0;

                treeNode.bucketMapper->bucketOf(chunks.begin(), chunks.end(), [&](auto chunks_it, size_t bucket) {
                    while (it != end
                        && readChunk(it->c_str(), it->length(), indexes.begin(), indexes.end()) < *chunks_it) {
                        ++it;
                    }
                    while (prev_bucket < bucket) {
                        constructChild(currentBucketBegin, it,
                                       lcpsBegin + (currentBucketBegin - begin), lcpsBegin + (it - begin),
                                       offset + bucketSizePrefixTemp, indexes[0],
                                       path + "/" + std::to_string(prev_bucket), layer);
                        bucketSizePrefixTemp += std::distance(currentBucketBegin, it);
                        currentBucketBegin = it;
                        prev_bucket++;
                    }
                });
                it = end;
                while (prev_bucket < treeNode.bucketMapper->numBuckets()) {
                    constructChild(currentBucketBegin, it,
                                   lcpsBegin + (currentBucketBegin - begin), lcpsBegin + (it - begin),
                                   offset + bucketSizePrefixTemp, indexes[0],
                                   path + "/" + std::to_string(prev_bucket), layer);
                    bucketSizePrefixTemp += std::distance(currentBucketBegin, it);
                    currentBucketBegin = it;
                    prev_bucket++;
                }
            }
            bucketOffsetsInput.at(layer).push_back(offset + nThisNode);
            treeNodes.insert(std::make_pair(path, treeNode));
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

        void constructChild(auto begin, auto end, auto lcpsBegin, auto lcpsEnd, size_t offset,
                            size_t minLCP, std::string path, size_t layer) {
            bucketOffsetsInput.at(layer).push_back(offset);
            size_t currentBucketSize = std::distance(begin, end);
            if (currentBucketSize <= 1) {
                // Nothing to do
            } else if (currentBucketSize < DIRECT_RANK_STORING_THRESHOLD) {
                // Perform direct rank storing
                uint32_t indexInBucket = 0;
                auto it = begin;
                while (it != end) {
                    retrieval.addInput(currentBucketSize, util::MurmurHash64(*it), indexInBucket);
                    it++;
                    indexInBucket++;
                }
            } else {
                // Recurse
                constructNode(begin, end, lcpsBegin, lcpsEnd, minLCP, offset, path, layer + 1);
            }
        }

    public:
        ~RecursiveDirectRankStoringV2Mmphf() {
            for (auto fano : bucketOffsets) {
                delete fano;
            }
            for (auto &node : treeNodes) {
                if (!node.second.directRankStoring) {
                    delete node.second.bucketMapper;
                }
            }
        }

        static std::string name() {
            return "RecursiveDirectRankStoringV2Mmphf";
        }

        size_t spaceBits() {
            std::cout<<"Retrieval:           "<<1.0*retrieval.spaceBits()/N<<std::endl;
            std::cout<<"Bucket mapper:       "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                                                    [] (size_t size, auto &node) { return size + (node.second.directRankStoring ? 0 : node.second.bucketMapper->size()); }) / N<<std::endl;
            std::cout<<"Bucket offsets:      "<<1.0 * std::accumulate(bucketOffsets.begin(), bucketOffsets.end(), 0,
                                 [] (size_t size, PartitionedEliasFano *fano) { return size + fano->bit_size(); }) / N<<std::endl;
            std::cout<<"Tree node data:      "<<(8.0 * treeNodes.size() * sizeof(TreeNode) + 3.0 * treeNodes.size())/N<<std::endl;
            std::cout<<"Height:              "<<height<<std::endl;

            return retrieval.spaceBits()
                        + 8 * (sizeof(*this)
                                + treeNodes.size() * sizeof(TreeNode)
                                + bucketOffsets.size() * sizeof(void*))
                        + std::accumulate(bucketOffsets.begin(), bucketOffsets.end(), 0,
                                  [] (size_t size, PartitionedEliasFano *fano) { return size + fano->bit_size(); })
                        + 8 * std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                              [] (size_t size, auto &node) { return size + (node.second.directRankStoring ? 0 : node.second.bucketMapper->size()); })
                        + 3 * treeNodes.size(); // Use an MPHF instead of std::unordered_map
        }

        void exportTreeStructure(std::ostream &os) {
            os<<"digraph {"<<std::endl;
            exportTreeStructureInternal(os, PATH_ROOT, 0);
            os<<"}"<<std::endl;
        }

        void exportTreeStructureInternal(std::ostream &os, std::string path, size_t layer) {
            float scaleY = 200;
            float scaleX = 1.6 * (bucketOffsets.size() * scaleY) / N;

            TreeNode node = treeNodes.at(path);
            size_t beginX = bucketOffsets.at(layer)->at(node.offsetsOffset);
            size_t endX = bucketOffsets.at(layer)->at(node.offsetsOffset + node.bucketMapper->numBuckets());
            os<<"  \""<<path<<"\" [ "<<std::endl;
            os<<"    pos = \""<<+scaleX*(beginX+endX)/2<<","<<scaleY*layer<<"\""<<std::endl;
            os<<"    layer = \""<<+layer<<"\""<<std::endl;
            os<<"    label = \""<<+(endX - beginX)<<"\""<<std::endl;
            os<<"  ]"<<std::endl;

            for (size_t i = 0; i < node.bucketMapper->numBuckets(); i++) {
                std::string childPath = path + "/" + std::to_string(i);
                os<<"  \""<<path<<"\" -> \""<<childPath<<"\""<<std::endl;
                size_t childBegin = bucketOffsets.at(layer)->at(node.offsetsOffset + i);
                size_t childSize = bucketOffsets.at(layer)->at(node.offsetsOffset + i + 1) - childBegin;
                if (childSize < DIRECT_RANK_STORING_THRESHOLD) {
                    os<<"  \""<<childPath<<"\" [ "<<std::endl;
                    os<<"    pos = \""<<+scaleX*childBegin<<","<<scaleY*(layer+1)<<"\""<<std::endl;
                    os<<"    layer = \""<<+(layer+1)<<"\""<<std::endl;
                    os<<"    label = \""<<+childSize<<"\""<<std::endl;
                    os<<"  ]"<<std::endl;
                } else {
                    exportTreeStructureInternal(os, childPath, layer + 1);
                }
            }
        }

        uint64_t operator ()(const std::string &string) {
            std::string path = PATH_ROOT;
            TreeNode node = treeNodes.at(path);
            size_t layer = 0;
            while (true) {
                auto indexes = node.getIndexes();
                uint64_t chunk = readChunk(string.c_str(), string.length(), indexes.begin(), indexes.end());
                size_t bucket;

                if (node.directRankStoring) {
                    std::string key = "Chunk" + path + "/" + std::to_string(chunk);
                    bucket = retrieval.query(node.numChildren, key);
                } else {
                    bucket = node.bucketMapper->bucketOf(chunk);
                }

                size_t bucketOffset = bucketOffsets.at(layer)->at(node.offsetsOffset + bucket);
                size_t nextBucketOffset = bucketOffsets.at(layer)->at(node.offsetsOffset + bucket + 1);
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
                    layer++;
                    path += "/" + std::to_string(bucket);
                    node = treeNodes.at(path);
                }
            }
        }
};
