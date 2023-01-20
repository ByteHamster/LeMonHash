#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include "bucket_mapping/SuccinctPgmBucketMapper.hpp"
#include "bucket_mapping/UnalignedPgmBucketMapper.hpp"
#include "DirectRankStoring.hpp"
#include "bucket_mapping/support/EliasFanoModified.hpp"
#include <MurmurHash64.h>
#include "support/PartitionedEliasFano.hpp"

using BucketMapperType = UnalignedPgmBucketMapper;

template <typename T>
T& maybe_deref(T &x) { return x; }

template <typename T>
T& maybe_deref(T *x) { return *x; }

template <typename T>
void maybe_delete(T &x) { x.~T(); }

template <typename T>
void maybe_delete(T *x) { delete x; }

template <typename T, typename ...Args>
void maybe_new(T *&x, Args&&... args) { x = new T(std::forward<Args>(args)...); }

template <typename T, typename ...Args>
void maybe_new(T &x, Args&&... args) { x = std::move(T(std::forward<Args>(args)...)); }

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
            static constexpr bool useIndirection = sizeof(BucketMapperType) > 8;
            union {
                std::conditional_t<useIndirection, BucketMapperType *, BucketMapperType> bucketMapper;
                size_t numChildren = 0;
            };
            size_t offsetsOffset : 48 = 0;
            size_t minLCP : 15 = ((1 << 15) - 1);
            bool directRankStoring : 1 = false;

            TreeNode() {}
            TreeNode(const TreeNode &other) = delete;
            TreeNode(TreeNode &&other) {
                offsetsOffset = other.offsetsOffset;
                minLCP = other.minLCP;
                directRankStoring = other.directRankStoring;
                if (directRankStoring)
                    numChildren = other.numChildren;
                else
                    bucketMapper = std::move(other.bucketMapper);
                other.numChildren = 0;
                other.offsetsOffset = 0;
                other.minLCP = ((1 << 15) - 1);
                other.directRankStoring = false;
            }
            TreeNode& operator=(const TreeNode &) = delete;
            TreeNode& operator=(TreeNode &&) = delete;

            const BucketMapperType &getBucketMapper() const {
                assert(!directRankStoring);
                return maybe_deref(bucketMapper);
            }

            size_t size() const {
                auto alreadyAccounted = useIndirection ? 0 : sizeof(BucketMapperType);
                return directRankStoring ? 0 : getBucketMapper().size() - alreadyAccounted;
            }

            template<typename RandomIt>
            void buildBucketMapper(RandomIt begin, RandomIt end,
                                   std::unordered_map<uint64_t, size_t> &chunkOccurrences) {
                if (!directRankStoring)
                    maybe_delete(bucketMapper);
                directRankStoring = false;
                maybe_new(bucketMapper, begin, end, chunkOccurrences);
            }

            ~TreeNode() {
                if (!directRankStoring)
                    maybe_delete(bucketMapper);
            }
        };

        MultiRetrievalDataStructure retrieval;
        size_t N = 0;
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
            std::unordered_map<uint64_t, size_t> chunkOccurrences;
            extractChunks(begin, end, chunks, chunkOccurrences, treeNode.minLCP);
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
                treeNode.buildBucketMapper(chunks.begin(), chunks.end(), chunkOccurrences);
                size_t numBuckets = treeNode.getBucketMapper().numBuckets();
                assert(numBuckets >= 2);

                auto it = begin;
                auto currentBucketBegin = begin;
                size_t prev_bucket = 0;
                size_t bucketSizePrefixTemp = 0;

                treeNode.getBucketMapper().bucketOf(chunks.begin(), chunks.end(), [&](auto chunks_it, size_t bucket) {
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
            treeNodes.emplace(path, std::move(treeNode));
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

        void extractChunks(const auto begin, const auto end, std::vector<uint64_t> &chunks,
                           std::unordered_map<uint64_t, size_t> &chunkOccurrences, size_t minLCP) {
            uint64_t previousChunk = 0;
            auto it = begin;
            while (it != end) {
                uint64_t chunk = readChunk((*it).c_str() + minLCP, (*it).length() - minLCP, 8);
                assert(chunk >= previousChunk);
                if (chunk != previousChunk || chunks.empty()) {
                    chunks.push_back(chunk);
                    previousChunk = chunk;
                }
                if (!chunkOccurrences.contains(chunk)) {
                    chunkOccurrences.insert(std::make_pair(chunk, 0));
                }
                chunkOccurrences.at(chunk)++;
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
        }

        static std::string name() {
            return "RecursiveDirectRankStoringMmphf";
        }

        size_t spaceBits() {
            std::cout<<"Retrieval:           "<<1.0*retrieval.spaceBits()/N<<std::endl;
            std::cout<<"Bucket mapper:       "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                                                    [] (size_t size, auto &node) { return size + node.second.size(); }) / N<<std::endl;
            std::cout<<"Bucket offsets:      "<<1.0 * std::accumulate(bucketOffsets.begin(), bucketOffsets.end(), 0,
                                                                      [] (size_t size, PartitionedEliasFano *fano) { return size + fano->bit_size(); }) / N<<std::endl;
            std::cout<<"Tree node data:      "<<(8.0 * treeNodes.size() * sizeof(TreeNode) + 3.0 * treeNodes.size())/N<<std::endl;
            std::cout<<"Height:              "<<bucketOffsets.size()<<std::endl;

            return retrieval.spaceBits()
                        + 8 * (sizeof(*this)
                            + treeNodes.size() * sizeof(TreeNode)
                            + bucketOffsets.size() * sizeof(void*))
                        + std::accumulate(bucketOffsets.begin(), bucketOffsets.end(), 0,
                                          [] (size_t size, PartitionedEliasFano *fano) { return size + fano->bit_size(); })
                        + 8 * std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                              [] (size_t size, auto &node) { return size + node.second.size(); })
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

            TreeNode &node = treeNodes.at(path);
            size_t beginX = bucketOffsets.at(layer)->at(node.offsetsOffset);
            size_t nodeSize = node.directRankStoring ? node.numChildren : node.getBucketMapper().numBuckets();
            size_t endX = bucketOffsets.at(layer)->at(node.offsetsOffset + nodeSize);
            os<<"  \""<<path<<"\" [ "<<std::endl;
            os<<"    pos = \""<<+scaleX*(beginX+endX)/2<<","<<scaleY*layer<<"\""<<std::endl;
            os<<"    layer = \""<<+layer<<"\""<<std::endl;
            os<<"    label = \""<<+(endX - beginX)<<"\""<<std::endl;
            os<<"  ]"<<std::endl;

            for (size_t i = 0; i < nodeSize; i++) {
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
            TreeNode *node = &treeNodes.at(path);
            size_t layer = 0;
            while (true) {
                uint64_t chunk = readChunk(string.c_str() + node->minLCP, string.length() - node->minLCP, 8);
                size_t bucket;
                if (node->directRankStoring) {
                    std::string key = "Chunk" + path + "/" + std::to_string(chunk);
                    bucket = retrieval.query(node->numChildren, key);
                } else {
                    bucket = node->getBucketMapper().bucketOf(chunk);
                }
                size_t bucketOffset = bucketOffsets.at(layer)->at(node->offsetsOffset + bucket);
                size_t nextBucketOffset = bucketOffsets.at(layer)->at(node->offsetsOffset + bucket + 1);
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
                    node = &treeNodes.at(path);
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
            static constexpr bool useIndirection = sizeof(BucketMapperType) > 8;

            union {
                std::conditional_t<useIndirection, BucketMapperType *, BucketMapperType> bucketMapper;
                size_t numChildren = 0;
            };
            bool directRankStoring : 1 = false;
            uint32_t offsetsOffset : 31 = 0;
            uint64_t indexes = 0;

            static constexpr uint8_t packBaseBits = 15;
            static constexpr uint8_t packDeltaBits = (sizeof(indexes) * 8 - packBaseBits) / 7;

            TreeNode() {}
            TreeNode(const TreeNode &other) = delete;

            TreeNode(TreeNode &&other) {
                directRankStoring = other.directRankStoring;
                offsetsOffset = other.offsetsOffset;
                indexes = other.indexes;
                if (directRankStoring)
                    numChildren = other.numChildren;
                else
                    bucketMapper = std::move(other.bucketMapper);
                other.numChildren = 0;
                other.directRankStoring = false;
                other.offsetsOffset = 0;
            }

            TreeNode& operator=(const TreeNode &) = delete;
            TreeNode& operator=(TreeNode &&) = delete;

            const BucketMapperType &getBucketMapper() const {
                assert(!directRankStoring);
                return maybe_deref(bucketMapper);
            }

            size_t size() const {
                auto alreadyAccounted = useIndirection ? 0 : sizeof(BucketMapperType);
                return directRankStoring ? 0 : getBucketMapper().size() - alreadyAccounted;
            }

            template<typename RandomIt>
            void buildBucketMapper(RandomIt begin, RandomIt end,
                                   std::unordered_map<uint64_t, size_t> &chunkOccurrences) {
                if (!directRankStoring)
                    maybe_delete(bucketMapper);
                directRankStoring = false;
                maybe_new(bucketMapper, begin, end, chunkOccurrences);
            }

            std::vector<uint32_t> storeIndexes(const auto &other) {
                if (other[0] > (1ull << packBaseBits))
                    throw std::runtime_error("Value too large for the current packing scheme");

                indexes = 0;
                indexes |= other[0];
                size_t previous = other[0];
                uint64_t maxDeltaValue = (1ull << packDeltaBits) - 1;
                std::vector<uint32_t> storedIndexes = { other[0] };
                storedIndexes.reserve(8);
                for (size_t packed = 1, i = 1; packed < 8 && i < other.size(); ++packed) {
                    if (other[i] - previous >= (1ull << packDeltaBits)) {
                        indexes |= maxDeltaValue << (packBaseBits + (packed - 1) * packDeltaBits);
                        previous += maxDeltaValue;
                        storedIndexes.push_back(previous);
                    } else {
                        indexes |= (other[i] - previous) << (packBaseBits + (packed - 1) * packDeltaBits);
                        previous = other[i];
                        ++i;
                        storedIndexes.push_back(previous);
                    }
                }
                return storedIndexes;
            }

            std::vector<uint32_t> getIndexes() const {
                std::vector<uint32_t> result;
                result.reserve(8);
                result.push_back(indexes & ((1ull << packBaseBits) - 1));
                for (size_t i = 1; i < 8; ++i) {
                    auto delta = (indexes >> (packBaseBits + (i - 1) * packDeltaBits)) & ((1ull << packDeltaBits) - 1);
                    if (delta == 0)
                        break;
                    result.push_back(result.back() + delta);
                }
                return result;
            }

            ~TreeNode() {
                if (!directRankStoring)
                    maybe_delete(bucketMapper);
            }
        };
        #pragma pack(pop)

        MultiRetrievalDataStructure retrieval;
        size_t N = 0;
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
            auto indexes = treeNode.storeIndexes(distinctMinima(lcpsBegin + 1, lcpsEnd, 8, knownCommonPrefixLength));
            if (bucketOffsetsInput.size() <= layer) {
                bucketOffsetsInput.resize(layer + 1);
            }
            treeNode.offsetsOffset = bucketOffsetsInput.at(layer).size();
            size_t nThisNode = std::distance(begin, end);

            std::vector<uint64_t> chunks;
            std::unordered_map<uint64_t, size_t> chunkOccurrences;
            extractChunks(begin, end, chunks, chunkOccurrences, indexes);
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
                treeNode.buildBucketMapper(chunks.begin(), chunks.end(), chunkOccurrences);
                size_t numBuckets = treeNode.getBucketMapper().numBuckets();
                assert(numBuckets >= 2);

                auto it = begin;
                auto currentBucketBegin = begin;
                size_t prev_bucket = 0;
                size_t bucketSizePrefixTemp = 0;

                treeNode.getBucketMapper().bucketOf(chunks.begin(), chunks.end(), [&](auto chunks_it, size_t bucket) {
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
                while (prev_bucket < numBuckets) {
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
            treeNodes.emplace(path, std::move(treeNode));
        }

        void extractChunks(const auto begin, const auto end, std::vector<uint64_t> &chunks,
                           std::unordered_map<uint64_t, size_t> &chunkOccurrences, const auto &indexes) {
            uint64_t previousChunk = 0;
            auto it = begin;
            while (it != end) {
                uint64_t chunk = readChunk(it->c_str(), it->length(), indexes.begin(), indexes.end());
                assert(chunk >= previousChunk);
                if (chunk != previousChunk || chunks.empty()) {
                    chunks.push_back(chunk);
                    previousChunk = chunk;
                }
                if (!chunkOccurrences.contains(chunk)) {
                    chunkOccurrences.insert(std::make_pair(chunk, 0));
                }
                chunkOccurrences.at(chunk)++;
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
        }

        static std::string name() {
            return "RecursiveDirectRankStoringV2Mmphf";
        }

        size_t spaceBits() {
            std::cout<<"Retrieval:           "<<1.0*retrieval.spaceBits()/N<<std::endl;
            std::cout<<"Bucket mapper:       "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                                                    [] (size_t size, auto &node) { return size + node.second.size(); }) / N<<std::endl;
            std::cout<<"Bucket offsets:      "<<1.0 * std::accumulate(bucketOffsets.begin(), bucketOffsets.end(), 0,
                                                                      [] (size_t size, PartitionedEliasFano *fano) { return size + fano->bit_size(); }) / N<<std::endl;
            std::cout<<"Tree node data:      "<<(8.0 * treeNodes.size() * sizeof(TreeNode) + 3.0 * treeNodes.size())/N<<std::endl;
            std::cout<<"Height:              "<<bucketOffsets.size()<<std::endl;

            return retrieval.spaceBits()
                        + 8 * (sizeof(*this)
                                + treeNodes.size() * sizeof(TreeNode)
                                + bucketOffsets.size() * sizeof(void*))
                        + std::accumulate(bucketOffsets.begin(), bucketOffsets.end(), 0,
                                          [] (size_t size, PartitionedEliasFano *fano) { return size + fano->bit_size(); })
                        + 8 * std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                              [] (size_t size, auto &node) { return size + node.second.size(); })
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

            TreeNode &node = treeNodes.at(path);
            size_t beginX = bucketOffsets.at(layer)->at(node.offsetsOffset);
            size_t nodeSize = node.directRankStoring ? node.numChildren : node.getBucketMapper().numBuckets();
            size_t endX = bucketOffsets.at(layer)->at(node.offsetsOffset + nodeSize);
            os<<"  \""<<path<<"\" [ "<<std::endl;
            os<<"    pos = \""<<+scaleX*(beginX+endX)/2<<","<<scaleY*layer<<"\""<<std::endl;
            os<<"    layer = \""<<+layer<<"\""<<std::endl;
            os<<"    label = \""<<+(endX - beginX)<<"\""<<std::endl;
            os<<"  ]"<<std::endl;

            for (size_t i = 0; i < nodeSize; i++) {
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
            TreeNode *node = &treeNodes.at(path);
            size_t layer = 0;
            while (true) {
                auto indexes = node->getIndexes();
                uint64_t chunk = readChunk(string.c_str(), string.length(), indexes.begin(), indexes.end());
                size_t bucket;

                if (node->directRankStoring) {
                    std::string key = "Chunk" + path + "/" + std::to_string(chunk);
                    bucket = retrieval.query(node->numChildren, key);
                } else {
                    bucket = node->getBucketMapper().bucketOf(chunk);
                }

                size_t bucketOffset = bucketOffsets.at(layer)->at(node->offsetsOffset + bucket);
                size_t nextBucketOffset = bucketOffsets.at(layer)->at(node->offsetsOffset + bucket + 1);
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
                    node = &treeNodes.at(path);
                }
            }
        }
};
