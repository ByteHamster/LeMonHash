#pragma once

#include "LeMonHashVL.hpp"

namespace lemonhash {
class LeMonHashVLIndexed {
    private:
        using BucketMapperType = UnalignedPGMBucketMapper<>;

        static constexpr size_t DIRECT_RANK_STORING_THRESHOLD = 128;

        static constexpr size_t CHUNK_DIRECT_RANK_STORING_THRESHOLD = 128;

        static constexpr size_t ALPHABET_MAPS_THRESHOLD = 512;

        static constexpr size_t LOG2_ALPHABET_MAPS_COUNT = 18;

        #pragma pack(push, 1)
        struct TreeNode {
            static constexpr bool useIndirection = sizeof(BucketMapperType) > 8;
            static constexpr uint8_t offsetsOffsetBits = 47 - LOG2_ALPHABET_MAPS_COUNT;
            union {
                std::conditional_t<useIndirection, BucketMapperType *, BucketMapperType> bucketMapper;
                size_t numChildren = 0;
            };
            bool directRankStoring : 1 = false;
            size_t offsetsOffset : offsetsOffsetBits = 0;
            size_t alphabetMapIndex : LOG2_ALPHABET_MAPS_COUNT;
            uint64_t packedIndexes = 0;

            static constexpr uint8_t packBaseBits = 15;
            static uint8_t getPackDeltaBits(uint8_t maxIndexesCount) {
                return (sizeof(packedIndexes) * 8 - packBaseBits) / (maxIndexesCount - 1);
            }

            TreeNode() {}
            TreeNode(const TreeNode &other) = delete;

            TreeNode(TreeNode &&other) {
                directRankStoring = other.directRankStoring;
                offsetsOffset = other.offsetsOffset;
                alphabetMapIndex = other.alphabetMapIndex;
                packedIndexes = other.packedIndexes;
                if (directRankStoring)
                    numChildren = other.numChildren;
                else
                    bucketMapper = std::move(other.bucketMapper);
                other.numChildren = 0;
                other.directRankStoring = false;
                other.offsetsOffset = 0;
                other.alphabetMapIndex = 0;
                other.bucketMapper = {};
            }

            TreeNode& operator=(const TreeNode &) = delete;
            TreeNode& operator=(TreeNode &&) = delete;

            const BucketMapperType &getBucketMapper() const {
                assert(!directRankStoring);
                return maybe_deref(bucketMapper);
            }

            size_t space() const {
                if (directRankStoring)
                    return 0;
                auto alreadyAccounted = useIndirection ? 0 : sizeof(BucketMapperType);
                return getBucketMapper().space() - alreadyAccounted;
            }

            template<typename RandomIt>
            void buildBucketMapper(RandomIt begin, RandomIt end) {
                destroyBucketMapper();
                directRankStoring = false;
                maybe_new(bucketMapper, begin, end);
            }

            std::vector<uint32_t> storeIndexes(const auto &indexes, uint8_t maxIndexesCount) {
                if (indexes[0] >= 1ull << packBaseBits)
                    throw std::runtime_error("Value too large for the current packing scheme");

                const auto packDeltaBits = getPackDeltaBits(maxIndexesCount);
                const auto maxDeltaValue = (uint64_t(1) << packDeltaBits) - 1;

                // store indexes[0] in packBaseBits, and indexes[i]-indexes[i-1] (i>0) in packDeltaBits
                packedIndexes = 0;
                packedIndexes |= indexes[0];
                size_t previous = indexes[0];
                std::vector<uint32_t> storedIndexes = { static_cast<uint32_t>(indexes[0]) };
                storedIndexes.reserve(maxIndexesCount);
                for (size_t deltaSlot = 0, i = 1; deltaSlot < maxIndexesCount - 1u && i < indexes.size(); ++deltaSlot) {
                    if (indexes[i] - previous >= maxDeltaValue) { // indexes[i] will span this slot and the next one(s)
                        packedIndexes |= maxDeltaValue << (packBaseBits + deltaSlot * packDeltaBits);
                        previous += maxDeltaValue;
                    } else {
                        packedIndexes |= (indexes[i] - previous) << (packBaseBits + deltaSlot * packDeltaBits);
                        previous = indexes[i];
                        ++i;
                        storedIndexes.push_back(previous);
                    }
                }
                return storedIndexes;
            }

            std::vector<uint32_t> getIndexes(uint8_t maxIndexesCount) const {
                const auto packDeltaBits = getPackDeltaBits(maxIndexesCount);
                const auto maxDeltaValue = (uint64_t(1) << packDeltaBits) - 1;
                std::vector<uint32_t> result;
                result.reserve(maxIndexesCount);
                result.push_back(packedIndexes & ((1ull << packBaseBits) - 1));
                uint64_t data = packedIndexes >> packBaseBits;
                uint64_t accumulator = 0;
                for (size_t i = 1; i < maxIndexesCount; ++i) {
                    if (data == 0) {
                        if (accumulator > 0)
                            result.push_back(result.back() + accumulator);
                        break;
                    }
                    auto delta = data & ((1ull << packDeltaBits) - 1);
                    if (delta == maxDeltaValue) {
                        accumulator += maxDeltaValue;
                    } else {
                        result.push_back(result.back() + delta + accumulator);
                        accumulator = 0;
                    }
                    data >>= packDeltaBits;
                }
                return result;
            }

            void destroyBucketMapper() {
                if (!directRankStoring)
                    maybe_delete(bucketMapper);
            }

            ~TreeNode() { destroyBucketMapper();}
        };
        #pragma pack(pop)

        MultiRetrievalDataStructure<> retrieval;
        size_t N = 0;
        std::vector<std::pair<uint64_t, TreeNode>> treeNodesInput;
        std::vector<TreeNode> treeNodes;
        using Mphf = pthash::single_phf<pthash::murmurhash2_64, pthash::compact_compact, true>;
        Mphf treeNodesMphf;
        std::vector<DuplicateFilterRank<PartitionedEliasFano>> bucketOffsets;
        AlphabetMapsCollection alphabetMaps;

    public:
        explicit LeMonHashVLIndexed(const std::vector<std::string> &strings) {
            N = strings.size();
            {
                auto lcps = computeLCPs(strings.begin(), strings.end());
                constructNode(strings.begin(), strings.end(), lcps.begin(), 0, 0, TreePath(), 0, 0);
            }
            retrieval.build();

            for (auto &v : bucketOffsets) {
                v.complete();
            }
            alphabetMaps.shrinkToFit();

            std::vector<uint64_t> mphfInput;
            for (auto &pair : treeNodesInput) {
                mphfInput.emplace_back(pair.first);
            }

            if (treeNodesInput.size() == 1) {
                treeNodes.resize(1);
                std::construct_at(&treeNodes[0], std::move(treeNodesInput.front().second));
            } else {
                pthash::build_configuration config;
                config.c = 7.0;
                config.alpha = 0.99;
                config.num_threads = 1;
                config.minimal_output = true;
                config.verbose_output = false;
                treeNodesMphf.build_in_internal_memory(mphfInput.begin(), mphfInput.size(), config);

                treeNodes.resize(treeNodesInput.size());
                for (auto &pair : treeNodesInput) {
                    size_t index = treeNodesMphf(pair.first);
                    std::construct_at(&treeNodes.at(index), std::move(pair.second));
                }
            }
            treeNodesInput.clear();
            treeNodesInput.shrink_to_fit();

            //std::ofstream myfile("tree2.dot");
            //exportTreeStructure(myfile);
            //myfile.close();
        }
    private:
        void constructNode(const auto begin, const auto end, const auto lcpsBegin, const size_t knownCommonPrefixLength,
                           size_t offset, const TreePath path, size_t layer, size_t ancestorAlphabetMapIndex) {
            size_t nThisNode = std::distance(begin, end);

            if (bucketOffsets.size() <= layer)
                bucketOffsets.resize(layer + 1);
            if (bucketOffsets.at(layer).size() > size_t(1) << TreeNode::offsetsOffsetBits)
                throw std::runtime_error("Increase offsetsOffsetBits");

            TreeNode treeNode;
            treeNode.offsetsOffset = bucketOffsets.at(layer).size();

            if ((!alphabetMaps.isFullForBits(LOG2_ALPHABET_MAPS_COUNT) && nThisNode > ALPHABET_MAPS_THRESHOLD) || alphabetMaps.empty()) {
                AlphabetMap am(lcpsBegin, nThisNode);
                if (!alphabetMaps.empty() && am.length64() == alphabetMaps.length64(ancestorAlphabetMapIndex)) {
                    treeNode.alphabetMapIndex = ancestorAlphabetMapIndex;
                } else {
                    treeNode.alphabetMapIndex = alphabetMaps.addOrFindSimilar(am);
                }
            } else {
                treeNode.alphabetMapIndex = ancestorAlphabetMapIndex;
            }

            auto maxIndexesCount = alphabetMaps.length64(treeNode.alphabetMapIndex);
            auto minima = distinctMinima(lcpsBegin + 1, lcpsBegin + nThisNode, maxIndexesCount, knownCommonPrefixLength);
            auto indexes = treeNode.storeIndexes(minima, maxIndexesCount);
            assert(indexes == treeNode.getIndexes(maxIndexesCount));

            auto [chunks, chunksOffsets] = extractChunks(begin, end, indexes, treeNode.alphabetMapIndex);
            assert(chunks.size() >= 2); // If all were the same, we would have not cut off enough

            auto tryPgmMapper = chunks.size() > CHUNK_DIRECT_RANK_STORING_THRESHOLD;
            if (tryPgmMapper) {
                treeNode.buildBucketMapper(chunks.begin(), chunks.end());
                auto &mapper = treeNode.getBucketMapper();
                if (mapper.bucketOf(chunks.front()) == mapper.bucketOf(chunks.back())) {
                    treeNode.destroyBucketMapper();
                    tryPgmMapper = false;
                } else {
                    assert(treeNode.getBucketMapper().numBuckets() >= 2);
                    auto currentBucketBegin = begin;
                    auto currentBucketEnd = begin;
                    size_t prevBucket = 0;
                    size_t bucketSizePrefixTemp = 0;

                    treeNode.getBucketMapper().bucketOf(chunks.begin(), chunks.end(), [&](auto chunks_it, size_t bucket) {
                        currentBucketEnd = begin + chunksOffsets.at(chunks_it - chunks.begin());
                        while (prevBucket < bucket) {
                            constructChild(currentBucketBegin, currentBucketEnd,
                                           lcpsBegin + (currentBucketBegin - begin), offset + bucketSizePrefixTemp,
                                           indexes[0], path, prevBucket, layer, treeNode.alphabetMapIndex);
                            bucketSizePrefixTemp += std::distance(currentBucketBegin, currentBucketEnd);
                            currentBucketBegin = currentBucketEnd;
                            prevBucket++;
                        }
                    });
                    if (currentBucketBegin != end) {
                        constructChild(currentBucketBegin, end, lcpsBegin + (currentBucketBegin - begin),
                                       offset + bucketSizePrefixTemp, indexes[0], path, prevBucket, layer,
                                       treeNode.alphabetMapIndex);
                    }
                }
            }

            if (!tryPgmMapper) {
                treeNode.directRankStoring = true;
                treeNode.numChildren = chunks.size();

                auto currentBucketEnd = begin;
                size_t bucketSizePrefixTemp = 0;
                for (size_t chunk = 0; chunk < chunks.size(); chunk++) {
                    uint64_t chunkValue = chunks.at(chunk);
                    auto currentBucketBegin = currentBucketEnd;
                    currentBucketEnd = begin + chunksOffsets[chunk + 1];
                    uint64_t key = path.getChild(chunkValue).alternativeHash();
                    retrieval.addInput(chunks.size(), key, chunk);
                    constructChild(currentBucketBegin, currentBucketEnd, lcpsBegin + (currentBucketBegin - begin),
                                   offset + bucketSizePrefixTemp, indexes[0], path, chunk, layer,
                                   treeNode.alphabetMapIndex);
                    bucketSizePrefixTemp += std::distance(currentBucketBegin, currentBucketEnd);
                }
            }
            bucketOffsets.at(layer).append(offset + nThisNode);
            treeNodesInput.emplace_back(path.currentNodeHash(), std::move(treeNode));
        }

        uint64_t readChunk(const std::string &s, const auto &indexes, size_t alphabetMapIndex) const {
            return alphabetMaps.readChunk(alphabetMapIndex, s.c_str(), s.length(), indexes);
        }

        auto extractChunks(const auto begin, const auto end, const auto &indexes, size_t alphabetMapIndex) {
            std::vector<uint64_t> chunks;
            std::vector<uint32_t> chunksOffsets;
            uint64_t previousChunk = 0;
            auto it = begin;
            while (it != end) {
                uint64_t chunk = readChunk(*it, indexes, alphabetMapIndex);
                assert(chunk >= previousChunk);
                if (chunk != previousChunk || chunks.empty()) {
                    chunks.push_back(chunk);
                    chunksOffsets.push_back(std::distance(begin, it));
                    previousChunk = chunk;
                }
                it++;
            }
            chunksOffsets.push_back(std::distance(begin, end));
            return std::make_pair(chunks, chunksOffsets);
        }

        void constructChild(auto begin, auto end, auto lcpsBegin, size_t offset, size_t minLCP,
                            const TreePath &path, size_t bucket, size_t layer, size_t ancestorAlphabetMapIndex) {
            bucketOffsets.at(layer).append(offset);
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
                constructNode(begin, end, lcpsBegin, minLCP, offset, path.getChild(bucket),
                              layer + 1, ancestorAlphabetMapIndex);
            }
        }

    public:
        ~LeMonHashVLIndexed() = default;

        static std::string name() {
            return "LeMonHashVLIndexed";
        }

        size_t spaceBits() {
            std::cout<<"Retrieval:           "<<1.0*retrieval.spaceBits()/N<<std::endl;
            std::cout<<"Bucket mapper:       "<<8.0*std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                                                    [] (size_t size, auto &node) { return size + node.space(); }) / N<<std::endl;
            std::cout<<"Bucket offsets:      "<<1.0 * std::accumulate(bucketOffsets.begin(), bucketOffsets.end(), 0,
                                                                      [] (size_t size, auto &fano) { return size + fano.bit_size(); }) / N<<std::endl;
            std::cout<<"Tree node data:      "<<(8.0 * (treeNodes.size() * sizeof(TreeNode) + sizeof(Mphf))
                                                + treeNodesMphf.num_bits())/N<<std::endl;
            std::cout<<"Alphabet maps:       "<<8.0 * alphabetMaps.sizeInBytes() / N<<std::endl;
            std::cout<<"7-bit alphabets:     "<<alphabetMaps.size().first<<std::endl;
            std::cout<<"8-bit alphabets:     "<<alphabetMaps.size().second<<std::endl;
            std::cout<<"Height:              "<<bucketOffsets.size()<<std::endl;
            std::cout<<"Nodes:               "<<treeNodes.size()<<std::endl;

            return retrieval.spaceBits()
                        + 8 * (sizeof(*this)
                                + treeNodes.size() * sizeof(TreeNode)
                                + sizeof(Mphf)
                                + bucketOffsets.size() * sizeof(void*))
                        + std::accumulate(bucketOffsets.begin(), bucketOffsets.end(), 0,
                                          [] (size_t size, auto &fano) { return size + fano.bit_size(); })
                        + 8 * std::accumulate(treeNodes.begin(), treeNodes.end(), 0,
                                              [] (size_t size, auto &node) { return size + node.space(); })
                        + 8 * alphabetMaps.sizeInBytes()
                        + treeNodesMphf.num_bits();
        }

        void exportTreeStructure(std::ostream &os) {
            os<<"digraph {"<<std::endl;
            exportTreeStructureInternal(os, TreePath(), 0);
            os<<"}"<<std::endl;
        }

        void exportTreeStructureInternal(std::ostream &os, TreePath path, size_t layer) {
            float scaleY = 200;
            float scaleX = 1.6 * (bucketOffsets.size() * scaleY) / N;

            TreeNode &node = treeNodes.at(treeNodesMphf(path.currentNodeHash()));
            size_t beginX = bucketOffsets.at(layer).at(node.offsetsOffset);
            size_t nodeSize = node.directRankStoring ? node.numChildren : node.getBucketMapper().numBuckets();
            size_t endX = bucketOffsets.at(layer).at(node.offsetsOffset + nodeSize);
            os<<"  \""<<path.currentNodeHash()<<"\" [ "<<std::endl;
            os<<"    pos = \""<<+scaleX*(beginX+endX)/2<<","<<scaleY*layer<<"\""<<std::endl;
            os<<"    layer = \""<<+layer<<"\""<<std::endl;
            os<<"    label = \""<<+(endX - beginX)<<"\""<<std::endl;
            os<<"  ]"<<std::endl;

            for (size_t i = 0; i < nodeSize; i++) {
                TreePath childPath = path.getChild(i);
                os<<"  \""<<path.currentNodeHash()<<"\" -> \""<<childPath.currentNodeHash()<<"\""<<std::endl;
                size_t childBegin = bucketOffsets.at(layer).at(node.offsetsOffset + i);
                size_t childSize = bucketOffsets.at(layer).at(node.offsetsOffset + i + 1) - childBegin;
                if (childSize < DIRECT_RANK_STORING_THRESHOLD) {
                    os<<"  \""<<childPath.currentNodeHash()<<"\" [ "<<std::endl;
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
            TreePath path;
            TreeNode *node = &treeNodes.at(treeNodes.size() == 1 ? 0 : treeNodesMphf(path.currentNodeHash()));
            size_t layer = 0;
            while (true) {
                auto indexes = node->getIndexes(alphabetMaps.length64(node->alphabetMapIndex));
                uint64_t chunk = readChunk(string, indexes, node->alphabetMapIndex);
                size_t bucket;

                if (node->directRankStoring) {
                    uint64_t key = path.getChild(chunk).alternativeHash();
                    bucket = retrieval.query(node->numChildren, key);
                } else {
                    bucket = node->getBucketMapper().bucketOf(chunk);
                }

                size_t bucketOffset = bucketOffsets.at(layer).at(node->offsetsOffset + bucket);
                size_t nextBucketOffset = bucketOffsets.at(layer).at(node->offsetsOffset + bucket + 1);
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
                    path = path.getChild(bucket);
                    node = &treeNodes.at(treeNodesMphf(path.currentNodeHash()));
                }
            }
        }
};
} // namespace lemonhash
