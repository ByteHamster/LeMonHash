#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include "chunking_strategy/FullChunkingStrategy.h"
#include "chunking_strategy/GreedyChunkingStrategy.h"
#include "chunking_strategy/SeparateChunkingStrategy.h"

template <typename ChunkingStrategy>
class ChunkedStringMmphf {
    private:
        std::vector<ChunkingStrategy> chunkingLayers;
        using Mmphf = DirectRankStoringMmphf<SuccinctPgmBucketMapper>;
        Mmphf *mmphf = nullptr;
        size_t N;
    public:
        explicit ChunkedStringMmphf(std::vector<std::string> strings)
                : N(strings.size()) {
            std::cout<<"Calculating LCPs"<<std::endl;
            size_t maxLCP = 0;
            for (size_t i = 1; i < strings.size(); i++) {
                maxLCP = std::max(maxLCP, LCP(strings.at(i), strings.at(i - 1)));
                assert(strings.at(i-1) < strings.at(i) && "Input needs to be sorted and unique");
            }

            while (maxLCP >= 8) {
                std::cout<<"Generating chunking layer "<<(chunkingLayers.size()+1)<<std::endl;
                chunkingLayers.emplace_back(maxLCP);
                ChunkingStrategy &chunkingLayer = chunkingLayers.back();
                for (std::string &string: strings) {
                    chunkingLayer.extractChunks(string);
                }
                chunkingLayer.build();
                maxLCP = 0;
                for (size_t i = 0; i < strings.size(); i++) {
                    std::string compressed = chunkingLayer.compress(strings.at(i));
                    //std::cout<<"Compressed string"<<toHexString(strings.at(i))<<" ----->"<<toHexString(compressed)<<std::endl;
                    strings.at(i) = compressed;
                    if (i > 0) {
                        maxLCP = std::max(maxLCP, LCP(strings.at(i), strings.at(i - 1)));
                        assert(strings.at(i-1) < strings.at(i));
                    }
                }
                //std::cout<<"------------------------"<<std::endl;
            }

            std::cout<<"Generating last level MMPHF"<<std::endl;
            std::vector<uint64_t> sortedChunks;
            for (size_t i = 0; i < strings.size(); i++) {
                sortedChunks.push_back(readChunk(strings.at(i).c_str(), strings.at(i).length()));
                assert((i == 0 || sortedChunks.at(i - 1) < sortedChunks.at(i)) && "All last-level strings must be unique and sorted");
            }
            mmphf = new Mmphf(sortedChunks);
        }

        ~ChunkedStringMmphf() {
            delete mmphf;
        }

        static std::string name() {
            return "ChunkedStringMmphf chunkingStrategy=" + ChunkingStrategy::name();
        }

        size_t spaceBits() {
            size_t bits = 8 * sizeof(*this);
            for (size_t i = 0; i < chunkingLayers.size(); i++) {
                std::cout<<"--------------- Chunking layer "<<i<<": ---------------"<<std::endl;
                size_t layerBits = chunkingLayers.at(i).spaceBits();
                bits += layerBits;
                std::cout<<"Layer total: "<<(1.0*layerBits/N)<<std::endl;
            }
            std::cout<<"--------------- Final layer: ---------------"<<std::endl;
            size_t layerBits = mmphf->spaceBits();
            bits += layerBits;
            std::cout<<"Layer total: "<<(1.0*layerBits/N)<<std::endl;
            return bits;
        }

        uint64_t operator ()(std::string string) {
            for (ChunkingStrategy &layer : chunkingLayers) {
                string = layer.compress(string);
            }
            return mmphf->operator()(readChunk(string.c_str(), string.length()));
        }
};