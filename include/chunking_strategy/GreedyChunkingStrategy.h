#pragma once

#include "chunking_strategy/support/util.h"
#include "ChunkingStrategy.hpp"

/**
 * Only converts the first 8 bytes to a compressible chunk.
 * All other bytes just stay the same.
 * Idea: Quickly "eat up" common prefixes without having full compression tables.
 */
struct GreedyChunkingStrategy : public ChunkingStrategy {
    using Mmphf = DirectRankStoringMmphf<SuccinctPgmBucketMapper>;
    Mmphf *mmphf = nullptr;
    std::unordered_set<uint64_t> chunks;

    std::string name() {
        return "GreedyChunkingStrategy";
    }

    explicit GreedyChunkingStrategy(size_t maxLCP, size_t chunkWidth)
            : ChunkingStrategy(maxLCP, chunkWidth) {
    }

    GreedyChunkingStrategy(GreedyChunkingStrategy && rhs) noexcept
            : ChunkingStrategy(rhs.maxLCP, rhs.chunkWidth) {
        maxLCP = rhs.maxLCP;
        mmphf = rhs.mmphf;
        chunkWidth = rhs.chunkWidth;
        rhs.mmphf = nullptr;
    }

    ~GreedyChunkingStrategy() noexcept {
        delete mmphf;
    }

    void extractChunks(std::string &string) {
        chunks.insert(readChunk(string.c_str(), string.length(), chunkWidth));
    }

    void build() {
        std::vector<uint64_t> sortedChunks;
        sortedChunks.insert(sortedChunks.end(), chunks.begin(), chunks.end());
        std::cout<<"Creating MMPHF of "<<sortedChunks.size()<<" different chunks."<<std::endl;
        std::sort(sortedChunks.begin(), sortedChunks.end());
        assert(!sortedChunks.empty());
        mmphf = new Mmphf(sortedChunks);
        chunks.clear();
        chunks.rehash(0);
    }

    std::string compress(std::string &string) const {
        StringBuilder builder;
        size_t chunkIndex = mmphf->operator()(readChunk(string.c_str(), string.length(), chunkWidth));
        builder.appendInt(chunkIndex + 1, BITS_NEEDED(mmphf->N + 1));
        if (string.length() >= chunkWidth) {
            return builder.toString() + string.substr(chunkWidth);
        }
        return builder.toString();
    }

    size_t spaceBits() {
        return 8*sizeof(*this) + mmphf->spaceBits(false);
    }
};
