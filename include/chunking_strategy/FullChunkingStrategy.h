#pragma once

#include "chunking_strategy/support/util.h"
#include "ChunkingStrategy.hpp"

/**
 * Converts the entire string to chunks.
 * Chunks of all positions are handled together.
 */
struct FullChunkingStrategy : public ChunkingStrategy {
    using Mmphf = DirectRankStoringMmphf<SuccinctPgmBucketMapper>;
    Mmphf *mmphf = nullptr;
    std::unordered_set<uint64_t> chunks;

    std::string name() {
        return "FullChunkingStrategy";
    }

    FullChunkingStrategy(size_t maxLCP, size_t chunkWidth)
            : ChunkingStrategy(maxLCP, chunkWidth) {
    }

    FullChunkingStrategy(FullChunkingStrategy && rhs) noexcept
            : ChunkingStrategy(rhs.maxLCP, rhs.chunkWidth) {
        maxLCP = rhs.maxLCP;
        mmphf = rhs.mmphf;
        chunkWidth = rhs.chunkWidth;
        rhs.mmphf = nullptr;
    }

    ~FullChunkingStrategy() {
        delete mmphf;
    }

    void extractChunks(const std::string &string) {
        const char *str = string.c_str();
        size_t length = std::min(maxLCP + 1, string.length());
        for (size_t i = 0; i < length; i += chunkWidth) {
            chunks.insert(readChunk(str + i, length - i, chunkWidth));
        }
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

    std::string compress(const std::string &string) const {
        StringBuilder builder;
        const char *str = string.c_str();
        size_t length = std::min(maxLCP + 1, string.length());
        for (size_t i = 0; i < length; i += chunkWidth) {
            size_t chunkIndex = mmphf->operator()(readChunk(str + i, length - i, chunkWidth));
            builder.appendInt(chunkIndex + 1, BITS_NEEDED(mmphf->N + 1));
        }
        return builder.toString();
    }

    size_t spaceBits() {
        return 8*sizeof(*this) + mmphf->spaceBits(false);
    }
};
