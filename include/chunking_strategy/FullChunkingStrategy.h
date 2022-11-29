#pragma once

#include "chunking_strategy/support/util.h"

/**
 * Converts the entire string to chunks.
 * Chunks of all positions are handled together.
 */
struct FullChunkingStrategy {
    using Mmphf = DirectRankStoringMmphf<SuccinctPgmBucketMapper>;
    Mmphf *mmphf = nullptr;
    size_t maxLCP;
    size_t chunkWidth;
    std::unordered_set<uint64_t> chunks;

    static std::string name() {
        return "FullChunkingStrategy";
    }

    FullChunkingStrategy(size_t maxLCP, size_t chunkWidth)
            : maxLCP(maxLCP), chunkWidth(chunkWidth) {
    }

    FullChunkingStrategy(FullChunkingStrategy && rhs) noexcept {
        maxLCP = rhs.maxLCP;
        mmphf = rhs.mmphf;
        chunkWidth = rhs.chunkWidth;
        rhs.mmphf = nullptr;
    }

    ~FullChunkingStrategy() {
        delete mmphf;
    }

    void extractChunks(std::string &string) {
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

    std::string compress(std::string &string) const {
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
