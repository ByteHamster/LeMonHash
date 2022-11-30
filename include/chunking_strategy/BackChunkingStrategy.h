#pragma once

#include "chunking_strategy/support/util.h"
#include "ChunkingStrategy.hpp"

/**
 * Keeps short strings as they are but compresses the back of long strings.
 */
struct BackChunkingStrategy : public ChunkingStrategy {
    using Mmphf = DirectRankStoringMmphf<SuccinctPgmBucketMapper>;
    Mmphf *mmphf = nullptr;
    std::unordered_set<uint64_t> chunks;

    std::string name() {
        return "BackChunkingStrategy";
    }

    BackChunkingStrategy(size_t maxLCP, size_t chunkWidth)
            : ChunkingStrategy(maxLCP, chunkWidth) {
    }

    BackChunkingStrategy(size_t maxLcp, size_t chunkWidth, BackChunkingStrategy &&rhs) noexcept
            : ChunkingStrategy(rhs.maxLCP, rhs.chunkWidth) {
        maxLCP = rhs.maxLCP;
        mmphf = rhs.mmphf;
        chunkWidth = rhs.chunkWidth;
        rhs.mmphf = nullptr;
    }

    ~BackChunkingStrategy() {
        delete mmphf;
    }

    void extractChunks(std::string &string) {
        const char *str = string.c_str();
        size_t length = std::min(maxLCP + 1, string.length());
        for (size_t i = maxLCP/2; i < length; i += chunkWidth) {
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
        for (size_t i = maxLCP/2; i < length; i += chunkWidth) {
            size_t chunkIndex = mmphf->operator()(readChunk(str + i, length - i, chunkWidth));
            builder.appendInt(chunkIndex + 1, BITS_NEEDED(mmphf->N + 1));
        }
        return string.substr(0, std::min(maxLCP/2, length)) + builder.toString();
    }

    size_t spaceBits() {
        return 8*sizeof(*this) + mmphf->spaceBits(false);
    }
};
