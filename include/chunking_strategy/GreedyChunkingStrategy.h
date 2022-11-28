#pragma once

#include "chunking_strategy/support/util.h"

/**
 * Only converts the first 8 bytes to a compressible chunk.
 * All other bytes just stay the same.
 * Idea: Quickly "eat up" common prefixes without having full compression tables.
 */
struct GreedyChunkingStrategy {
    using Mmphf = DirectRankStoringMmphf<PgmBucketMapper>;
    Mmphf *mmphf = nullptr;
    std::unordered_set<uint64_t> chunks;
    size_t maxLCP;

    static std::string name() {
        return "GreedyChunkingStrategy";
    }

    explicit GreedyChunkingStrategy(size_t maxLCP)
            : maxLCP(maxLCP) {
    }

    GreedyChunkingStrategy(GreedyChunkingStrategy && rhs) noexcept {
        maxLCP = rhs.maxLCP;
        mmphf = rhs.mmphf;
        rhs.mmphf = nullptr;
    }

    ~GreedyChunkingStrategy() noexcept {
        delete mmphf;
    }

    void extractChunks(std::string &string) {
        chunks.insert(readChunk(string.c_str(), string.length()));
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
        size_t chunkIndex = mmphf->operator()(readChunk(string.c_str(), string.length()));
        builder.appendInt(chunkIndex + 1, BITS_NEEDED(mmphf->N + 1));
        if (string.length() >= 8) {
            return builder.toString() + string.substr(8);
        }
        return builder.toString();
    }

    size_t spaceBits() {
        return 8*sizeof(*this) + mmphf->spaceBits();
    }
};
