#pragma once

#include "util.h"

/**
 * Converts the entire string to chunks.
 * Chunks of all positions are handled together.
 */
struct FullChunkingStrategy {
    using Mmphf = DirectRankStoringMmphf<PgmBucketMapper<1.0f>>;
    Mmphf *mmphf = nullptr;
    size_t maxLCP;
    std::unordered_set<uint64_t> chunks;

    static std::string name() {
        return "FullChunkingStrategy";
    }

    explicit FullChunkingStrategy(size_t maxLCP)
            : maxLCP(maxLCP) {
    }

    FullChunkingStrategy(FullChunkingStrategy && rhs) {
        maxLCP = rhs.maxLCP;
        mmphf = rhs.mmphf;
        rhs.mmphf = nullptr;
    }

    ~FullChunkingStrategy() {
        delete mmphf;
    }

    void extractChunks(std::string &string) {
        const char *str = string.c_str();
        size_t length = std::min(maxLCP + 1, string.length());
        for (size_t i = 0; i < length; i += 8) {
            chunks.insert(readChunk(str + i, length - i));
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
        size_t encodedChunkWidth = bytesNeeded(mmphf->N + 1);
        std::string newString;
        const char *str = string.c_str();
        size_t length = std::min(maxLCP + 1, string.length());
        for (size_t i = 0; i < length; i += 8) {
            size_t chunkIndex = mmphf->operator()(readChunk(str + i, length - i));
            appendChunkIndex(newString, chunkIndex, encodedChunkWidth);
        }
        return newString;
    }

    size_t spaceBits() {
        return 8*sizeof(*this) + mmphf->spaceBits();
    }
};
