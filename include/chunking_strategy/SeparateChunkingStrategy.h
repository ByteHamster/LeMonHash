#pragma once

#include "chunking_strategy/support/util.h"

/**
 * Converts the chunks of the entire string.
 * Each chunk position is stored in its own retrieval data structure.
 */
struct SeparateChunkingStrategy {
    using Mmphf = DirectRankStoringMmphf<SuccinctPgmBucketMapper>;
    size_t maxLCP;
    size_t chunkWidth;
    std::vector<Mmphf *> mmphfs;
    std::vector<std::unordered_set<uint64_t>> chunks;

    static std::string name() {
        return "SeparateChunkingStrategy";
    }

    explicit SeparateChunkingStrategy(size_t maxLCP, size_t chunkWidth)
            : maxLCP(maxLCP), chunkWidth(chunkWidth) {
        mmphfs.resize(maxLCP / chunkWidth + 1, nullptr);
        chunks.resize(maxLCP / chunkWidth + 1);
    }

    SeparateChunkingStrategy(SeparateChunkingStrategy && rhs) noexcept {
        maxLCP = rhs.maxLCP;
        mmphfs = rhs.mmphfs;
        chunkWidth = rhs.chunkWidth;
        rhs.mmphfs.clear();
    }

    ~SeparateChunkingStrategy() {
        for (Mmphf *mmphf : mmphfs) {
            delete mmphf;
        }
    }

    void extractChunks(std::string &string) {
        const char *str = string.c_str();
        size_t length = std::min(maxLCP + 1, string.length());
        for (size_t i = 0; i < length; i += chunkWidth) {
            chunks.at(i / chunkWidth).insert(readChunk(str + i, length - i, chunkWidth));
        }
    }

    void build() {
        for (size_t i = 0; i < chunks.size(); i++) {
            if (chunks.at(i).empty()) {
                continue;
            }
            std::vector<uint64_t> sortedChunks;
            sortedChunks.insert(sortedChunks.end(), chunks.at(i).begin(), chunks.at(i).end());
            std::cout<<"Creating MMPHF of "<<sortedChunks.size()<<" different chunks."<<std::endl;
            std::sort(sortedChunks.begin(), sortedChunks.end());
            mmphfs.at(i) = new Mmphf(sortedChunks);
            chunks.at(i).clear();
            chunks.at(i).rehash(0);
        }
    }

    std::string compress(std::string &string) const {
        StringBuilder builder;
        const char *str = string.c_str();
        size_t length = std::min(maxLCP + 1, string.length());
        for (size_t i = 0; i < length; i += chunkWidth) {
            size_t chunkIndex = mmphfs.at(i / chunkWidth)->operator()(readChunk(str + i, length - i, chunkWidth));
            builder.appendInt(chunkIndex + 1, BITS_NEEDED(mmphfs.at(i / chunkWidth)->N + 1));
        }
        return builder.toString();
    }

    size_t spaceBits() {
        size_t size = 8 * sizeof(*this);
        for (size_t i = 0; i < mmphfs.size(); i += 1) {
            if (mmphfs.at(i) == nullptr) {
                continue;
            }
            size += mmphfs.at(i)->spaceBits(false);
        }
        return size;
    }
};
