#pragma once

#include "chunking_strategy/support/util.h"

/**
 * Converts the chunks of the entire string.
 * Each chunk position is stored in its own retrieval data structure.
 */
struct SeparateChunkingStrategy {
    using Mmphf = DirectRankStoringMmphf<PgmBucketMapper<1.0f>>;
    size_t maxLCP;
    std::vector<Mmphf *> mmphfs;
    std::vector<std::unordered_set<uint64_t>> chunks;

    static std::string name() {
        return "SeparateChunkingStrategy";
    }

    explicit SeparateChunkingStrategy(size_t maxLCP)
            : maxLCP(maxLCP) {
        mmphfs.resize(maxLCP / 8 + 1, nullptr);
        chunks.resize(maxLCP / 8 + 1);
    }

    SeparateChunkingStrategy(SeparateChunkingStrategy && rhs) {
        maxLCP = rhs.maxLCP;
        mmphfs = rhs.mmphfs;
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
        for (size_t i = 0; i < length; i += 8) {
            chunks.at(i / 8ul).insert(readChunk(str + i, length - i));
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
        std::string newString;
        const char *str = string.c_str();
        size_t length = std::min(maxLCP + 1, string.length());
        size_t i = 0;
        for (; i < length; i += 8) {
            size_t encodedChunkWidth = bytesNeeded(mmphfs.at(i / 8ul)->N + 1);
            size_t chunkIndex = mmphfs.at(i / 8ul)->operator()(readChunk(str + i, length - i));
            appendChunkIndex(newString, chunkIndex, encodedChunkWidth);
        }
        if (i < length) {
            newString.append(string.data() + i, length - i);
        }
        return newString;
    }

    size_t spaceBits() {
        size_t size = 8 * sizeof(*this);
        for (size_t i = 0; i < mmphfs.size(); i += 1) {
            if (mmphfs.at(i) == nullptr) {
                continue;
            }
            size += mmphfs.at(i)->spaceBits();
        }
        return size;
    }
};
