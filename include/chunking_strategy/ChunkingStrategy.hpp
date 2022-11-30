#pragma once

#include <string>

class ChunkingStrategy {
    public:
        size_t maxLCP;
        size_t chunkWidth;

        ChunkingStrategy(size_t maxLCP, size_t chunkWidth)
            : maxLCP(maxLCP), chunkWidth(chunkWidth) {
        }

        virtual void extractChunks(std::string &string) = 0;
        virtual void build() = 0;
        virtual std::string compress(std::string &string) const = 0;
        virtual size_t spaceBits() = 0;
        virtual std::string name() = 0;
};