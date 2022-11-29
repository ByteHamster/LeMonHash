#pragma once

#include <string>

class ChunkingStrategy {
    public:
        virtual void extractChunks(std::string &string) = 0;
        virtual void build() = 0;
        virtual std::string compress(std::string &string) const = 0;
        virtual size_t spaceBits() = 0;
};