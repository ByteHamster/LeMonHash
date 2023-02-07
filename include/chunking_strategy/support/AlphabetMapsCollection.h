#pragma once
#include <vector>
#include "AlphabetMap.h"

class AlphabetMapsCollection {
        std::vector<AlphabetMap<true>> maps7;
        std::vector<AlphabetMap<false>> maps8;
    public:

        AlphabetMapsCollection() = default;

        uint64_t pushBack(const AlphabetMap<false> &map) {
            if (map.isAscii()) {
                maps7.emplace_back(map.toAscii());
                return (maps7.size() - 1) << 1;
            }

            maps8.emplace_back(map);
            return ((maps8.size() - 1) << 1) | 1;
        }

        bool isFullForBits(uint8_t mapIndexBits) const {
            auto maxIndex = (1ull << (mapIndexBits - 1)) - 1;
            return maps7.size() >= maxIndex || maps8.size() >= maxIndex;
        }

        uint8_t length64(uint64_t mapIndex) const {
            if (mapIndex & 1)
                return maps8[mapIndex >> 1].length64();
            return maps7[mapIndex >> 1].length64();
        }

        uint64_t readChunk(uint64_t mapIndex, const char *string, size_t stringLength) const {
            if (mapIndex & 1)
                return maps8[mapIndex >> 1].readChunk(string, stringLength);
            return maps7[mapIndex >> 1].readChunk(string, stringLength);
        }

        uint64_t readChunk(uint64_t mapIndex, const char *string, size_t stringLength, const auto &indexes) const {
            if (mapIndex & 1)
                return maps8[mapIndex >> 1].readChunk(string, stringLength, indexes);
            return maps7[mapIndex >> 1].readChunk(string, stringLength, indexes);
        }

        std::pair<size_t, size_t> size() const {
            return {maps7.size(), maps8.size()};
        }

        bool empty() const {
            return maps7.empty() && maps8.empty();
        }

        size_t sizeInBytes() const {
            return maps7.size() * sizeof(maps7[0]) + maps8.size() * sizeof(maps8[0]);
        }

        void shrinkToFit() {
            maps7.shrink_to_fit();
            maps8.shrink_to_fit();
        }
};
