#pragma once
#include <vector>
#include "AlphabetMap.h"

class AlphabetMapsCollection {
        std::vector<AlphabetMap<true>> maps7;
        std::vector<AlphabetMap<false>> maps8;
    public:
        static constexpr size_t MAP_INDEX_FULL_BYTE = 1;

        AlphabetMapsCollection() {
            uint64_t fullAscii[] = {~0ul, ~0ul, 0ul, 0ul};
            addOrFindSimilar(AlphabetMap<false>(fullAscii));
            uint64_t fullByte[] = {~0ul, ~0ul, ~0ul, ~0ul};
            addOrFindSimilar(AlphabetMap<false>(fullByte));
            uint64_t lettersAndNumbers[] = {0xAFFFEFFF00000000, 0x17FFFFFE87FFFFFF, 0ul, 0ul};
            addOrFindSimilar(AlphabetMap<false>(lettersAndNumbers));
            uint64_t lowercaseLettersAndNumbers[] = {0x3FFA00100000000, 0x7FFFFFE80000000, 0ul, 0ul};
            addOrFindSimilar(AlphabetMap<false>(lowercaseLettersAndNumbers));
            uint64_t uppercaseLettersAndNumbers[] = {0x3FFA00100000000, 0x87FFFFFE, 0ul, 0ul};
            addOrFindSimilar(AlphabetMap<false>(uppercaseLettersAndNumbers));
        }

        uint64_t addOrFindSimilar(const AlphabetMap<false> &map) {
            if (map.isAscii()) {
                AlphabetMap<true> asciiMap = map.toAscii();
                size_t length64 = asciiMap.length64();
                for (size_t i = 0; i < maps7.size(); i++) {
                    if (maps7[i].contains(asciiMap) && maps7[i].length64() == length64) {
                        return i << 1;
                    }
                }
                maps7.push_back(asciiMap);
                return (maps7.size() - 1) << 1;
            }

            size_t length64 = map.length64();
            for (size_t i = 0; i < maps8.size(); i++) {
                if (maps8[i].contains(map) && maps8[i].length64() == length64) {
                    return (i << 1) | 1;
                }
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
