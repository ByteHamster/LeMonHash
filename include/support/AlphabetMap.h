#pragma once

#include <stdint.h>
#include <cassert>
#include <bit>
#include <string>
#include <numeric>
#include "util.h"

/** Stores a subset of the alphabet of 7-bit or 8-bit values (depending on whether ASCII is true). */
template<bool ASCII=false>
struct AlphabetMap {
        static constexpr uint8_t words = ASCII ? 2 : 4;
        uint64_t bitMap[words] = {0};

        void set(uint8_t c) {
            if constexpr (ASCII)
                assert(c < 128);
            bitMap[c / 64] |= uint64_t(1) << (c % 64);
        }

    public:

        AlphabetMap() = default;

        AlphabetMap(uint64_t other[words]) {
            for (uint8_t i = 0; i < words; i++)
                bitMap[i] = other[i];
        }

        /** Constructs from the alphabet of the strings in the given range. If terminator is true, then '\0' is part of
         * the alphabet. */
        AlphabetMap(const auto begin, const auto end, bool terminator = true) {
            if (terminator)
                set(0);
            for (auto it = begin; it != end; ++it)
                for (uint8_t c : *it)
                    set(c);
        }

        /** Constructs from the alphabet of the strings in the given range, where only a suffix of each string is
         * considered.
         * If branchingCharacters is true, then consider only the alphabet of branching characters in the compacted trie
         * of the (suffixes of the) strings in the given range, which must be sorted lexicographically.
         * If terminator is true, then '\0' is part of the alphabet. */
        AlphabetMap(const auto begin, const auto end, size_t fromIndex, bool branchingCharacters, bool terminator = true) {
            if (terminator)
                set(0);

            if (branchingCharacters) {
                for (auto it = begin + 1; it != end; ++it) {
                    auto lcp = LCP(std::string_view(*it).substr(fromIndex), std::string_view(*(it - 1)).substr(fromIndex));
                    uint8_t c0 = (*std::prev(it))[fromIndex + lcp];
                    uint8_t c1 = (*it)[fromIndex + lcp];
                    set(c0);
                    set(c1);
                }
            } else {
                for (auto it = begin; it != end; ++it)
                    for (uint8_t c: std::string_view(*it).substr(fromIndex))
                        set(c);
            }
        }

        /** Constructs from the alphabet of branching characters in the compacted trie of the strings in the given range,
         * which must be sorted lexicographically, by using a precomputed LCP array.
         * If terminator is true, then '\0' is part of the alphabet. */
        AlphabetMap(const auto begin, const auto end, const auto lcpsBegin, bool terminator = true) {
            if (terminator)
                set(0);

            auto itLcps = lcpsBegin + 1;
            for (auto it = begin + 1; it != end; ++it, ++itLcps) {
                auto lcp = *itLcps;
                uint8_t c0 = (*std::prev(it))[lcp];
                uint8_t c1 = (*it)[lcp];
                set(c0);
                set(c1);
            }
        }

        /** Returns the rank of a given character in the alphabet. */
        size_t rank(uint8_t c) const {
            if constexpr (ASCII)
                c = c < 128 ? c : 127;
            auto rank = 0;
            for (auto i = 0; i < c / 64; i++)
                rank += std::popcount(bitMap[i]);
            return rank + std::popcount(bitMap[c / 64] & ((1ull << (c % 64)) - 1));
        }

        /** Returns true iff the given character is in the alphabet. */
        bool contains(uint8_t c) const {
            if constexpr (ASCII) {
                if (c >= 128)
                    return false;
            }
            return bitMap[c / 64] & (1ull << (c % 64));
        }

        /** Returns true iff the given alphabet map is contained in this alphabet map. */
        bool contains(const AlphabetMap<ASCII> &other) const {
            for (auto i = 0; i < words; i++)
                if ((bitMap[i] & other.bitMap[i]) != other.bitMap[i])
                    return false;
            return true;
        }

        /** Returns true iff the alphabet is a subset of ASCII. */
        bool isAscii() const {
            if constexpr (ASCII)
                return true;
            return bitMap[2] == 0 && bitMap[3] == 0;
        }

        /** Returns the ASCII version of this alphabet map. */
        AlphabetMap<true> toAscii() const {
            assert(isAscii());
            uint64_t a[2] = {bitMap[0], bitMap[1]};
            return {a};
        }

        /** Returns the number of characters in the alphabet. */
        size_t size() const {
            return std::accumulate(bitMap, bitMap + words, 0, [](auto a, auto b) { return a + std::popcount(b); });
        }

        /** Returns the length of the longest string from this alphabet that can be fit in a 64-bit integer. */
        size_t length64() const {
            return size2length64(size());
        }

        /** Creates a uint64_t from a prefix of the given string. */
        uint64_t readChunk(const char *string, size_t stringLength) const {
            size_t sigma = size();
            size_t characters = size2length64(sigma);
            uint64_t chunk = 0;
            size_t i = 0;
            for (; i < characters && i < stringLength; i++)
                chunk = chunk * sigma + std::min<size_t>(rank(string[i]), sigma - 1);
            for (; i < characters; i++)
                chunk *= sigma;
            return chunk;
        }

        /** Creates a uint64_t from the characters at the given indexes of the given string. */
        uint64_t readChunk(const char *string, size_t stringLength, const auto &indexes) const {
            size_t sigma = size();
            size_t characters = size2length64(sigma);
            uint64_t chunk = 0;
            size_t i = 0;
            for (; i < characters && i < indexes.size() && indexes[i] < stringLength; i++)
                chunk = chunk * sigma + std::min<size_t>(rank(string[indexes[i]]), sigma - 1);
            for (; i < characters; i++)
                chunk *= sigma;
            return chunk;
        }

    private:

        size_t size2length64(size_t size) const {
            static const uint8_t lookup[257] = {
                64, 64, 64, 40, 32, 27, 24, 22, 21, 20, 19, 18, 17, 17, 16, 16, 16, 15, 15, 15, 14,
                14, 14, 14, 13, 13, 13, 13, 13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 11,
                11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 10,
                10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
                9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                8, 8, 8, 8, 8, 8};
            return lookup[size]; // return size_t(64 / std::log2(std::max<size_t>(2, size)));
            /* Lookup table generation code:
                std::string line;
                for (size_t alphabetSize = 0; alphabetSize <= 256; alphabetSize++) {
                    size_t formula = 64 / std::log2(std::max<size_t>(2, alphabetSize));
                    line += std::to_string(formula) + ",";
                    if (line.length() > 80) {
                        std::cout<<line<<std::endl;
                        line = "";
                    } else {
                        line += " ";
                    }
                }
                std::cout<<line<<std::endl;*/
        }
};
