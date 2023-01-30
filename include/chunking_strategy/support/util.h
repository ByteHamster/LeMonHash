#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <bit>
#include "StringBuilder.hpp"

std::string toHexString(const char *string, size_t length) {
    char* str = new char[3*length+1];
    for (size_t i = 0; i < length; i++) {
        std::snprintf(str + i * 3, 4, " %02x", (unsigned char) string[i]);
    }
    std::string cppstr(str);
    delete[] str;
    return cppstr;
}

std::string toHexString(std::string &string) {
    return toHexString(string.c_str(), string.length());
}

std::string toBinaryString(const char *string, size_t length) {
    char* str = new char[9*length+1];
    for (size_t i = 0; i < length; i++) {
        std::snprintf(str + i * 9, 10, " %c%c%c%c%c%c%c%c",
                  (string[i] & 0x80 ? '1' : '0'),
                  (string[i] & 0x40 ? '1' : '0'),
                  (string[i] & 0x20 ? '1' : '0'),
                  (string[i] & 0x10 ? '1' : '0'),
                  (string[i] & 0x08 ? '1' : '0'),
                  (string[i] & 0x04 ? '1' : '0'),
                  (string[i] & 0x02 ? '1' : '0'),
                  (string[i] & 0x01 ? '1' : '0'));
    }
    std::string cppstr(str);
    delete[] str;
    return cppstr;
}

std::string toBinaryString(std::string &string) {
    return toBinaryString(string.c_str(), string.length());
}

size_t LCP(const auto &s1, const auto &s2) {
    size_t lcp = 0;
    size_t minLength = std::min(s1.length(), s2.length());
    auto s1ptr = s1.data();
    auto s2ptr = s2.data();
    while (lcp < minLength && s1ptr[lcp] == s2ptr[lcp]) {
        lcp++;
    }
    return lcp;
}

std::vector<uint32_t> computeLCPs(const auto begin, const auto end) {
    std::vector<uint32_t> lcps(std::distance(begin, end));
    for (auto it = begin + 1; it != end; ++it)
        lcps[std::distance(begin, it)] = LCP(*it, *std::prev(it));
    return lcps;
}

/** Finds k distinct elements >= lowerBound in the range [begin, end) and returns them in a sorted vector. */
auto distinctMinima(const auto begin, const auto end, size_t k, auto lowerBound) {
    // TODO: The algorithm below costs O((end-begin)log(k))
    using value_type = typename decltype(begin)::value_type;
    std::set<value_type> s;

    auto it = begin;
    while (it != end && s.size() < k) {
        if (*it >= lowerBound)
            s.insert(*it);
        ++it;
    }

    for (; it != end; ++it) {
        if (*it >= lowerBound && (s.empty() || *it < *std::prev(s.end())) && !s.contains(*it)) {
            if (!s.empty())
                s.erase(std::prev(s.end()));
            s.insert(*it);
        }
    }

    std::vector<value_type> result;
    result.reserve(s.size());
    std::copy(s.begin(), s.end(), std::back_inserter(result));
    return result;
}

#define BITS_NEEDED(x) (64 - __builtin_clzll(x))

/**
 * Reads an 8-byte prefix of a string to a uint64_t.
 * By reversing the bytes, this ensures that sorting the integer
 * has the same effect as sorting the string lexicograpically.
 */
uint64_t readChunk(const char *string, size_t maxLength, size_t chunkWidth) {
    assert(chunkWidth <= 8);
    assert(chunkWidth >= 1);
    if (chunkWidth == 8 && maxLength >= 8)
        return __builtin_bswap64(*((uint64_t*) string));
    uint64_t chunk = 0;
    char *chunkRaw = (char*) &chunk;
    for (size_t i = 0; i < chunkWidth && i < maxLength; i++) {
        chunkRaw[chunkWidth - 1 - i] = string[i];
    }
    return chunk;
}

/** Creates a uint64_t with the characters at the given indexes of the given string. */
template<typename Iterator>
uint64_t readChunk(const char *string, size_t stringLength, Iterator indexesBegin, Iterator indexesEnd) {
    assert(indexesEnd - indexesBegin <= 8);
    assert(indexesEnd - indexesBegin >= 1);
    uint64_t chunk = 0;
    auto chunkRaw = (char*) &chunk;
    auto chunkWidth = std::distance(indexesBegin, indexesEnd);
    for (auto it = indexesBegin; it != indexesEnd && *it < stringLength; ++it)
        chunkRaw[chunkWidth - 1 - (it - indexesBegin)] = string[*it];
    return chunk;
}

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
    uint8_t rank(uint8_t c) const {
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
    uint8_t size() const {
        return std::accumulate(bitMap, bitMap + words, 0, [](auto a, auto b) { return a + std::popcount(b); });
    }

    uint8_t size2length64(uint8_t size) const {
        static const uint8_t lookup[256] = {
            0, 0, 63, 39, 31, 27, 24, 22, 21, 19, 18, 18, 17, 17, 16, 16, 15, 15, 15, 14, 14, 14, 14, 13, 13, 13, 13,
            13, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9,
            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
            8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
            8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
            8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};
        return lookup[size];
        //return uint8_t(63 / std::log2(size));
    }

    /** Returns the length of the longest string from this alphabet that can be fit in a 64-bit integer. */
    uint8_t length64() const {
        return size2length64(size());
    }

    /** Creates a uint64_t from a prefix of the given string. */
    uint64_t readChunk(const char *string, size_t stringLength) const {
        auto sigma = size();
        auto characters = size2length64(sigma);
        uint64_t chunk = 0;
        size_t i = 0;
        for (; i < characters && i < stringLength; i++)
            chunk = chunk * sigma + std::min<uint8_t>(rank(string[i]), sigma - 1);
        for (; i < characters; i++)
            chunk *= sigma;
        return chunk;
    }

    /** Creates a uint64_t from the characters at the given indexes of the given string. */
    uint64_t readChunk(const char *string, size_t stringLength, const auto &indexes) const {
        auto sigma = size();
        auto characters = size2length64(sigma);
        uint64_t chunk = 0;
        size_t i = 0;
        for (; i < characters && i < indexes.size() && indexes[i] < stringLength; i++)
            chunk = chunk * sigma + std::min<uint8_t>(rank(string[indexes[i]]), sigma - 1);
        for (; i < characters; i++)
            chunk *= sigma;
        return chunk;
    }
};

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
