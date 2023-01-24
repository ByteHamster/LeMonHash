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

/** Stores a subset of the alphabet of 1-byte values [0, 1, ..., 255]. */
class AlphabetMap {
    uint64_t bitMap[4] = {0, 0, 0, 0};

public:

    AlphabetMap() = default;

    /** Constructs from the alphabet of the strings in the given range. If terminator is true, then '\0' is part of
     * the alphabet. */
    AlphabetMap(const auto begin, const auto end, bool terminator = true) {
        if (terminator)
            bitMap[0] = 1;
        for (auto it = begin; it != end; ++it)
            for (uint8_t c : *it)
                bitMap[c / 64] |= 1ull << (c % 64);
    }

    /** Constructs from the alphabet of the strings in the given range, where only a suffix of each string is
     * considered.
     * If branchingCharacters is true, then consider only the alphabet of branching characters in the compacted trie
     * of the (suffixes of the) strings in the given range, which must be sorted lexicographically.
     * If terminator is true, then '\0' is part of the alphabet. */
    AlphabetMap(const auto begin, const auto end, size_t fromIndex, bool branchingCharacters, bool terminator = true) {
        if (terminator)
            bitMap[0] = 1;

        if (branchingCharacters) {
            for (auto it = begin + 1; it != end; ++it) {
                auto lcp = LCP(std::string_view(*it).substr(fromIndex), std::string_view(*(it - 1)).substr(fromIndex));
                uint8_t c0 = (*std::prev(it))[fromIndex + lcp];
                uint8_t c1 = (*it)[fromIndex + lcp];
                bitMap[c0 / 64] |= 1ull << (c0 % 64);
                bitMap[c1 / 64] |= 1ull << (c1 % 64);
            }
        } else {
            for (auto it = begin; it != end; ++it)
                for (uint8_t c: std::string_view(*it).substr(fromIndex))
                    bitMap[c / 64] |= 1ull << (c % 64);
        }
    }

    /** Returns the rank of a given character in the alphabet. */
    uint8_t rank(uint8_t c) const {
        auto rank = 0;
        for (auto i = 0; i < c / 64; i++)
            rank += std::popcount(bitMap[i]);
        return rank + std::popcount(bitMap[c / 64] & ((1ull << (c % 64)) - 1));
    }

    /** Returns true iff the given character is in the alphabet. */
    bool contains(uint8_t c) const {
        return bitMap[c / 64] & (1ull << (c % 64));
    }

    /** Returns true iff the given alphabet map is contained in this alphabet map. */
    bool contains(const AlphabetMap &other) const {
        for (auto i = 0; i < 4; i++)
            if ((bitMap[i] & other.bitMap[i]) != other.bitMap[i])
                return false;
        return true;
    }

    /** Returns the number of characters in the alphabet. */
    uint8_t size() const {
        return std::accumulate(bitMap, bitMap + 4, 0, [](auto a, auto b) { return a + std::popcount(b); });
    }

    /** Returns the length of the longest string from this alphabet that can be fit in a 64-bit integer. */
    uint8_t length64() const {
        return uint8_t(63 / std::log2(size()));
    }

    /** Creates a uint64_t from a prefix of the given string. */
    uint64_t readChunk(const char *string, size_t stringLength) const {
        auto sigma = size();
        auto characters = length64();
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
        auto characters = length64();
        uint64_t chunk = 0;
        size_t i = 0;
        for (; i < characters && i < indexes.size() && indexes[i] < stringLength; i++)
            chunk = chunk * sigma + std::min<uint8_t>(rank(string[indexes[i]]), sigma - 1);
        for (; i < characters; i++)
            chunk *= sigma;
        return chunk;
    }
};
