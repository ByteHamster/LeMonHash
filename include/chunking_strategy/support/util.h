#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
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

size_t LCP(const std::string &s1, const std::string &s2) {
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