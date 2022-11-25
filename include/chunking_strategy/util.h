#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>

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

size_t LCP(std::string &s1, std::string &s2) {
    size_t lcp = 0;
    size_t minLength = std::min(s1.length(), s2.length());
    char* s1ptr = s1.data();
    char* s2ptr = s2.data();
    while (lcp < minLength && s1ptr[lcp] == s2ptr[lcp]) {
        lcp++;
    }
    return lcp;
}

size_t bytesNeeded(uint64_t i) {
    size_t bytesNeeded = 8;
    uint64_t mask = (~0ul) >> 8;
    while (i == (i & mask)) {
        bytesNeeded--;
        mask = mask >> 8;
    }
    return bytesNeeded;
}

/**
 * Reads an 8-byte prefix of a string to a uint64_t.
 * By reversing the bytes, this ensures that sorting the integer
 * has the same effect as sorting the string lexicograpically.
 */
uint64_t readChunk(const char *string, size_t maxLength) {
    uint64_t chunk = 0;
    char *chunkRaw = (char*) &chunk;
    for (size_t i = 0; i < 8 && i < maxLength; i++) {
        chunkRaw[7 - i] = string[i];
    }
    return chunk;
}

void appendChunkIndex(std::string &string, size_t chunkIndex, size_t width) {
    chunkIndex++; // String end (\0) must be smaller than all chunk indices
    for (size_t k = 0; k < width; k++) {
        string += ((char*) &chunkIndex)[width - k - 1];
    }
}
