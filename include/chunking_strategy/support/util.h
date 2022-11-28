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

#define BITS_NEEDED(x) (64 - __builtin_clzll(x))

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
