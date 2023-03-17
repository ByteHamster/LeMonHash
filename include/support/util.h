#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <bit>
#include <sdsl/bits.hpp>
#include <set>

namespace lemonhash {
struct LcpDetails {
    uint32_t lcp;
    char branchingCharacter0;
    char branchingCharacter1;
};

LcpDetails LCP(const auto &s1, const auto &s2) {
    size_t lcp = 0;
    size_t minLength = std::min(s1.length(), s2.length());
    const char* s1ptr = s1.data();
    const char* s2ptr = s2.data();
    while (lcp < minLength - 8 && *(uint64_t*)(&s1ptr[lcp]) == *(uint64_t*)(&s2ptr[lcp])) {
        lcp += 8;
    }
    if (lcp < minLength - 8) {
        lcp += std::countr_zero(*(uint64_t*)(&s1ptr[lcp]) ^ *(uint64_t*)(&s2ptr[lcp])) / 8;
    } else {
        while (lcp < minLength && s1ptr[lcp] == s2ptr[lcp]) {
            lcp++;
        }
    }
    char branchingCharacter0 = lcp < s1.length() ? s1ptr[lcp] : '\0';
    char branchingCharacter1 = lcp < s2.length() ? s2ptr[lcp] : '\0';
    return {static_cast<uint32_t>(lcp), branchingCharacter0, branchingCharacter1};
}

std::vector<LcpDetails> computeLCPs(const auto begin, const auto end) {
    std::vector<LcpDetails> lcps(std::distance(begin, end));
    for (auto it = begin + 1; it != end; ++it)
        lcps[std::distance(begin, it)] = LCP(*std::prev(it), *it);
    return lcps;
}

/** Finds k distinct elements >= lowerBound in the range [begin, end) and returns them in a sorted vector. */
auto distinctMinima(const auto begin, const auto end, size_t k, auto lowerBound) {
    // TODO: The algorithm below costs O((end-begin)log(k))
    std::set<size_t> s;

    auto it = begin;
    while (it != end && s.size() < k) {
        if ((*it).lcp >= lowerBound)
            s.insert((*it).lcp);
        ++it;
    }

    for (; it != end; ++it) {
        if ((*it).lcp >= lowerBound && (s.empty() || (*it).lcp < *std::prev(s.end())) && !s.contains((*it).lcp)) {
            if (!s.empty())
                s.erase(std::prev(s.end()));
            s.insert((*it).lcp);
        }
    }

    std::vector<size_t> result;
    result.reserve(s.size());
    std::copy(s.begin(), s.end(), std::back_inserter(result));
    return result;
}

/** Finds the next set bit after a given position. */
static uint64_t nextOne(size_t i, const uint64_t *data) {
    auto wordIdx = i / 64;
    auto word = data[wordIdx] & sdsl::bits::lo_unset[i % 64];
    word &= word - 1;
    while (word == 0)
        word = data[++wordIdx];
    return wordIdx * 64 + __builtin_ctzll(word);
}

/** Finds the previous set bit before a given position. */
static uint64_t prevOne(size_t i, const uint64_t *data) {
    auto wordIdx = i / 64;
    auto word = data[wordIdx] & sdsl::bits::lo_set[i % 64];
    while (word == 0)
        word = data[--wordIdx];
    return wordIdx * 64 + 63 - __builtin_clzll(word);
}
} // namespace lemonhash
