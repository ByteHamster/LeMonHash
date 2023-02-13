#pragma once

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <bit>

size_t LCP(const auto &s1, const auto &s2) {
    size_t lcp = 0;
    size_t minLength = std::min(s1.length(), s2.length());
    const char* s1ptr = s1.data();
    const char* s2ptr = s2.data();
    while (lcp < minLength - 8 && *(uint64_t*)(&s1ptr[lcp]) == *(uint64_t*)(&s2ptr[lcp])) {
        lcp += 8;
    }
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
