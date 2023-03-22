#pragma once

#include <string>
#include "support/sequence/EliasFanoModified.hpp"

namespace lemonhash {
/**
 * Like LinearBucketMapper, but uses multiple segments, each of the same width.
 */
 template <size_t keysPerSegment = 256>
struct SegmentedLinearBucketMapper {
    const size_t n;
    EliasFanoM eliasFano;

    template<typename RandomIt>
    SegmentedLinearBucketMapper(RandomIt begin, RandomIt end)
            : n(end - begin),
              eliasFano(n / keysPerSegment + 2, *std::prev(end)) {

        auto it = begin;
        while (it < end) {
            eliasFano.push_back(*it);
            it += keysPerSegment;
        }
        eliasFano.push_back(*std::prev(end));
        eliasFano.buildRankSelect();
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const {
        auto segmentPtr = eliasFano.predecessorPosition(key);
        size_t segment = segmentPtr.index();
        uint64_t offset = *segmentPtr;
        uint64_t nextOffset = UINT64_MAX;
        if (segment < eliasFano.size() - 1) {
            ++segmentPtr;
            nextOffset = *segmentPtr;
        }
        double slope = (float) (nextOffset - offset) / keysPerSegment;
        size_t rankInSegment = std::min<size_t>((double) (key - offset) / slope, keysPerSegment);
        size_t estimatedRank = segment * keysPerSegment + rankInSegment;
        return std::min(n - 1, estimatedRank);
    }

    template<typename Iterator, typename Func>
    void bucketOf(Iterator first, Iterator last, Func f) const {
        while (first != last) {
            f(first, bucketOf(*first));
            ++first;
        }
    }

    [[nodiscard]] size_t space() const {
        return sizeof(*this) + eliasFano.space();
    }

    [[nodiscard]] size_t numBuckets() const {
        return n;
    }

    static std::string name() {
        return "SegmentedLinearBucketMapper keysPerSegment=" + std::to_string(keysPerSegment);
    }

    [[nodiscard]] std::string info() const {
        return "numBuckets=" + std::to_string(n);
    }
};
} // namespace lemonhash
