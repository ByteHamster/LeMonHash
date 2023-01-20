#pragma once

#include "support/UnalignedPGM.hpp"
#include "support/util.hpp"

/**
 * Uses a PGM-index to get an approximate rank, which we use as bucket index.
 */
struct UnalignedPgmBucketMapper {
    pgm::UnalignedPGMIndex pgm;

    UnalignedPgmBucketMapper() = default;

    template<typename RandomIt>
    UnalignedPgmBucketMapper(RandomIt begin, RandomIt end,
                             std::unordered_map<uint64_t, size_t> &chunkOccurrences) : pgm() {
        auto bestCost = std::numeric_limits<size_t>::max();

        size_t cost;
        size_t previousBucket;
        size_t bucketSize = 0;
        RandomIt bucketBegin;
        auto updateCost = [&](auto it, size_t bucket) {
            bucketSize += (it == end) ? 0 : chunkOccurrences.at(*it);
            if (bucket != previousBucket) {
                size_t bucketCostEstimate = bucketSize <= 1 ? 0 : bucketSize * BIT_WIDTH(bucketSize - 1);
                cost += std::min(bucketCostEstimate, 6 * bucketSize); // 6 * bucketSize is an estimate for recursion
                previousBucket = bucket;
                bucketBegin = it;
                bucketSize = 0;
            }
        };

        for (auto epsilon : {3, 7, 15, 31, 63}) {
            pgm::UnalignedPGMIndex p(begin, end, epsilon);

            cost = p.size_in_bytes() * 8;
            previousBucket = 0;
            bucketBegin = begin;
            p.for_each(begin, end, updateCost);
            updateCost(end, std::numeric_limits<uint64_t>::max());

            auto segmentsCount = p.segments_count();
            if (cost < bestCost) {
                pgm = std::move(p);
                bestCost = cost;
            }
            if (segmentsCount == 1)
                break;
        }
        assert(bucketOf(*std::prev(end)) != bucketOf(*begin) && "Need at least 2 buckets, otherwise we get a loop");
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const  {
        return pgm.approximate_rank(key);
    }

    template<typename Iterator, typename Func>
    void bucketOf(Iterator first, Iterator last, Func f) const {
        pgm.for_each(first, last, f);
    }

    [[nodiscard]] size_t size() const {
        return pgm.size_in_bytes();
    }

    [[nodiscard]] size_t numBuckets() const {
        return pgm.size();
    }

    static std::string name() {
        return std::string("UnalignedPgmBucketMapper");
    }

};