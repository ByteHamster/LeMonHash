#pragma once

#include "support/PolymorphicPGM.hpp"
#include "support/util.hpp"

/**
 * Uses a PGM-index to get an approximate rank, which we use as bucket index.
 */
struct PolymorphicPgmBucketMapper {
    pgm::PolymorphicPGMIndex pgm;

    PolymorphicPgmBucketMapper() = default;

    template<typename RandomIt>
    PolymorphicPgmBucketMapper(RandomIt begin, RandomIt end) : pgm() {
        auto bestCost = std::numeric_limits<size_t>::max();

        size_t cost;
        size_t previousBucket;
        RandomIt bucketBegin;
        auto updateCost = [&](auto it, size_t bucket) {
            if (bucket != previousBucket) {
                auto bucketSize = std::distance(bucketBegin, it);
                cost += bucketSize <= 1 ? 0 : bucketSize * BIT_WIDTH(bucketSize - 1);
                previousBucket = bucket;
                bucketBegin = it;
            }
        };

        for (auto epsilon : {3, 7, 15, 31, 63}) {
            pgm::PolymorphicPGMIndex p(begin, end, epsilon);

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

    [[nodiscard]] size_t segments_count() const {
        return pgm.segments_count();
    }

    [[nodiscard]] size_t numBuckets() const {
        return pgm.size();
    }

    static std::string name() {
        return std::string("PolymorphicPgmBucketMapper");
    }

};