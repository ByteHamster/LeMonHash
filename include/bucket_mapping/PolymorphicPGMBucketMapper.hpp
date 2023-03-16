#pragma once

#include "support/PolymorphicPGM.hpp"
#include <bit>

/**
 * Uses a PGM-index to get an approximate rank, which we use as bucket index.
 */
struct PolymorphicPGMBucketMapper {
    pgm::PolymorphicPGMIndex pgm;

    PolymorphicPGMBucketMapper() = default;

    template<typename RandomIt>
    PolymorphicPGMBucketMapper(RandomIt begin, RandomIt end) : pgm() {
        auto bestCost = std::numeric_limits<size_t>::max();

        size_t cost;
        size_t previousBucket;
        RandomIt bucketBegin;
        auto updateCost = [&](auto it, size_t bucket) {
            if (bucket != previousBucket) {
                auto bucketSize = (size_t) std::distance(bucketBegin, it);
                cost += bucketSize <= 1 ? 0 : bucketSize * std::bit_width(bucketSize - 1);
                previousBucket = bucket;
                bucketBegin = it;
            }
        };

        for (auto epsilon : {63, 31, 15, 7, 3}) {
            pgm::PolymorphicPGMIndex p(begin, end, epsilon);

            cost = p.size_in_bytes() * 8;
            if (cost >= bestCost) {
                // If PGM alone already is larger, additional epsilon parameters will be larger as well.
                break;
            }
            previousBucket = 0;
            bucketBegin = begin;
            p.for_each(begin, end, updateCost);
            updateCost(end, std::numeric_limits<uint64_t>::max());

            if (cost < bestCost) {
                pgm = std::move(p);
                bestCost = cost;
            }
        }
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const  {
        return pgm.approximate_rank(key);
    }

    template<typename Iterator, typename Func>
    void bucketOf(Iterator first, Iterator last, Func f) const {
        pgm.for_each(first, last, f);
    }

    [[nodiscard]] size_t space() const {
        return pgm.size_in_bytes();
    }

    [[nodiscard]] size_t segments_count() const {
        return pgm.segments_count();
    }

    [[nodiscard]] size_t numBuckets() const {
        return pgm.size();
    }

    static std::string name() {
        return "PolymorphicPGMBucketMapper";
    }

};