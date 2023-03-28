#pragma once

#include "support/SuccinctPGM.hpp"

namespace lemonhash {
/**
 * Uses a succinct version of the PGM Index to get an approximate rank, which we use as bucket index.
 */
template <bool autoTune = true>
struct SuccinctPGMBucketMapper {
    pgm::SuccinctPGMIndex<> pgm;

    SuccinctPGMBucketMapper() = default;

    template<typename RandomIt>
    SuccinctPGMBucketMapper(RandomIt begin, RandomIt end) {
        if constexpr (!autoTune) {
            pgm = std::move(decltype(pgm)(begin, end, 31));
            return;
        }

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

        // Evaluate PGM
        for (auto epsilon : {63, 31, 15}) {
            decltype(pgm) p(begin, end, epsilon);

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

    [[nodiscard]] size_t numBuckets() const {
        return pgm.size();
    }

    static std::string name() {
        return std::string("SuccinctPGMBucketMapper") + (autoTune ? "AutoTune" : "FixedEps");
    }

    [[nodiscard]] std::string info() const {
        return "epsilon=" + std::to_string(pgm.epsilon_value());
    }
};
} // namespace lemonhash
