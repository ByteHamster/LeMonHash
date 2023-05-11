#pragma once

#include <cstdint>
#include <bit>
#include "support/PGM.hpp"

namespace lemonhash {
/**
 * Uses the PGM Index to get an approximate rank, which we use as bucket index.
 * Has small space overhead for uniform distribution, but enables using other distributions.
 */
struct PGMBucketMapper {
    pgm::PGMIndex<uint64_t> pgmIndex;
    size_t numBuckets_;

    template<typename RandomIt>
    PGMBucketMapper(RandomIt begin, RandomIt end) {
        auto best_space = std::numeric_limits<size_t>::max();

        for (auto epsilon : {63, 31, 15, 7, 3}) {
            pgm::PGMIndex<uint64_t> pgm(begin, end, epsilon);

            size_t space = pgm.size_in_bytes();
            if (space >= best_space) {
                // If PGM alone already is larger, additional epsilon parameters will be larger as well.
                break;
            }
            std::vector<size_t> bucket_sizes(std::distance(begin, end) + 1);
            pgm.for_each(begin, end, [&] (auto, auto approx_rank) {
                ++bucket_sizes[approx_rank];
            });

            size_t ranks_bits = 0;
            for (auto b: bucket_sizes)
                ranks_bits += b <= 1 ? 0ull : b * std::bit_width(b - 1);

            space += ranks_bits / 8;
            if (space < best_space) {
                pgmIndex = std::move(pgm);
                best_space = space;
            }
        }

        numBuckets_ = bucketOf(*std::prev(end)) + 1;
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const {
        return pgmIndex.approximate_rank(key);
    }

    template<typename Iterator, typename Func>
    void bucketOf(Iterator first, Iterator last, Func f) const {
        pgmIndex.for_each(first, last, f);
    }

    [[nodiscard]] size_t space() const {
        return sizeof(*this) + pgmIndex.size_in_bytes();
    }

    [[nodiscard]] size_t numBuckets() const {
        return numBuckets_;
    }

    static std::string name() {
        return "PGMBucketMapper";
    }

    [[nodiscard]] std::string info() const {
        return "epsilon=" + std::to_string(pgmIndex.epsilon_value());
    }
};
} // namespace lemonhash
