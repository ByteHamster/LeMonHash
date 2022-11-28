#pragma once

#include "support/PGM.hpp"
#include "support/util.hpp"

/**
 * Uses the PGM Index to get an approximate rank, which we use as bucket index.
 * Has small space overhead for uniform distribution, but enables using other distributions.
 */
struct PgmBucketMapper {
    pgm::PGMIndex<uint64_t> pgmIndex;
    size_t numBuckets;

    template<typename RandomIt>
    PgmBucketMapper(RandomIt begin, RandomIt end) {
        auto best_space = std::numeric_limits<size_t>::max();

        for (auto epsilon : {3, 7, 15, 31, 63}) {
            pgm::PGMIndex<uint64_t> pgm(begin, end, epsilon);

            std::vector<size_t> bucket_sizes(std::distance(begin, end) + 1);
            pgm.for_each(begin, end, [&] (auto, auto approx_rank) {
                ++bucket_sizes[approx_rank];
            });

            size_t ranks_bits = 0;
            for (auto b: bucket_sizes)
                ranks_bits += b <= 1 ? 0ull : b * BIT_WIDTH(b - 1);

            auto space = pgm.size_in_bytes() + ranks_bits / 8;
            if (space < best_space) {
                pgmIndex = pgm;
                best_space = space;
            }
        }

        numBuckets = bucketOf(*std::prev(end)) + 1;
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const {
        return pgmIndex.approximate_rank(key);
    }

    [[nodiscard]] size_t size() const {
        return sizeof(*this) + pgmIndex.size_in_bytes();
    }

    [[nodiscard]] constexpr static float elementsPerBucket() {
        return 1.0;
    }

    static std::string name() {
        return std::string("PgmBucketMapper");
    }

    [[nodiscard]] std::string info() const {
        return "epsilon=" + std::to_string(pgmIndex.epsilon_value());
    }
};
