#pragma once

#include "support/SuccinctPGM.hpp"
#include "support/util.hpp"

/**
 * Uses a succinct version of the PGM Index to get an approximate rank, which we use as bucket index.
 */
struct SuccinctPgmBucketMapper {
    pgm::SuccinctPGMIndex<uint64_t> *pgmIndex = nullptr;
    size_t numBuckets;

    template<typename RandomIt>
    SuccinctPgmBucketMapper(RandomIt begin, RandomIt end) {
        auto best_space = std::numeric_limits<size_t>::max();

        for (auto epsilon : {3, 7, 15, 31, 63}) {
            auto *pgm = new pgm::SuccinctPGMIndex<uint64_t>(begin, end, epsilon);

            std::vector<size_t> bucket_sizes(std::distance(begin, end) + 1);
            pgm->for_each(begin, end, [&] (auto, auto approx_rank) {
                ++bucket_sizes[approx_rank];
            });

            size_t ranks_bits = 0;
            for (auto b: bucket_sizes)
                ranks_bits += b <= 1 ? 0ull : b * BIT_WIDTH(b - 1);

            auto space = pgm->size_in_bytes() + ranks_bits / 8;
            if (space < best_space) {
                delete pgmIndex;
                pgmIndex = pgm;
                best_space = space;
            } else {
                delete pgm;
            }
        }

        numBuckets = bucketOf(*std::prev(end)) + 1;
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const {
        return pgmIndex->approximate_rank(key);
    }

    [[nodiscard]] size_t size() const {
        return sizeof(*this) + sizeof(*pgmIndex) + pgmIndex->size_in_bytes();
    }

    [[nodiscard]] constexpr static float elementsPerBucket() {
        return 1.0;
    }

    static std::string name() {
        return std::string("SuccinctPgmBucketMapper");
    }

    [[nodiscard]] std::string info() const {
        return "epsilon=" + std::to_string(pgmIndex->epsilon_value());
    }

    ~SuccinctPgmBucketMapper() { delete pgmIndex; }
};
