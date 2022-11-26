#pragma once

#include "support/SuccinctPGM.hpp"

/**
 * Uses a succinct version of the PGM Index to get an approximate rank, which we use as bucket index.
 */
template <size_t Epsilon = 31>
struct SuccinctPgmBucketMapper {
    pgm::SuccinctPGMIndex<uint64_t, Epsilon> pgmIndex;
    size_t numBuckets;

    template<typename RandomIt>
    SuccinctPgmBucketMapper(RandomIt begin, RandomIt end)
            : pgmIndex(begin, end), numBuckets(bucketOf(*std::prev(end)) + 1) {
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
        return std::string("SuccinctPgmBucketMapper")
                + " epsilon=" + std::to_string(Epsilon);
    }
};
