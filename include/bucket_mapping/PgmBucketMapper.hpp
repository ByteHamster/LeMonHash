#pragma once

#include "support/PGM.hpp"

/**
 * Uses the PGM Index to get an approximate rank, which we use as bucket index.
 * Has small space overhead for uniform distribution, but enables using other distributions.
 */
struct PgmBucketMapper {
    pgm::PGMIndex<uint64_t> pgmIndex;
    size_t numBuckets;

    template<typename RandomIt>
    PgmBucketMapper(RandomIt begin, RandomIt end)
            : pgmIndex(begin, end, 31, 8), numBuckets(bucketOf(*std::prev(end)) + 1) {
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
};
