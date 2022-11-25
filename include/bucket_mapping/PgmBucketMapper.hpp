#pragma once

#include <pgm/pgm_index.hpp>

/**
 * Uses the PGM Index to get an approximate rank, which we use as bucket index.
 * Has small space overhead for uniform distribution, but enables using other distributions.
 */
template <float elsPerBucket, size_t Epsilon=31>
struct PgmBucketMapper {
    pgm::PGMIndex<uint64_t, Epsilon, 8> pgmIndex;
    size_t numBuckets;

    template<typename RandomIt>
    PgmBucketMapper(RandomIt begin, RandomIt end)
            : pgmIndex(begin, end), numBuckets(bucketOf(*std::prev(end)) + 1) {
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const {
        return std::floor((float) pgmIndex.search(key).pos / elsPerBucket);
    }

    [[nodiscard]] size_t size() const {
        return sizeof(*this) + pgmIndex.size_in_bytes();
    }

    [[nodiscard]] constexpr static float elementsPerBucket() {
        return elsPerBucket;
    }

    static std::string name() {
        return std::string("PgmBucketMapper")
               + " elementsPerBucket=" + std::to_string(elementsPerBucket())
               + " epsilon=" + std::to_string(Epsilon);
    }
};
