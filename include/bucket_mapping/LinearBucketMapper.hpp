#pragma once

/**
 * Each object is mapped linearly to its bucket. This only works well for uniform distributed inputs.
 */
template <float elsPerBucket>
struct LinearBucketMapper {
    size_t numBuckets;
    uint64_t u;

    template<typename RandomIt>
    LinearBucketMapper(RandomIt begin, RandomIt end)
            : numBuckets((end - begin) / elsPerBucket),
              u(*std::prev(end)) {
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const {
        return std::min(numBuckets - 1, key / (u / (numBuckets - 1)));
    }

    [[nodiscard]] size_t size() const {
        return sizeof(*this);
    }

    [[nodiscard]] constexpr static float elementsPerBucket() {
        return elsPerBucket;
    }

    static std::string name() {
        return std::string("LinearBucketMapper")
               + " elementsPerBucket=" + std::to_string(elementsPerBucket());
    }

    [[nodiscard]] std::string info() const {
        return "numBuckets=" + std::to_string(numBuckets) + " u=" + std::to_string(u);
    }
};
