#pragma once

#include "BucketMapper.hpp"

/**
 * Each object is mapped linearly to its bucket. This only works well for uniform distributed inputs.
 */
struct LinearBucketMapper : public BucketMapper {
    static constexpr float elsPerBucket = 1.0f;
    size_t numBuckets_;
    uint64_t u;

    template<typename RandomIt>
    LinearBucketMapper(RandomIt begin, RandomIt end)
            : numBuckets_((end - begin) / elsPerBucket),
              u(*std::prev(end)) {
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const final {
        return std::min(numBuckets_ - 1, key / (u / (numBuckets_ - 1)));
    }

    void bucketOf(Iterator first, Iterator last, Func f) const final {
        while (first != last) {
            f(first, bucketOf(*first));
            ++first;
        }
    }

    [[nodiscard]] size_t size() const final {
        return sizeof(*this);
    }

    [[nodiscard]] size_t numBuckets() const final {
        return numBuckets_;
    }

    [[nodiscard]] constexpr static float elementsPerBucket() {
        return elsPerBucket;
    }

    static std::string name() {
        return std::string("LinearBucketMapper")
               + " elementsPerBucket=" + std::to_string(elementsPerBucket());
    }

    [[nodiscard]] std::string info() const {
        return "numBuckets=" + std::to_string(numBuckets_) + " u=" + std::to_string(u);
    }

    ~LinearBucketMapper() override = default;
};
