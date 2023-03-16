#pragma once

/**
 * Each object is mapped linearly to its bucket. This only works well for uniform distributed inputs.
 */
struct LinearBucketMapper {
    static constexpr float elsPerBucket = 1.0f;
    size_t numBuckets_;
    uint64_t u;

    template<typename RandomIt>
    LinearBucketMapper(RandomIt begin, RandomIt end)
            : numBuckets_((end - begin) / elsPerBucket),
              u(*std::prev(end)) {
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const {
        return std::min(numBuckets_ - 1, key / (u / (numBuckets_ - 1)));
    }

    template<typename Iterator, typename Func>
    void bucketOf(Iterator first, Iterator last, Func f) const {
        while (first != last) {
            f(first, bucketOf(*first));
            ++first;
        }
    }

    [[nodiscard]] size_t size() const {
        return sizeof(*this);
    }

    [[nodiscard]] size_t numBuckets() const {
        return numBuckets_;
    }

    [[nodiscard]] constexpr static float elementsPerBucket() {
        return elsPerBucket;
    }

    static std::string name() {
        return "LinearBucketMapper elementsPerBucket=" + std::to_string(elementsPerBucket());
    }

    [[nodiscard]] std::string info() const {
        return "numBuckets=" + std::to_string(numBuckets_) + " u=" + std::to_string(u);
    }
};
