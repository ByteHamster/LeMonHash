#pragma once

#include "support/ShortPGM.hpp"
#include "support/util.hpp"

#pragma pack(push, 1)

/**
 * Uses a succinct version of the PGM Index to get an approximate rank, which we use as bucket index.
 */
struct ShortPgmBucketMapper {
    static constexpr uint8_t typesBitWidth = BIT_WIDTH(std::tuple_size_v<pgm::ShortPGMIndex<>::all_types> - 1);
    void *ptr = nullptr;
    size_t numBuckets_ : 48 - typesBitWidth;
    size_t type : typesBitWidth;

    ShortPgmBucketMapper() = default;

    template<typename RandomIt>
    ShortPgmBucketMapper(RandomIt begin, RandomIt end) {
        auto n = (uint64_t) std::distance(begin, end);
        auto bestCost = std::numeric_limits<size_t>::max();

        size_t cost;
        size_t previousBucket;
        RandomIt bucketBegin;
        auto updateCost = [&](auto it, size_t bucket) {
            if (bucket != previousBucket) {
                auto bucketSize = std::distance(bucketBegin, it);
                cost += bucketSize <= 1 ? 0 : bucketSize * BIT_WIDTH(bucketSize - 1);
                previousBucket = bucket;
                bucketBegin = it;
            }
        };

        for (auto epsilon : {3, 7, 15, 31, 63}) {
            auto [p, t] = pgm::ShortPGMIndex<>::make(begin, end, epsilon);

            pgm::ShortPGMIndex<>::visit([&](auto pgm) { cost = pgm->size_in_bytes() * 8; }, p, t);
            previousBucket = 0;
            bucketBegin = begin;
            pgm::ShortPGMIndex<>::visit([&](auto pgm) { pgm->for_each(begin, end, updateCost, n); }, p, t);
            updateCost(end, 0);

            if (cost < bestCost) {
                if (ptr)
                    visit([&](auto pgm) { delete pgm; });
                ptr = p;
                type = t;
                bestCost = cost;
            } else {
                pgm::ShortPGMIndex<>::visit([&](auto pgm) { delete pgm; }, p, t);
            }
        }

        numBuckets_ = n;
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const  {
        size_t result;
        visit([&](auto pgm) { result = pgm->approximate_rank(key, numBuckets_); });
        return result;
    }

    template<typename Iterator, typename Func>
    void bucketOf(Iterator first, Iterator last, Func f) const {
        visit([&](auto pgm) {  pgm->for_each(first, last, f, numBuckets_); });
    }

    [[nodiscard]] size_t size() const {
        size_t bytes;
        visit([&](auto pgm) { bytes = pgm->size_in_bytes(); });
        return sizeof(*this) + bytes;
    }

    [[nodiscard]] size_t numBuckets() const {
        return numBuckets_;
    }

    static std::string name() {
        return std::string("ShortPgmBucketMapper");
    }

    ~ShortPgmBucketMapper() {
        visit([](auto pgm) { delete pgm; });
    }

private:

    template <typename F>
    void visit(F f) {
        pgm::ShortPGMIndex<>::visit(f, ptr, type);
    }

    template <typename F>
    void visit(F f) const {
        pgm::ShortPGMIndex<>::visit(f, ptr, type);
    }

};

#pragma pack(pop)