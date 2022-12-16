#pragma once

#include "BucketMapper.hpp"
#include "support/SuccinctPGM.hpp"
#include "support/util.hpp"

/**
 * Uses a succinct version of the PGM Index to get an approximate rank, which we use as bucket index.
 */
struct SuccinctPgmBucketMapper {
    union {
        pgm::SuccinctPGMIndex<uint64_t> *pgmIndex = nullptr;
        uint64_t uOverN;
    };
    size_t numBuckets_ : 63;
    bool usesPgmIndex : 1;

    SuccinctPgmBucketMapper() = default;

    template<typename RandomIt>
    SuccinctPgmBucketMapper(RandomIt begin, RandomIt end) {
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
            auto *pgm = new pgm::SuccinctPGMIndex<uint64_t>(begin, end, epsilon);

            cost = pgm->size_in_bytes() * 8;
            previousBucket = 0;
            bucketBegin = begin;
            pgm->for_each(begin, end, updateCost);
            updateCost(end, 0);

            if (cost < bestCost) {
                delete pgmIndex;
                pgmIndex = pgm;
                bestCost = cost;
            } else {
                delete pgm;
            }
        }

        // Evaluate if linear mapper is better
        auto tmpUOverN = *(end - 1) / (n - 1);
        cost = 0;
        previousBucket = 0;
        bucketBegin = begin;
        for (auto it = begin; it != end; ++it)
            updateCost(it, linearMapper(*it, n, tmpUOverN));
        updateCost(end, 0);

        if (cost < bestCost) {
            delete pgmIndex;
            usesPgmIndex = false;
            uOverN = tmpUOverN;
            numBuckets_ = n;
        } else {
            usesPgmIndex = true;
            numBuckets_ = bucketOf(*std::prev(end)) + 1;
        }
    }

    [[nodiscard]] size_t bucketOf(uint64_t key) const  {
        if (usesPgmIndex)
            return pgmIndex->approximate_rank(key);
        return linearMapper(key, numBuckets_, uOverN);
    }

    template<typename Iterator, typename Func>
    void bucketOf(Iterator first, Iterator last, Func f) const {
        if (usesPgmIndex) {
            pgmIndex->for_each(first, last, f);
        } else {
            while (first != last) {
                f(first, bucketOf(*first));
                ++first;
            }
        }
    }

    [[nodiscard]] size_t size() const {
        return sizeof(*this) + (usesPgmIndex ? sizeof(*pgmIndex) + pgmIndex->size_in_bytes() : 0);
    }

    [[nodiscard]] size_t numBuckets() const {
        return numBuckets_;
    }

    [[nodiscard]] constexpr static float elementsPerBucket() {
        return 1.0;
    }

    static std::string name() {
        return std::string("SuccinctPgmBucketMapper");
    }

    [[nodiscard]] std::string info() const {
        if (usesPgmIndex)
            return "epsilon=" + std::to_string(pgmIndex->epsilon_value());
        return "u/n=" + std::to_string(uOverN);
    }

    ~SuccinctPgmBucketMapper() {
        if (usesPgmIndex)
            delete pgmIndex;
    }

private:

    [[nodiscard]] size_t linearMapper(uint64_t key, size_t numBuckets, uint64_t uOverN) const {
        return std::min(numBuckets - 1, key / uOverN);
    }
};
