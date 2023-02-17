#pragma once

#include "BucketMapper.hpp"
#include "support/SuccinctPGM.hpp"

/**
 * Uses a succinct version of the PGM Index to get an approximate rank, which we use as bucket index.
 */
struct SuccinctPGMBucketMapper {
    union {
        pgm::SuccinctPGMIndex<> *pgmIndex = nullptr;
        uint64_t uOverN;
    };
    size_t numBuckets_ : 63;
    bool usesPgmIndex : 1;

    SuccinctPGMBucketMapper() = default;

    template<typename RandomIt>
    SuccinctPGMBucketMapper(RandomIt begin, RandomIt end) {
        auto n = (uint64_t) std::distance(begin, end);
        auto bestCost = std::numeric_limits<size_t>::max();

        size_t cost;
        size_t previousBucket;
        RandomIt bucketBegin;
        auto updateCost = [&](auto it, size_t bucket) {
            if (bucket != previousBucket) {
                auto bucketSize = (size_t) std::distance(bucketBegin, it);
                cost += bucketSize <= 1 ? 0 : bucketSize * std::bit_width(bucketSize - 1);
                previousBucket = bucket;
                bucketBegin = it;
            }
        };

        // Evaluate simple linear mapper
        auto tmpUOverN = *(end - 1) / (n - 1);
        cost = 0;
        previousBucket = 0;
        bucketBegin = begin;
        for (auto it = begin; it != end; ++it)
            updateCost(it, linearMapper(*it, n, tmpUOverN));
        updateCost(end, std::numeric_limits<uint64_t>::max());
        usesPgmIndex = false;
        uOverN = tmpUOverN;
        numBuckets_ = n;
        bestCost = cost;

        // Evaluate PGM
        for (auto epsilon : {3, 7, 15, 31, 63}) {
            auto *pgm = new pgm::SuccinctPGMIndex<>(begin, end, epsilon);

            cost = pgm->size_in_bytes() * 8;
            // If the cost of the PGM alone already is larger, we do not need to run for_each
            if (cost < bestCost) {
                previousBucket = 0;
                bucketBegin = begin;
                pgm->for_each(begin, end, updateCost);
                updateCost(end, std::numeric_limits<uint64_t>::max());
            }

            auto segmentsCount = pgm->segments_count();
            if (cost < bestCost) {
                if (usesPgmIndex) {
                    delete pgmIndex;
                }
                pgmIndex = pgm;
                usesPgmIndex = true;
                numBuckets_ = bucketOf(*std::prev(end)) + 1;
                bestCost = cost;
            } else {
                delete pgm;
            }
            if (segmentsCount == 1)
                break;
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
        return sizeof(*this) + (usesPgmIndex ? pgmIndex->size_in_bytes() : 0);
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

    ~SuccinctPGMBucketMapper() {
        if (usesPgmIndex)
            delete pgmIndex;
    }

private:

    [[nodiscard]] size_t linearMapper(uint64_t key, size_t numBuckets, uint64_t uOverN) const {
        return std::min(numBuckets - 1, key / uOverN);
    }
};
