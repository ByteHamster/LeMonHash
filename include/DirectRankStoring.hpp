#include <vector>
#include <algorithm>
#include <iostream>
#include <EliasFano.h>
#include <Function.h>
#include <pgm/pgm_index.hpp>
#include "MultiRetrievalDataStructure.hpp"

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
};

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

/**
 * Monotone Minimal Perfect Hash Function (MMPHF) using the Direct Rank Storing (DRS) technique.
 * Each object is mapped to a bucket (using different techniques).
 * Within the buckets, a retrieval data structure explicitly stores the ranks of all objects.
 * Given that we expect only one object per bucket, we can usually use a 1-bit retrieval data structure.
 * The bucket sizes are stored with Elias-Fano.
 */
template <typename BucketMapper>
class DirectRankStoringMmphf {
    public:
        size_t N;
    private:
        BucketMapper bucketMapper;
        MultiRetrievalDataStructure retrievalDataStructure;
        util::EliasFano<util::floorlog2(std::max(1.0f, BucketMapper::elementsPerBucket()))> bucketSizePrefix;
    public:
        static std::string name() {
            return "DirectRankStoringMmphf bucketMapper=" + BucketMapper::name();
        }

        explicit DirectRankStoringMmphf(std::vector<uint64_t> &data)
                : N(data.size()),
                  bucketMapper(data.begin(), data.end()),
                  bucketSizePrefix(bucketMapper.numBuckets + 1, data.size() + 1) {
            std::vector<std::vector<uint64_t>> buckets;
            buckets.resize(bucketMapper.numBuckets);

            for (uint64_t key : data) {
                // In a final version, this would be solved by a linear scan
                // instead of actually hashing objects to bucket data structures
                buckets.at(bucketMapper.bucketOf(key)).push_back(key);
            }
            size_t bucketSizePrefixTemp = 0;
            for (auto & bucket : buckets) {
                size_t bucketSize = bucket.size();
                for (size_t j = 0; j < bucketSize; j++) {
                    if (bucketSize > 1) {
                        retrievalDataStructure.addInput(bucketSize, bucket.at(j), j);
                    }
                }
                bucketSizePrefix.push_back(bucketSizePrefixTemp);
                bucketSizePrefixTemp += bucketSize;
            }
            bucketSizePrefix.push_back(bucketSizePrefixTemp);
            bucketSizePrefix.buildRankSelect();

            retrievalDataStructure.build();
        }

        size_t operator()(uint64_t key) {
            auto ptr = bucketSizePrefix.at(bucketMapper.bucketOf(key));
            size_t bucketOffset = *ptr;
            ++ptr;
            size_t nextBucketOffset = *ptr;
            size_t bucketSize = nextBucketOffset - bucketOffset;
            if (bucketSize <= 1) {
                return bucketOffset;
            } else {
                return bucketOffset + retrievalDataStructure.query(bucketSize, key);
            }
        }

        size_t spaceBits() {
            std::cout<<"EliasFano:    "<<(8.0*bucketSizePrefix.space()/N)<<std::endl;
            std::cout<<"BucketMapper: "<<(8.0*bucketMapper.size()/N)<<std::endl;
            size_t bytes = bucketMapper.size() + sizeof(*this) + bucketSizePrefix.space();
            return 8 * bytes + retrievalDataStructure.spaceBits(N);
        }
};
