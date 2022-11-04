#include <vector>
#include <algorithm>
#include <iostream>
#include <EliasFano.h>
#include <SimpleRibbon.h>
#include <Function.h>
#include <pgm/pgm_index.hpp>

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
    private:
        size_t N;
        BucketMapper bucketMapper;
        SimpleRibbon<1> *retrieval1;
        SimpleRibbon<2> *retrieval2;
        SimpleRibbon<3> *retrieval3;
        SimpleRibbon<4> *retrieval4;
        SimpleRibbon<5> *retrieval5;
        SimpleRibbon<6> *retrieval6;
        SimpleRibbon<7> *retrieval7;
        SimpleRibbon<32, 64, uint32_t> *retrieval8Plus;
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
            std::vector<std::pair<uint64_t, uint8_t>> retrieval1Input;
            std::vector<std::pair<uint64_t, uint8_t>> retrieval2Input;
            std::vector<std::pair<uint64_t, uint8_t>> retrieval3Input;
            std::vector<std::pair<uint64_t, uint8_t>> retrieval4Input;
            std::vector<std::pair<uint64_t, uint8_t>> retrieval5Input;
            std::vector<std::pair<uint64_t, uint8_t>> retrieval6Input;
            std::vector<std::pair<uint64_t, uint8_t>> retrieval7Input;
            std::vector<std::pair<uint64_t, uint32_t>> retrieval8PlusInput;

            buckets.resize(bucketMapper.numBuckets);
            retrieval1Input.reserve(0.3 * N);
            retrieval2Input.reserve(0.24 * N);
            retrieval3Input.reserve(0.01 * N);

            for (uint64_t key : data) {
                buckets.at(bucketMapper.bucketOf(key)).push_back(key);
            }
            size_t bucketSizePrefixTemp = 0;
            for (auto & bucket : buckets) {
                size_t bucketSize = bucket.size();
                for (size_t j = 0; j < bucketSize; j++) {
                    if (bucketSize <= 1) {
                        // Done
                    } else if (bucketSize <= 2) {
                        retrieval1Input.emplace_back(std::make_pair(bucket.at(j), j));
                    } else if (bucketSize <= 4) {
                        retrieval2Input.emplace_back(std::make_pair(bucket.at(j), j));
                    } else if (bucketSize <= 8) {
                        retrieval3Input.emplace_back(std::make_pair(bucket.at(j), j));
                    } else if (bucketSize <= 16) {
                        retrieval4Input.emplace_back(std::make_pair(bucket.at(j), j));
                    } else if (bucketSize <= 32) {
                        retrieval5Input.emplace_back(std::make_pair(bucket.at(j), j));
                    } else if (bucketSize <= 64) {
                        retrieval6Input.emplace_back(std::make_pair(bucket.at(j), j));
                    } else if (bucketSize <= 128) {
                        retrieval7Input.emplace_back(std::make_pair(bucket.at(j), j));
                    } else {
                        retrieval8PlusInput.emplace_back(std::make_pair(bucket.at(j), j));
                    }
                }
                bucketSizePrefix.push_back(bucketSizePrefixTemp);
                bucketSizePrefixTemp += bucketSize;
            }
            bucketSizePrefix.push_back(bucketSizePrefixTemp);
            bucketSizePrefix.buildRankSelect();

            retrieval1 = new SimpleRibbon<1>(retrieval1Input);
            retrieval2 = new SimpleRibbon<2>(retrieval2Input);
            retrieval3 = new SimpleRibbon<3>(retrieval3Input);
            retrieval4 = new SimpleRibbon<4>(retrieval4Input);
            retrieval5 = new SimpleRibbon<5>(retrieval5Input);
            retrieval6 = new SimpleRibbon<6>(retrieval6Input);
            retrieval7 = new SimpleRibbon<7>(retrieval7Input);
            retrieval8Plus = new SimpleRibbon<32, 64, uint32_t>(retrieval8PlusInput);
        }

        ~DirectRankStoringMmphf() {
            delete retrieval1;
            delete retrieval2;
            delete retrieval3;
            delete retrieval4;
            delete retrieval5;
            delete retrieval6;
            delete retrieval7;
            delete retrieval8Plus;
        }

        size_t operator()(uint64_t key) {
            auto ptr = bucketSizePrefix.at(bucketMapper.bucketOf(key));
            size_t bucketOffset = *ptr;
            ++ptr;
            size_t nextBucketOffset = *ptr;
            size_t bucketSize = nextBucketOffset - bucketOffset;
            if (bucketSize <= 1) {
                return bucketOffset;
            } else if (bucketSize <= 2) {
                return bucketOffset + retrieval1->retrieve(key);
            } else if (bucketSize <= 4) {
                return bucketOffset + retrieval2->retrieve(key);
            } else if (bucketSize <= 8) {
                return bucketOffset + retrieval3->retrieve(key);
            } else if (bucketSize <= 16) {
                return bucketOffset + retrieval4->retrieve(key);
            } else if (bucketSize <= 32) {
                return bucketOffset + retrieval5->retrieve(key);
            } else if (bucketSize <= 64) {
                return bucketOffset + retrieval6->retrieve(key);
            } else if (bucketSize <= 128) {
                return bucketOffset + retrieval7->retrieve(key);
            } else {
                return bucketOffset + retrieval8Plus->retrieve(key);
            }
        }

        size_t spaceBits() {
            std::cout<<"Space usage (bits/key) of the retrieval data structures:"<<std::endl;
            std::cout<<"1-bit: "<<(8.0*retrieval1->size()/N)<<", ";
            std::cout<<"2-bit: "<<(8.0*retrieval2->size()/N)<<", ";
            std::cout<<"3-bit: "<<(8.0*retrieval3->size()/N)<<", ";
            std::cout<<"4-bit: "<<(8.0*retrieval4->size()/N)<<std::endl;
            std::cout<<"5-bit: "<<(8.0*retrieval5->size()/N)<<", ";
            std::cout<<"6-bit: "<<(8.0*retrieval6->size()/N)<<", ";
            std::cout<<"7-bit: "<<(8.0*retrieval7->size()/N)<<", ";
            std::cout<<"8+-bit: "<<(8.0*retrieval8Plus->size()/N)<<std::endl;
            std::cout<<"EliasFano:    "<<(8.0*bucketSizePrefix.space()/N)<<std::endl;
            std::cout<<"BucketMapper: "<<(8.0*bucketMapper.size()/N)<<std::endl;
            size_t bytes = retrieval1->size() + retrieval2->size() + retrieval3->size()
                           + retrieval4->size() + retrieval5->size()
                           + retrieval6->size() + retrieval7->size()
                           + retrieval8Plus->size()
                           + bucketMapper.size();
            bytes += sizeof(*this);
            bytes += bucketSizePrefix.space();
            return 8 * bytes;
        }
};
