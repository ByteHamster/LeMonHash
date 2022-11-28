#include <vector>
#include <algorithm>
#include <iostream>
#include <EliasFano.h>
#include <Function.h>
#include "MultiRetrievalDataStructure.hpp"
#include "bucket_mapping/LinearBucketMapper.hpp"
#include "bucket_mapping/PgmBucketMapper.hpp"
#include "bucket_mapping/SuccinctPgmBucketMapper.hpp"

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

            for (size_t i = 0; i < data.size(); i++) {
                // In a final version, this would be solved by a linear scan
                // instead of actually hashing objects to bucket data structures
                uint64_t &key = data.at(i);
                buckets.at(bucketMapper.bucketOf(key)).push_back(key);
                assert(i == 0 || data.at(i) > data.at(i - 1));
                assert(i == 0 || bucketMapper.bucketOf(data.at(i)) >= bucketMapper.bucketOf(data.at(i - 1)));
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
            std::cout<<"BucketMapper: "<<(8.0*bucketMapper.size()/N)
                     <<" ("<<bucketMapper.info()<<")"<<std::endl;
            size_t bytes = bucketMapper.size() + sizeof(*this) + bucketSizePrefix.space();
            return 8 * bytes + retrievalDataStructure.spaceBits(N);
        }
};
