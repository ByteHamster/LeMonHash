#pragma once

#include <cassert>
#include <utility>
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>
#include "support/util.h"

namespace util {

/** A randomly-accessible sequence of bucket offsets/sizes. */
class BucketOffsets {
private:
    pasta::BitVector bv; ///< Concatenates the unary representation of the bucket sizes.
    pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES> select;
    size_t position = 0;

public:

    /** Creates a new data structure for the given number of buckets and elements. */
    BucketOffsets(size_t buckets, size_t elements) : bv(buckets + elements + 1, false) {
        bv[0] = true;
    }

    /** Adds a new bucket with the given size. */
    void push(size_t bucketSize) {
        position += bucketSize + 1;
        assert(position < bv.size());
        bv[position] = true;
    }

    /** Finalizes the data structure. */
    void done() {
        select = {bv};
    }

    /** Returns a pair with the bucket offset and the bucket size. */
    [[nodiscard]] std::pair<uint64_t, uint64_t> at(size_t bucket) const {
        auto a = select.select1(bucket + 1);
        auto b = nextOne(a, bv.data().data());
        return {a - bucket, b - a - 1};
    }

    /** Returns the space usage of this data structure, in bytes. */
    [[nodiscard]] size_t space() const {
        return bv.size() / 8 + select.space_usage() + sizeof(*this);
    }
};

} // namespace util