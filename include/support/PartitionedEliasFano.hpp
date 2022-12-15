#pragma once

#include <compact_elias_fano.hpp>
#include <partitioned_sequence.hpp>
#include <vector>

class PartitionedEliasFano {
    private:
        succinct::bit_vector bv;
        ds2i::partitioned_sequence<ds2i::compact_elias_fano>::enumerator enumerator;
    public:
        explicit PartitionedEliasFano(std::vector<size_t> &vector) {
            size_t n = vector.size();
            size_t universe = vector.back() + 1;

            succinct::bit_vector_builder bvb;
            ds2i::global_parameters params;
            ds2i::partitioned_sequence<ds2i::compact_elias_fano>::write(bvb, vector.begin(), universe, n, params);
            succinct::bit_vector(&bvb).swap(bv);
            enumerator = ds2i::partitioned_sequence<ds2i::compact_elias_fano>::enumerator(bv, 0, universe, n, params);
        }

        size_t at(size_t i) {
            return enumerator.move(i).second;
        }

        size_t bit_size() {
            return bv.size() + 8 * sizeof(*this);
        }
};
