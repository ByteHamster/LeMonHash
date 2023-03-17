#pragma once

#include <compact_elias_fano.hpp>
#include <partitioned_sequence.hpp>
#include <strict_sequence.hpp>
#include <vector>

namespace lemonhash {
/**
 * Stores a monotonic list of integers using Elias-Fano coding.
 * Uses a partitioned Elias-Fano implementation that optimizes irregularities.
 */
class PartitionedEliasFano {
    private:
        std::vector<size_t> inputData;
        succinct::bit_vector bv;
        ds2i::partitioned_sequence<ds2i::compact_elias_fano>::enumerator enumerator;
        size_t n = 0;
    public:
        PartitionedEliasFano() = default;

        PartitionedEliasFano(PartitionedEliasFano &&other) noexcept {
            inputData = std::move(other.inputData);
            other.bv.swap(bv);
            enumerator = other.enumerator;
            n = other.n;
        }

        void append(size_t x) {
            assert(bv.size() == 0); // Not already completed
            n++;
            inputData.push_back(x);
        }

        size_t size() const {
            return n;
        }

        void complete() {
            if (n == 0) {
                return;
            }
            size_t universe = inputData.back() + 1;

            succinct::bit_vector_builder bvb;
            ds2i::global_parameters params;
            ds2i::partitioned_sequence<ds2i::compact_elias_fano>::write(
                    bvb, inputData.begin(), universe, inputData.size(), params);
            succinct::bit_vector(&bvb).swap(bv);
            enumerator = ds2i::partitioned_sequence<ds2i::compact_elias_fano>::enumerator(
                    bv, 0, universe, inputData.size(), params);

            inputData.clear();
            inputData.shrink_to_fit();
        }

        size_t at(size_t i) {
            assert(i < n);
            return enumerator.move(i).second;
        }

        size_t bit_size() {
            return bv.size() + 8 * sizeof(*this);
        }
};
} // namespace lemonhash
