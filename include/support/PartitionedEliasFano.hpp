#pragma once

#include <compact_elias_fano.hpp>
#include <partitioned_sequence.hpp>
#include <strict_sequence.hpp>
#include <vector>

/**
 * Stores a monotonic list of integers using Elias-Fano coding.
 * Uses a partitioned Elias-Fano implementation that optimizes irregularities.
 * Additionally, filters out runs of repeated same values using rank/select.
 */
class PartitionedEliasFano {
    private:
        std::vector<size_t> inputData;
        succinct::bit_vector bv;
        ds2i::partitioned_sequence<ds2i::compact_elias_fano>::enumerator enumerator;
        size_t n = 0;
        static constexpr size_t DUPLICATE_FILTER_GRANULARITY = 3;
        std::vector<size_t> insertionBuffer;
        pasta::BitVector duplicateFilterIsAreaStored;
        pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES> duplicateFilterRank;
    public:
        PartitionedEliasFano() {
            duplicateFilterIsAreaStored.resize(100, false);
            duplicateFilterIsAreaStored[0] = true;
        }

        PartitionedEliasFano(PartitionedEliasFano &&other) noexcept {
            inputData = std::move(other.inputData);
            other.bv.swap(bv);
            enumerator = other.enumerator;
            n = other.n;
            insertionBuffer = std::move(other.insertionBuffer);
            duplicateFilterIsAreaStored = std::move(other.duplicateFilterIsAreaStored);
            duplicateFilterRank = std::move(other.duplicateFilterRank);
        }

        void append(size_t x) {
            assert(bv.size() == 0); // Not already completed
            n++;
            if (n <= DUPLICATE_FILTER_GRANULARITY) {
                inputData.push_back(x);
            } else {
                insertionBuffer.push_back(x);
                flushInsertionBuffer(false);
            }
        }

        void flushInsertionBuffer(bool forceFlush) {
            if (insertionBuffer.size() < DUPLICATE_FILTER_GRANULARITY && !forceFlush) {
                return;
            }
            if (insertionBuffer.empty()) {
                return;
            }
            bool allSame = true;
            for (size_t i : insertionBuffer) {
                allSame = allSame && i == inputData.back();
            }
            size_t area = (n - 1) / DUPLICATE_FILTER_GRANULARITY;
            if (!allSame || forceFlush) {
                for (size_t i : insertionBuffer) {
                    inputData.push_back(i);
                }
                if (duplicateFilterIsAreaStored.size() < area + 1) {
                    duplicateFilterIsAreaStored.resize(area * 2, false);
                }
                duplicateFilterIsAreaStored[area] = true;
            }
            insertionBuffer.clear();
        }

        size_t size() const {
            return n;
        }

        void complete() {
            if (n == 0) {
                return;
            }
            flushInsertionBuffer(true);
            size_t universe = inputData.back() + 1;

            succinct::bit_vector_builder bvb;
            ds2i::global_parameters params;
            ds2i::partitioned_sequence<ds2i::compact_elias_fano>::write(
                    bvb, inputData.begin(), universe, inputData.size(), params);
            succinct::bit_vector(&bvb).swap(bv);
            enumerator = ds2i::partitioned_sequence<ds2i::compact_elias_fano>::enumerator(
                    bv, 0, universe, inputData.size(), params);

            duplicateFilterIsAreaStored.resize(n / DUPLICATE_FILTER_GRANULARITY + 1, false);
            duplicateFilterRank = pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES>(duplicateFilterIsAreaStored);

            inputData.clear();
            inputData.shrink_to_fit();
        }

        size_t at(size_t i) {
            assert(i < n);
            size_t area = i / DUPLICATE_FILTER_GRANULARITY;
            size_t areaLocation = DUPLICATE_FILTER_GRANULARITY * duplicateFilterRank.rank1(area);
            size_t index = areaLocation;
            bool isAreaStored = duplicateFilterIsAreaStored[area];
            if (isAreaStored) {
                index += i % DUPLICATE_FILTER_GRANULARITY;
            } else {
                // Take last key of previous area
                index--;
            }
            assert(index <= i);
            assert(index < enumerator.size());
            return enumerator.move(index).second;
        }

        size_t bit_size() {
            return bv.size() + 8 * sizeof(*this) + duplicateFilterIsAreaStored.size() + 8 * duplicateFilterRank.space_usage();
        }
};

class DuplicatesAwareEliasFano {
    struct InputData {
        std::vector<size_t> distinct;
        std::vector<size_t> duplicates;
    };

    using pefi = ds2i::partitioned_sequence<ds2i::indexed_sequence>;
    using pefs = ds2i::partitioned_sequence<ds2i::strict_sequence>;

    succinct::bit_vector bv;
    pefs::enumerator distinctEnumerator;   ///< the values of the unique elements
    pefi::enumerator duplicatesEnumerator; ///< the indexes of the repeated (or unique) elements (if complement=true)
    InputData *input = nullptr;
    size_t n : 30 = 0;
    bool hasDuplicates : 1;
    bool complement : 1;

public:

    DuplicatesAwareEliasFano() {
        input = new InputData();
        input->distinct.reserve(128);
        input->duplicates.reserve(128);
    }

    DuplicatesAwareEliasFano(DuplicatesAwareEliasFano &&other) noexcept {
        other.bv.swap(bv);
        distinctEnumerator = other.distinctEnumerator;
        duplicatesEnumerator = other.duplicatesEnumerator;
        input = other.input;
        other.input = nullptr;
        n = other.n;
        hasDuplicates = other.hasDuplicates;
    }

    void complete() {
        if (n == 0)
            return;

        hasDuplicates = !input->duplicates.empty();

        succinct::bit_vector_builder bvb;
        ds2i::global_parameters params;
        pefs::write(bvb, input->distinct.begin(), input->distinct.back() + 1, input->distinct.size(), params);
        auto offset = bvb.size();
        if (hasDuplicates) {
            if (n - input->duplicates.size() < input->duplicates.size()) {
                complement = true;
                std::vector<size_t> complemented = {0};
                complemented.reserve(n - input->duplicates.size());
                uint64_t prev = 0;
                for (auto x : input->duplicates) {
                    for (auto i = prev + 1; i < x; ++i)
                        complemented.push_back(i);
                    prev = x;
                }
                for (size_t i = prev + 1; i < n; ++i)
                    complemented.push_back(i);
                assert(complemented.size() == n - input->duplicates.size());
                input->duplicates = std::move(complemented);
            } else {
                complement = false;
            }
            pefi::write(bvb, input->duplicates.begin(), input->duplicates.back() + 1, input->duplicates.size(), params);
        }
        succinct::bit_vector(&bvb).swap(bv);
        distinctEnumerator = {bv, 0, input->distinct.back() + 1, input->distinct.size(), params};

        if (hasDuplicates)
            duplicatesEnumerator = {bv, offset, input->duplicates.back() + 1, input->duplicates.size(), params};

        delete input;
        input = nullptr;
    }

    void append(size_t x) {
        if (n == 0 || input->distinct.back() != x)
            input->distinct.push_back(x);
        else if (input->distinct.back() == x)
            input->duplicates.push_back(n);
        ++n;
    }

    size_t size() const {
        return n;
    }

    size_t at(size_t i) {
        assert(i < n);
        size_t j = i;
        if (hasDuplicates) {
            auto [pos, val] = duplicatesEnumerator.next_geq(i);
            if (complement)
                j = pos - ((val > i && pos > 0));
            else
                j = i - pos - (val == i && pos < duplicatesEnumerator.size());
        }
        return distinctEnumerator.move(j).second;
    }

    size_t bit_size() {
        return bv.size() + 8 * sizeof(*this);
    }

    ~DuplicatesAwareEliasFano() {
        delete input;
    }
};