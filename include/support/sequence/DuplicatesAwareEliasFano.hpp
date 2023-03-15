#pragma once

#include <compact_elias_fano.hpp>
#include <partitioned_sequence.hpp>
#include <strict_sequence.hpp>
#include <vector>

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
