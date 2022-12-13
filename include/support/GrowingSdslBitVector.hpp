#pragma once

#include <pasta/bit_vector/bit_vector.hpp>

class GrowingSdslBitVector {
    public:
        size_t size = 0;
        sdsl::bit_vector data;

        void append(bool x) {
            if (data.size() < size + 1) {
                data.resize(2 * data.size() + 1);
            }
            data[size] = x;
            size++;
        }

        void shrinkToFit() {
            data.resize(size);
        }
};