#pragma once

class BitVectorBuilder {
    private:
        pasta::BitVector bitVector;
        size_t size = 0;
    public:
        BitVectorBuilder() {
            bitVector.resize(100);
        }

        bool append(bool x) {
            if (bitVector.size() <= size) {
                bitVector.resize(bitVector.size() * 2);
            }
            bitVector[size] = x;
            size++;
        }

        pasta::BitVector build() {
            bitVector.resize(size);
            return std::move(bitVector);
        }
};