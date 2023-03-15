#pragma once

#include <vector>
#include <cassert>
#include <sdsl/bit_vectors.hpp>
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>

namespace util {

// Returns the width of the integers stored in the low part of an Elias-Fano-coded sequence.
//
// This is based on:
//
//   Ma, Puglisi, Raman, Zhukova:
//   On Elias-Fano for Rank Queries in FM-Indexes.
//   DCC 2021.
//
// Implementation credit: Jouni Siren, https://github.com/vgteam/sdsl-lite
static uint8_t elias_fano_lo_width(size_t universe, size_t ones) {
    uint8_t low_width = 1;
    // Multisets with too many ones will have width 1.
    if (ones > 0 && ones <= universe) {
        double ideal_width = std::log2((static_cast<double>(universe) * std::log(2.0)) / static_cast<double>(ones));
        low_width = (uint8_t) std::round(std::max(ideal_width, 1.0));
    }
    return low_width;
}

/**
 * Compressed, monotone integer array.
 * Commonly used to store the positions of 1-bits in sparse bit vectors.
 */
class EliasFanoM {
private:
    sdsl::int_vector<> L;
    pasta::BitVector H;
    size_t count = 0;
    size_t universeSize = 0;
    pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES> *rankSelect = nullptr;

#ifndef NDEBUG
    uint64_t previousInsert = 0;
#endif

    [[nodiscard]] uint8_t lowerBits() const { return L.width(); }

    [[nodiscard]] uint64_t maskLowerBits() const { return sdsl::bits::lo_set[lowerBits()]; }

public:

    /**
     * Efficient pointer into an Elias-Fano coded sequence.
     * When incrementing/decrementing and reading, no additional select query is performed.
     */
    struct ElementPointer {
    private:
        size_t positionL;
        size_t positionH;
        const EliasFanoM *fano;

    public:
        ElementPointer(size_t positionH, size_t positionL, const EliasFanoM &fano)
            : positionL(positionL), positionH(positionH), fano(&fano) {
            assert(fano.H[positionH] == 1);
        }

        ElementPointer &operator++() {
            if (positionL >= fano->count - 1) {
                // Incremented more than the number of elements in the sequence.
                // Dereferencing it now is undefined behavior but decrementing again makes it usable again.
                positionL++;
                return *this;
            }
            assert(fano->H[positionH] == 1);
            positionL++;
            positionH++;
            while (fano->H[positionH] == 0) {
                positionH++;
            }
            assert(fano->H[positionH] == 1);
            return *this;
        }

        ElementPointer &operator--() {
            if (positionL >= fano->count) {
                // Was incremented more than the number of elements in the sequence.
                // Will be dereferenceable again if decremented to be inside the bounds.
                positionL--;
                return *this;
            }
            assert(positionL > 0);
            assert(fano->H[positionH] == 1);
            positionL--;
            positionH--;
            while (positionH > 0 && fano->H[positionH] == 0) {
                positionH--;
            }
            assert(fano->H[positionH] == 1);
            return *this;
        }

        uint64_t operator*() const {
            assert(positionL < fano->count);
            uint64_t l = fano->lowerBits() ? fano->L[positionL] : 0;
            return ((positionH - positionL) << fano->lowerBits()) | l;
        }

        size_t operator-(const ElementPointer &pointer) const {
            return index() - pointer.index();
        }

        size_t index() const {
            return positionL;
        }
    };

    ~EliasFanoM() {
        invalidateSelectDatastructure();
    }

    EliasFanoM(size_t num, uint64_t universeSize)
        : L(elias_fano_lo_width(universeSize, num) == 0 ? 0 : num, 0, elias_fano_lo_width(universeSize, num)),
          H((universeSize >> lowerBits()) + num + 1, false),
          universeSize(universeSize) {
    }

    /**
     * Each index MUST be added exactly once but they can be added without ordering.
     * Either push_back OR add can be called. Combining them is not supported.
     */
    void add(size_t index, uint64_t element) {
        assert(index < L.size() || lowerBits() == 0);
        assert(element < universeSize);
        uint64_t l = element & maskLowerBits();
        uint64_t h = element >> lowerBits();
        assert(element == h * (1l << lowerBits()) + l);
        if (lowerBits() != 0) {
            L[index] = l;
        }
        assert(h + index < H.size());
        H[h + index] = true;
        invalidateSelectDatastructure();
        count++;
    }

    void push_back(uint64_t element) {
        #ifndef NDEBUG
        assert(element >= previousInsert);
        previousInsert = element;
        #endif
        add(count, element);
    }

    void invalidateSelectDatastructure() {
        delete rankSelect;
    }

    /**
     * Returns an ElementPointer to the last stored element that is <= the parameter.
     * When multiple duplicate elements are stored, returns the first occurrence.
     */
    [[nodiscard]] ElementPointer predecessorPosition(uint64_t element) const {
        if (rankSelect == nullptr) {
            throw std::logic_error("Rank/Select not initialized yet. Missing call to buildRankSelect");
        }
        assert(element >= *at(0));

        const uint64_t elementH = element >> lowerBits();
        const uint64_t elementL = element & maskLowerBits();
        uint64_t positionH;
        uint64_t positionL;
        if (elementH == 0) {
            positionH = 0;
            positionL = 0;
        } else {
            positionH = rankSelect->select0(elementH) + 1;
            assert(positionH <= H.size());
            positionL = positionH - elementH;
            assert(positionL <= L.size());
        }
        if (H[positionH] == 0 || positionL == L.size()) {
            // No item with same upper bits stored
            if (positionL > 0) {
                // Return previous item
                positionL--;
                positionH--; // positionH >= positionL, so no underflow
            }
        } else if (lowerBits() != 0) {
            // Look through elements with the same upper bits
            while (true) {
                const uint64_t lower = L[positionL];
                if (lower > elementL) {
                    // Return previous item
                    if (positionL > 0) {
                        positionL--;
                        positionH--; // positionH >= positionL, so no underflow
                    }
                    break;
                } else if (lower == elementL) {
                    // Return first equal item
                    break;
                } else if (H[positionH + 1] == 0) {
                    // End of section. Next item will be larger, so return this.
                    break;
                }
                positionH++;
                positionL++;
            }
        }
        // In case we returned the last item of the previous block, we need to find out its upper bits.
        while (positionH > 0 && H[positionH] == 0) {
            positionH--;
        }
        assert(*at(positionL) <= element);
        assert(positionL == count - 1 || *at(positionL + 1) >= element);
        assert(positionL == 0 || *at(positionL - 1) < element);

        ElementPointer ptr(positionH, positionL, *this);
        #ifndef NDEBUG
        assert(*ptr <= element);
        if (positionL < count - 1) {
            ++ptr;
            assert(*ptr >= element);
            --ptr;
        }
        if (positionL > 0) {
            --ptr;
            assert(*ptr < element);
            ++ptr;
        }
        #endif
        return ptr;
    }

    ElementPointer begin() const {
        size_t positionH = 0;
        while (H[positionH] == 0) {
            positionH++;
        }
        return ElementPointer(positionH, 0, *this);
    }

    [[nodiscard]] ElementPointer at(size_t position) const {
        if (rankSelect == nullptr) {
            throw std::logic_error("Rank/Select not initialized yet. Missing call to buildRankSelect");
        }
        uint64_t positionH = rankSelect->select1(position + 1);
        return ElementPointer(positionH, position, *this);
    }

    [[nodiscard]] uint64_t universe_size() {
        return universeSize;
    }

    void buildRankSelect() {
        if (rankSelect == nullptr) {
            rankSelect = new pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES>(H);
        }
    }

    [[nodiscard]] size_t size() const {
        return L.size();
    }

    /**
     * Space usage of this data structure, in bytes.
     */
    [[nodiscard]] size_t space() const {
        return L.capacity() / 8 + H.size() / 8 + selectStructureOverhead() + sizeof(*this);
    }

    [[nodiscard]] int selectStructureOverhead() const {
        return rankSelect->space_usage();
    }
};

/**
 * Elias-Fano implementation without efficient rank/select.
 */
struct SequentialEliasFano {
    static constexpr uint8_t bit_size_bits = 6;

    /** Returns the bit-width of the low part of an Elias-Fano encoding of the given elements, and the bit offsets to
     * the end of the various components of the encoding: metadata, low part, high part. */
    static std::tuple<uint8_t, size_t, size_t, size_t> encodingOffsets(const std::vector<uint64_t> &elements) {
        auto u = elements.back() + 1;
        auto n = elements.size();
        auto l = elias_fano_lo_width(u, n);
        auto metadata = bit_size_bits;
        auto loBits = l * n;
        auto hiBits = (u >> l) + n + 1;
        return {l, metadata, metadata + loBits, metadata + loBits + hiBits};
    }

    /** Returns the number of bits needed to encode the given elements. */
    static size_t encodingSize(const std::vector<uint64_t> &elements) { return std::get<3>(encodingOffsets(elements)); }

    /** Encodes the given elements into the given data array, starting at the given bit offset. The area between
     * offset and offset + encodingSize(elements) must be already zeroed. */
    static void write(uint64_t *data, size_t offset, const std::vector<uint64_t> &elements) {
        auto [l, metadataEnd, loEnd, hiEnd] = encodingOffsets(elements);
        auto loMask = (uint64_t(1) << l) - 1;

        auto wordPtr = data + offset / 64;
        auto wordOffset = uint8_t(offset % 64);
        sdsl::bits::write_int_and_move(wordPtr, l, wordOffset, bit_size_bits);

        for (size_t i = 0; i < elements.size(); i++) {
            if (l > 0)
                sdsl::bits::write_int_and_move(wordPtr, elements[i] & loMask, wordOffset, l);
            auto hiBitPos = offset + loEnd + (elements[i] >> l) + i;
            data[hiBitPos / 64] |= uint64_t(1) << (hiBitPos % 64);
        }
    }

    class Iterator {
        const uint64_t *data; ///< Pointer to the data array.
        uint8_t l;            ///< Number of bits used for each integer stored in the low part.
        size_t metadataEnd;   ///< Bit offset of the end of the metadata (i.e. start of the low part).
        size_t loEnd;         ///< Bit offset of the end of the low part (i.e. start of the high part).
        size_t n;             ///< Number of elements.
        size_t rank;          ///< Rank/position of element in [0, n) currently pointed by the iterator.
        size_t hiPos;         ///< loEnd + position in the high part of the element currently pointed by the iterator.

    public:

        Iterator(const uint64_t *data, size_t offset, size_t n) : data(data), n(n) {
            auto wordPtr = data + offset / 64;
            auto wordOffset = uint8_t(offset % 64);
            l = sdsl::bits::read_int_and_move(wordPtr, wordOffset, bit_size_bits);
            metadataEnd = offset + bit_size_bits;
            loEnd = offset + bit_size_bits + l * n;
            reset();
        }

        /** Returns the element pointed by the iterator. */
        uint64_t operator*() {
            auto loBitPos = metadataEnd + rank * l;
            auto lo = sdsl::bits::read_int(data + loBitPos / 64, loBitPos % 64, l);
            return (hiPos - loEnd - rank) << l | lo;
        }

        /** Advances the iterator to the next element. */
        Iterator &operator++() {
            ++rank;
            nextHi();
            return *this;
        }

        /** Moves the iterator to the previous element. */
        void operator--() {
            --rank;
            prevHi();
        }

        /** Returns the position of the element pointed by the iterator. */
        [[nodiscard]] size_t index() const { return rank; }

        /** Returns the element at the given position. */
        uint64_t at(size_t pos) {
            reset();
            for (; rank < pos; ++rank)
                nextHi();
            return **this;
        }

        /** Returns the number of elements in the sequence. */
        [[nodiscard]] size_t size() const { return n; }

    private:

        /** Initializes the iterator to the first element. */
        void reset() {
            rank = 0;
            hiPos = loEnd;
            if ((data[hiPos / 64] & uint64_t(1) << (hiPos % 64)) == 0)
                nextHi();
        }

        /** Moves hiPos to the next set bit in the high part */
        void nextHi() {
            auto wordIdx = hiPos / 64;
            uint64_t word = data[wordIdx] & uint64_t(-1) << (hiPos % 64);
            word &= word - 1;
            while (word == 0)
                word = data[++wordIdx];
            hiPos = wordIdx * 64 + __builtin_ctzll(word);
        }

        /** Moves hiPos to the previous set bit in the high part */
        void prevHi() {
            auto wordIdx = hiPos / 64;
            uint64_t word = data[wordIdx] & ((uint64_t(1) << hiPos % 64) - 1);
            while (word == 0)
                word = data[--wordIdx];
            hiPos = wordIdx * 64 + 63 - __builtin_clzll(word);
        }
    };
};

} // Namespace util
