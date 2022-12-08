#pragma once

#include <vector>
#include <cassert>
#include <sdsl/bit_vectors.hpp>
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>

namespace util {
/**
 * Compressed, monotone integer array.
 * Commonly used to store the positions of 1-bits in sparse bit vectors.
 */
class EliasFanoM {
private:
    sdsl::int_vector<> L;
    using ConstIntVector = const sdsl::int_vector<>;
    pasta::BitVector H;
    size_t count = 0;
    size_t universeSize = 0;
    pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES> *rankSelect = nullptr;

#ifndef NDEBUG
    uint64_t previousInsert = 0;
#endif

    [[nodiscard]] uint8_t lowerBits() const { return L.width(); }

    [[nodiscard]] uint64_t maskLowerBits() const { return sdsl::bits::lo_set[lowerBits()]; }

    // Returns `L.width()`.
    //
    // This is based on:
    //
    //   Ma, Puglisi, Raman, Zhukova:
    //   On Elias-Fano for Rank Queries in FM-Indexes.
    //   DCC 2021.
    //
    // Implementation credit: Jouni Siren, https://github.com/vgteam/sdsl-lite
    static uint8_t get_params(size_t universe, size_t ones) {
        uint8_t low_width = 1;
        // Multisets with too many ones will have width 1.
        if (ones > 0 && ones <= universe) {
            double ideal_width = std::log2((static_cast<double>(universe) * std::log(2.0)) / static_cast<double>(ones));
            low_width = (uint8_t) std::round(std::max(ideal_width, 1.0));
        }
        return low_width;
    }

public:

    /**
     * Efficient pointer into an Elias-Fano coded sequence.
     * When incrementing/decrementing and reading, no additional select query is performed.
     */
    struct ElementPointer {
    private:
        size_t positionL;
        size_t positionH;
        size_t h;
        const EliasFanoM *fano;

    public:
        ElementPointer(size_t h, size_t positionH, size_t positionL, const EliasFanoM &fano)
            : positionL(positionL), positionH(positionH), h(h), fano(&fano) {
            assert(fano.H[positionH] == 1);
        }

        ElementPointer& operator++() {
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
                h++;
            }
            assert(fano->H[positionH] == 1);
            return *this;
        }

        ElementPointer& operator--() {
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
                h--;
            }
            assert(fano->H[positionH] == 1);
            return *this;
        }

        uint64_t operator *() const {
            assert(positionL < fano->count);
            if (fano->lowerBits() == 0) {
                return h;
            }
            uint64_t l = static_cast<ConstIntVector&>(fano->L)[positionL];
            return (h << fano->lowerBits()) + l;
        }

        size_t operator -(const ElementPointer &pointer) const {
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
        : L(get_params(universeSize, num) == 0 ? 0 : num, 0, get_params(universeSize, num)),
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
        assert(element == h*(1l << lowerBits()) + l);
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
                const uint64_t lower = static_cast<ConstIntVector&>(L)[positionL];
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
        uint64_t resultH = elementH;
        while (positionH > 0 && H[positionH] == 0) {
            positionH--;
            resultH--;
        }
        assert(*at(positionL) <= element);
        assert(positionL == count - 1 || *at(positionL + 1) >= element);
        assert(positionL == 0 || *at(positionL - 1) < element);

        ElementPointer ptr(resultH, positionH, positionL, *this);
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
        size_t h = 0;
        size_t positionH = 0;
        while (H[positionH] == 0) {
            positionH++;
            h++;
        }
        return ElementPointer(h, positionH, 0, *this);
    }

    [[nodiscard]] ElementPointer at(size_t position) const {
        if (rankSelect == nullptr) {
            throw std::logic_error("Rank/Select not initialized yet. Missing call to buildRankSelect");
        }
        uint64_t positionH = rankSelect->select1(position + 1);
        uint64_t h = positionH - position;
        return ElementPointer(h, positionH, position, *this);
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
        return L.capacity()/8 + H.size()/8 + selectStructureOverhead() + sizeof(*this);
    }

    [[nodiscard]] int selectStructureOverhead() const {
        return rankSelect->space_usage();
    }
};

} // Namespace util
