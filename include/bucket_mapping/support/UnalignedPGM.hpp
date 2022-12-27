#pragma once

#include "OptimalPiecewiseLinearModel.hpp"
#include "ShortPGM.hpp"
#include "util.hpp"
#include <sdsl/bits.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <cassert>
#include <vector>

namespace pgm {

class UnalignedPGMIndex {
    static constexpr uint8_t first_key_bits = 64;
    static constexpr uint8_t bit_width_bits = 6;
    static constexpr uint8_t slope_bits = 32;
    using K = uint64_t;

    uint64_t *data;

    // data layout:
    //   first_key in first_key_bits
    //   key_bits in bit_width_bits
    //   size_bits in bit_width_bits
    //   size in size_bits
    //   n_segments in size_bits
    //   remaining segments (key - first_key, slope, intercept) in (key_bits, slope_bits, size_bits)

    [[nodiscard]] std::tuple<uint64_t, uint8_t, uint8_t, size_t, size_t, size_t> metadata() const {
        const uint64_t *ptr = data + 1;
        uint8_t offset = 0;
        auto key_bits = sdsl::bits::read_int_and_move(ptr, offset, bit_width_bits) + 1;
        auto size_bits = sdsl::bits::read_int_and_move(ptr, offset, bit_width_bits) + 1;
        auto size = sdsl::bits::read_int_and_move(ptr, offset, size_bits);
        auto n_segments = sdsl::bits::read_int_and_move(ptr, offset, size_bits);
        return {data[0], key_bits, size_bits, size, n_segments, 64 * (ptr - data) + offset};
    }

    [[nodiscard]] uint64_t segment_key_delta(size_t i, size_t segments_offset,
                                             uint8_t segment_bits, uint8_t key_bits) const {
        segments_offset += i * segment_bits;
        return sdsl::bits::read_int(data + segments_offset / 64, segments_offset % 64, key_bits);
    }

    [[nodiscard]] std::tuple<uint64_t, float, int64_t, int64_t> segment(size_t i) const {
        auto [first_key, key_bits, size_bits, size, n_segments, segments_offset] = metadata();
        segments_offset += i * (key_bits + slope_bits + size_bits);
        const uint64_t *ptr = data + segments_offset / 64;
        uint8_t offset = segments_offset % 64;

        auto key = sdsl::bits::read_int_and_move(ptr, offset, key_bits) + first_key;
        auto slope = as_float(sdsl::bits::read_int_and_move(ptr, offset, slope_bits));
        auto intercept = int64_t(sdsl::bits::read_int_and_move(ptr, offset, size_bits));

        auto next_intercept = int64_t(size);
        if (i + 1 < n_segments) {
            sdsl::bits::move_right(ptr, offset, key_bits);
            sdsl::bits::move_right(ptr, offset, slope_bits);
            next_intercept = int64_t(sdsl::bits::read_int(ptr, offset, size_bits));
        }

        return {key, slope, intercept, next_intercept};
    }

    static float as_float(uint32_t i) {
        float f;
        memcpy(&f, &i, sizeof(uint32_t));
        return f;
    }

    static uint32_t as_uint32(float f) {
        uint32_t i;
        memcpy(&i, &f, sizeof(float));
        return i;
    }

public:

    UnalignedPGMIndex() : data(nullptr) {}
    UnalignedPGMIndex(const UnalignedPGMIndex&) = delete;
    UnalignedPGMIndex& operator=(const UnalignedPGMIndex &) = delete;
    UnalignedPGMIndex(UnalignedPGMIndex&&) = delete;

    UnalignedPGMIndex& operator=(UnalignedPGMIndex &&other) {
        if (this != &other) {
            delete[] data;
            data = other.data;
            other.data = nullptr;
        }
        return *this;
    };

    ~UnalignedPGMIndex() {
        delete[] data;
    }

    template<typename RandomIt>
    UnalignedPGMIndex(RandomIt first, RandomIt last, uint8_t epsilon) {
        auto n = (size_t) std::distance(first, last);
        assert(n > 0);

        std::vector<Segment> segments;
        segments.reserve(n / std::max(2, epsilon * epsilon));

        auto in_fun = [&](auto i) { return std::pair<K, size_t>(first[i], i); };
        auto out_fun = [&](auto cs, auto) { segments.emplace_back(cs); };
        internal::make_segmentation_mod(n, epsilon, in_fun, out_fun, false);

        auto key_bits = BIT_WIDTH(std::max<uint64_t>(1, segments.back().key - segments.front().key));
        auto size_bits = BIT_WIDTH(n);
        auto bit_size = first_key_bits + 2 * bit_width_bits + 2 * size_bits + (key_bits + slope_bits + size_bits) * segments.size();
        auto words = (bit_size + 63) / 64;
        data = new uint64_t[words];

        auto ptr = data;
        uint8_t offset = 0;
        sdsl::bits::write_int_and_move(ptr, segments.front().key, offset, first_key_bits);
        sdsl::bits::write_int_and_move(ptr, key_bits - 1, offset, bit_width_bits);
        sdsl::bits::write_int_and_move(ptr, size_bits - 1, offset, bit_width_bits);
        sdsl::bits::write_int_and_move(ptr, n, offset, size_bits);
        sdsl::bits::write_int_and_move(ptr, segments.size(), offset, size_bits);

        for (auto &s : segments) {
            assert(s.intercept >= 0);
            sdsl::bits::write_int_and_move(ptr, s.key - segments[0].key, offset, key_bits);
            sdsl::bits::write_int_and_move(ptr, as_uint32(s.slope), offset, slope_bits);
            sdsl::bits::write_int_and_move(ptr, s.intercept, offset, size_bits);
        }
    }


    /** Returns the approximate rank of @p key. */
    [[nodiscard]] size_t approximate_rank(const K &key) const {
        auto [first_key, key_bits, size_bits, size, n_segments, segments_offset] = metadata();
        if (key < first_key)
            return 0;

        auto segment_bits = key_bits + slope_bits + size_bits;
        auto key_delta = key - first_key;
        size_t i = 1;

        if (n_segments > 8) {
            auto count = n_segments - 1;
            while (count > 0) {
                auto step = count / 2;
                auto mid = i + step;
                if (key_delta >= segment_key_delta(mid, segments_offset, segment_bits, key_bits)) {
                    i = mid + 1;
                    count -= step + 1;
                } else {
                    count = step;
                }
            }
        } else {
            while (i < n_segments && key_delta >= segment_key_delta(i, segments_offset, segment_bits, key_bits))
                ++i;
        }

        auto [segment_key, slope, intercept, next_intercept] = segment(i - 1);
        auto pos = int64_t(slope * (key - segment_key)) + intercept;
        return std::min<size_t>(pos > 0 ? size_t(pos) : 0ull, next_intercept - 1);
    }

    /** Execute a given function for each key in the sorted range [first, last). The function takes as the argument
     * an iterator to the current key and the corresponding approximate rank computed by the index. */
    template<typename ForwardIt, typename F>
    void for_each(ForwardIt first, ForwardIt last, F f) const {
        auto [first_key, key_bits, size_bits, size, n_segments, segments_offset] = metadata();
        auto segment_bits = key_bits + slope_bits + size_bits;
        auto [segment_key, slope, intercept, next_intercept] = segment(0);
        auto next_segment_key = first_key + segment_key_delta(1, segments_offset, segment_bits, key_bits);
        size_t segment_i = 0;
        auto it = first;
        while (it != last) {
            if (segment_i + 1 != n_segments && *it >= next_segment_key) {
                ++segment_i;
                std::tie(segment_key, slope, intercept, next_intercept) = segment(segment_i);
                if (segment_i + 1 != n_segments)
                    next_segment_key = first_key + segment_key_delta(segment_i + 1, segments_offset, segment_bits, key_bits);
            }
            auto pos = int64_t(slope * (*it - segment_key)) + intercept;
            f(it, std::min<size_t>(pos > 0 ? size_t(pos) : 0ull, next_intercept - 1));
            ++it;
        }
    }

    /** Returns the number of elements the index was built on. */
    [[nodiscard]] size_t size() const {
        return std::get<3>(metadata());
    }

    /**  Returns the size of the index in bytes. */
    [[nodiscard]] size_t size_in_bytes() const {
        auto [first_key, key_bits, size_bits, size, n_segments, segments_offset] = metadata();
        auto bit_size = first_key_bits + 2 * bit_width_bits + 2 * size_bits + (key_bits + slope_bits + size_bits) * n_segments;
        auto words = (bit_size + 63) / 64;
        return sizeof(*this) + words * sizeof(uint64_t);
    }
};

}