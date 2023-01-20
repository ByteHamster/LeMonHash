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
    static constexpr auto align_val = std::align_val_t(std::max<size_t>(__STDCPP_DEFAULT_NEW_ALIGNMENT__, 2));
    static constexpr uint8_t first_key_bits = 64;
    static constexpr uint8_t bit_width_bits = 6;
    static constexpr uint8_t slope_bits = 32;
    using K = uint64_t;

    bool one_segment : 1;
    uint64_t raw_ptr : 63;

    // if !one_segment, then raw_ptr << 1 points to a memory area containing:
    //   first_key in first_key_bits
    //   key_bits-1 in bit_width_bits
    //   size_bits-1 in bit_width_bits
    //   size-1 in size_bits
    //   n_segments-1 in size_bits
    //   first segment (slope, intercept) in (slope_bits, size_bits)
    //   remaining segments, where ith segment: (key-first_key-i, slope, intercept) in (key_bits, slope_bits, size_bits)

    [[nodiscard]] uint64_t *data() const {
        assert(!one_segment);
        return reinterpret_cast<uint64_t *>(raw_ptr << 1);
    }

    void free_data() {
        if (!one_segment)
            operator delete[] (data(), align_val);
    }

    [[nodiscard]] std::pair<float, uint32_t> get_one_segment() const {
        assert(one_segment);
        uint32_t n = raw_ptr >> 32;
        float slope = as_float(uint32_t(raw_ptr & 0xFFFFFFFF));
        return {slope, n};
    }

    [[nodiscard]] std::tuple<uint64_t, uint8_t, uint8_t, size_t, size_t, size_t> metadata() const {
        assert(!one_segment);
        const uint64_t *ptr = data() + 1;
        uint8_t offset = 0;
        auto key_bits = sdsl::bits::read_int_and_move(ptr, offset, bit_width_bits) + 1;
        auto size_bits = sdsl::bits::read_int_and_move(ptr, offset, bit_width_bits) + 1;
        auto size = sdsl::bits::read_int_and_move(ptr, offset, size_bits) + 1;
        auto n_segments = sdsl::bits::read_int_and_move(ptr, offset, size_bits) + 1;
        return {data()[0], key_bits, size_bits, size, n_segments, 64 * (ptr - data()) + offset};
    }

    static size_t words_needed(uint8_t key_bits, uint8_t size_bits, size_t n_segments) {
        auto bits = first_key_bits + 2 * bit_width_bits + 2 * size_bits
            + (key_bits + slope_bits + size_bits) * n_segments - key_bits;
        return (bits + 63) / 64;
    }

    [[nodiscard]] uint64_t segment_key_delta(size_t i, size_t segments_offset,
                                             uint8_t segment_bits, uint8_t key_bits) const {
        assert(i > 0);
        segments_offset += i * segment_bits - key_bits;
        return i + sdsl::bits::read_int(data() + segments_offset / 64, segments_offset % 64, key_bits);
    }

    [[nodiscard]] std::tuple<uint64_t, float, int64_t, int64_t> segment(size_t i) const {
        auto [first_key, key_bits, size_bits, size, n_segments, segments_offset] = metadata();
        segments_offset += i * (key_bits + slope_bits + size_bits) - (i > 0 ? key_bits : 0);
        const uint64_t *ptr = data() + segments_offset / 64;
        uint8_t offset = segments_offset % 64;

        auto key = first_key;
        if (i > 0)
            key += i + sdsl::bits::read_int_and_move(ptr, offset, key_bits);
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

    UnalignedPGMIndex() : one_segment(false), raw_ptr(0) {}
    UnalignedPGMIndex(const UnalignedPGMIndex &other) = delete;
    UnalignedPGMIndex(UnalignedPGMIndex &&other) {
        one_segment = other.one_segment;
        raw_ptr = other.raw_ptr;
        other.one_segment = false;
        other.raw_ptr = 0;
    }

    UnalignedPGMIndex& operator=(const UnalignedPGMIndex &other) = delete;

    UnalignedPGMIndex& operator=(UnalignedPGMIndex &&other) {
        if (this != &other) {
            free_data();
            one_segment = other.one_segment;
            raw_ptr = other.raw_ptr;
            other.one_segment = false;
            other.raw_ptr = 0;
        }
        return *this;
    };

    ~UnalignedPGMIndex() {
        free_data();
    }

    template<typename RandomIt>
    UnalignedPGMIndex(RandomIt first, RandomIt last, uint8_t epsilon) {
        auto n = (size_t) std::distance(first, last);
        assert(n >= 2);

        std::vector<Segment> segments;
        segments.reserve(n / std::max(2, epsilon * epsilon));

        internal::OptimalPiecewiseLinearModel<uint64_t, size_t>::CanonicalSegment first_segment;
        auto in_fun = [&](auto i) { return std::pair<K, size_t>(first[i], i); };
        auto out_fun = [&](auto cs, auto) {
            if (segments.empty())
                first_segment = cs;
            segments.emplace_back(cs);
        };
        internal::make_segmentation_mod(n, epsilon, in_fun, out_fun, false);

        if (segments.size() == 1 && n < (1ull << 31)) {
            auto [flag, slope] = first_segment.get_segment_through_zero();
            if (flag) {
                one_segment = true;
                raw_ptr = n << 32 | as_uint32(slope);
                assert(get_one_segment() == std::make_pair((float) slope, (uint32_t) n));
                return;
            }
        }

        auto key_bits = BIT_WIDTH(std::max<uint64_t>(1, segments.back().key - segments.front().key - (segments.size() - 1)));
        auto size_bits = BIT_WIDTH(n - 1);
        auto ptr = new(align_val) uint64_t[words_needed(key_bits, size_bits, segments.size())];
        one_segment = false;
        raw_ptr = reinterpret_cast<uint64_t>(ptr) >> 1;
        assert(reinterpret_cast<uint64_t>(ptr) == raw_ptr << 1);

        uint8_t offset = 0;
        sdsl::bits::write_int_and_move(ptr, segments.front().key, offset, first_key_bits);
        sdsl::bits::write_int_and_move(ptr, key_bits - 1, offset, bit_width_bits);
        sdsl::bits::write_int_and_move(ptr, size_bits - 1, offset, bit_width_bits);
        sdsl::bits::write_int_and_move(ptr, n - 1, offset, size_bits);
        sdsl::bits::write_int_and_move(ptr, segments.size() - 1, offset, size_bits);

        for (size_t i = 0; i < segments.size(); ++i) {
            assert(segments[i].intercept >= 0);
            if (i > 0)
                sdsl::bits::write_int_and_move(ptr, segments[i].key - segments[0].key - i, offset, key_bits);
            sdsl::bits::write_int_and_move(ptr, as_uint32(segments[i].slope), offset, slope_bits);
            sdsl::bits::write_int_and_move(ptr, segments[i].intercept, offset, size_bits);
        }
    }

    /** Returns the approximate rank of @p key. */
    [[nodiscard]] size_t approximate_rank(const K &key) const {
        if (one_segment) {
            auto [slope, n] = get_one_segment();
            auto p = int64_t(std::round(key * slope));
            return std::min<size_t>(p > 0 ? size_t(p) : 0ull, n - 1);
        }

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
        auto pos = int64_t(std::round(slope * (key - segment_key))) + intercept;
        return std::min<size_t>(pos > 0 ? size_t(pos) : 0ull, next_intercept - 1);
    }

    /** Execute a given function for each key in the sorted range [first, last). The function takes as the argument
     * an iterator to the current key and the corresponding approximate rank computed by the index. */
    template<typename ForwardIt, typename F>
    void for_each(ForwardIt first, ForwardIt last, F f) const {
        if (one_segment) {
            auto [slope, n] = get_one_segment();
            while (first != last) {
                auto p = int64_t(std::round(*first * slope));
                f(first, std::min<size_t>(p > 0 ? size_t(p) : 0ull, n - 1));
                ++first;
            }
            return;
        }

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
            auto pos = int64_t(std::round(slope * (*it - segment_key))) + intercept;
            f(it, std::min<size_t>(pos > 0 ? size_t(pos) : 0ull, next_intercept - 1));
            ++it;
        }
    }

    /** Returns the number of elements the index was built on. */
    [[nodiscard]] size_t size() const {
        return one_segment ? get_one_segment().second : std::get<3>(metadata());
    }

    /** Returns the number of segments in the index. */
    [[nodiscard]] size_t segments_count() const {
        return one_segment ? 1 : std::get<4>(metadata());
    }

    /**  Returns the size of the index in bytes. */
    [[nodiscard]] size_t size_in_bytes() const {
        if (one_segment)
            return sizeof(*this);
        auto [_, key_bits, size_bits, _1, n_segments, _2] = metadata();
        return sizeof(*this) + words_needed(key_bits, size_bits, n_segments) * sizeof(uint64_t);
    }
};

}