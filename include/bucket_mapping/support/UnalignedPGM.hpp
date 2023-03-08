#pragma once

#include "OptimalPiecewiseLinearModel.hpp"
#include <sdsl/bits.hpp>
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <cassert>
#include <vector>

namespace pgm {

/**
 * Variant of the PGM that stores its data in an uncompressed bit vector.
 * Also compresses "rank space" just before a segment start that no key is mapped to.
 * So this does not actually return rank estimates.
 */
class UnalignedPGMIndex {
    static constexpr uint8_t first_key_bits = 64;
    static constexpr uint8_t bit_width_bits = 6;
    static constexpr uint8_t slope_bits = 32;

    static constexpr uint8_t one_segment_size_bits = 10;
    static constexpr uint8_t one_segment_intercept_bits = 21;
    static constexpr int64_t max_one_segment_intercept = int64_t(1) << (one_segment_intercept_bits - 1);
    static_assert(slope_bits + one_segment_size_bits + one_segment_intercept_bits == 63);

    using K = uint64_t;

    bool one_segment = false;
    uint64_t *raw_ptr;
    bool createdFromRawPtr = false;

    struct OneSegmentData {
        uint32_t _one_segment : 1 = 1;
        uint32_t n : one_segment_size_bits;
        uint32_t intercept : one_segment_intercept_bits;
        float slope;
    };
    static_assert(sizeof(OneSegmentData) == sizeof(uint64_t));
    OneSegmentData one_segment_data;

    // if !one_segment, then raw_ptr << 1 points to a memory area containing:
    //   first_key in first_key_bits
    //   key_bits-1 in bit_width_bits
    //   size_bits-1 in bit_width_bits
    //   intercept_bits-1 in bit_width_bits
    //   size-1 in size_bits
    //   n_segments-1 in size_bits
    //   first segment (slope, intercept) in (slope_bits, intercept_bits)
    //   remaining segments, where ith segment: (key-first_key-i, slope, intercept-i) in (key_bits, slope_bits, intercept_bits)

    [[nodiscard]] uint64_t *data() const {
        assert(!one_segment);
        return raw_ptr;
    }

    void free_data() {
        if (!one_segment && !createdFromRawPtr)
            delete[] data();
    }

    [[nodiscard]] std::tuple<uint64_t, uint8_t, uint8_t, uint8_t, size_t, size_t, size_t> metadata() const {
        assert(!one_segment);
        const uint64_t *ptr = data();
        uint8_t offset = 1;
        auto first_key = sdsl::bits::read_int_and_move(ptr, offset, first_key_bits);
        auto key_bits = sdsl::bits::read_int_and_move(ptr, offset, bit_width_bits) + 1;
        auto size_bits = sdsl::bits::read_int_and_move(ptr, offset, bit_width_bits) + 1;
        auto intercept_bits = sdsl::bits::read_int_and_move(ptr, offset, bit_width_bits) + 1;
        auto size = sdsl::bits::read_int_and_move(ptr, offset, size_bits) + 1;
        auto n_segments = sdsl::bits::read_int_and_move(ptr, offset, size_bits) + 1;
        return {first_key, key_bits, size_bits, intercept_bits, size, n_segments, 64 * (ptr - data()) + offset};
    }

    static size_t bytes_needed(uint8_t key_bits, uint8_t size_bits, uint8_t intercept_bits, size_t n_segments) {
        auto bits = 1 + first_key_bits + 3 * bit_width_bits + 2 * size_bits
            + (key_bits + slope_bits + intercept_bits) * n_segments - key_bits;
        return (bits + 7) / 8;
    }

    [[nodiscard]] uint64_t segment_key_delta(size_t i, size_t segments_offset,
                                             uint8_t segment_bits, uint8_t key_bits) const {
        assert(i > 0);
        segments_offset += i * segment_bits - key_bits;
        return i + sdsl::bits::read_int(data() + segments_offset / 64, segments_offset % 64, key_bits);
    }

    [[nodiscard]] std::tuple<uint64_t, float, int64_t, int64_t> segment(size_t i) const {
        auto [first_key, key_bits, size_bits, intercept_bits, size, n_segments, segments_offset] = metadata();
        segments_offset += i * (key_bits + slope_bits + intercept_bits) - (i > 0 ? key_bits : 0);
        const uint64_t *ptr = data() + segments_offset / 64;
        uint8_t offset = segments_offset % 64;

        auto key = first_key;
        if (i > 0)
            key += i + sdsl::bits::read_int_and_move(ptr, offset, key_bits);
        auto slope = as_float(sdsl::bits::read_int_and_move(ptr, offset, slope_bits));
        auto intercept = int64_t(sdsl::bits::read_int_and_move(ptr, offset, intercept_bits)) + i;

        auto next_intercept = int64_t(size);
        if (i + 1 < n_segments) {
            sdsl::bits::move_right(ptr, offset, key_bits);
            sdsl::bits::move_right(ptr, offset, slope_bits);
            next_intercept = int64_t(sdsl::bits::read_int_and_move(ptr, offset, intercept_bits)) + i + 1;
        }
        assert((ptr - data()) * 8 < size_in_bytes()
                || ((ptr - data()) * 8 == size_in_bytes() && offset <= 8 * (size_in_bytes() % 8)));

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
            createdFromRawPtr = other.createdFromRawPtr;
            one_segment_data = other.one_segment_data;
            raw_ptr = other.raw_ptr;
            other.one_segment = false;
            other.raw_ptr = nullptr;
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

        bool skip_first = n > 2; // the first key will always be mapped to rank 0, so we skip it in the segmentation

        using segment_t = internal::OptimalPiecewiseLinearModel<uint64_t, size_t>::CanonicalSegment;
        segment_t first_segment;
        int32_t after_last_point_of_previous_segment = 0;
        auto in_fun = [&](auto i) { return std::pair<K, size_t>(first[i + skip_first], i + skip_first); };
        auto out_fun = [&](segment_t cs, auto lastPointOfSegment) {
            if (segments.empty())
                first_segment = cs;
            segments.emplace_back(cs);
            Segment &new_segment = segments.back();
            new_segment.intercept = std::min(after_last_point_of_previous_segment, new_segment.intercept);
            after_last_point_of_previous_segment = int64_t(std::round(new_segment.slope * (lastPointOfSegment.first - new_segment.key))) + new_segment.intercept + 1;
        };
        internal::make_segmentation_mod(n - skip_first, epsilon, in_fun, out_fun, false);

        if (segments.size() > 1 && segments.back().key == *(last - 1)) {
            // if last segment covers only one key we can remove it, as the previous one can map that key correctly
            segments.pop_back();
        }

        if (segments.size() == 1 && n < 1 << one_segment_size_bits) {
            auto [slope, intercept] = first_segment.get_floating_point_segment(0);
            if (intercept > -max_one_segment_intercept && intercept < max_one_segment_intercept) {
                one_segment = true;
                one_segment_data.n = n;
                one_segment_data.intercept = intercept;
                one_segment_data.slope = slope;
                return;
            }
        }

        auto key_bits = std::bit_width(std::max<uint64_t>(1, segments.back().key - segments.front().key - (segments.size() - 1)));
        auto size_bits = std::bit_width(n - 1);
        auto intercept_bits = std::bit_width(std::max<uint64_t>(1, segments.back().intercept - (segments.size() - 1)));
        size_t bytes = bytes_needed(key_bits, size_bits, intercept_bits, segments.size());
        raw_ptr = new uint64_t[(bytes + 7) / 8];
        auto ptr = raw_ptr;
        one_segment = false;

        uint8_t offset = 0;
        sdsl::bits::write_int_and_move(ptr, false, offset, 1); // Not one_segment
        sdsl::bits::write_int_and_move(ptr, segments.front().key, offset, first_key_bits);
        sdsl::bits::write_int_and_move(ptr, key_bits - 1, offset, bit_width_bits);
        sdsl::bits::write_int_and_move(ptr, size_bits - 1, offset, bit_width_bits);
        sdsl::bits::write_int_and_move(ptr, intercept_bits - 1, offset, bit_width_bits);
        sdsl::bits::write_int_and_move(ptr, std::min(n, (size_t) after_last_point_of_previous_segment) - 1, offset, size_bits);
        sdsl::bits::write_int_and_move(ptr, segments.size() - 1, offset, size_bits);

        for (size_t i = 0; i < segments.size(); ++i) {
            assert(segments[i].intercept >= 0);
            if (i > 0)
                sdsl::bits::write_int_and_move(ptr, segments[i].key - segments[0].key - i, offset, key_bits);
            sdsl::bits::write_int_and_move(ptr, as_uint32(segments[i].slope), offset, slope_bits);
            sdsl::bits::write_int_and_move(ptr, segments[i].intercept - i, offset, intercept_bits);
        }
        assert((ptr - raw_ptr) * 8 < size_in_bytes()
                || (ptr - raw_ptr) * 8 == size_in_bytes() && offset <= 8 * (size_in_bytes() % 8));
    }

    /** Returns the approximate rank of @p key. */
    [[nodiscard]] size_t approximate_rank(const K &key) const {
        if (one_segment) {
            auto p = int64_t(std::round(key * one_segment_data.slope)) + one_segment_data.intercept;
            return std::min<size_t>(p > 0 ? size_t(p) : 0ull, one_segment_data.n - 1);
        }

        auto [first_key, key_bits, size_bits, intercept_bits, size, n_segments, segments_offset] = metadata();
        if (key < first_key)
            return 0;

        auto segment_bits = key_bits + slope_bits + intercept_bits;
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
            while (first != last) {
                auto p = int64_t(std::round(*first * one_segment_data.slope)) + one_segment_data.intercept;
                f(first, std::min<size_t>(p > 0 ? size_t(p) : 0ull, one_segment_data.n - 1));
                ++first;
            }
            return;
        }

        auto [first_key, key_bits, size_bits, intercept_bits, size, n_segments, segments_offset] = metadata();
        auto segment_bits = key_bits + slope_bits + intercept_bits;
        auto [segment_key, slope, intercept, next_intercept] = segment(0);
        auto next_segment_key = first_key + segment_key_delta(1, segments_offset, segment_bits, key_bits);
        size_t segment_i = 0;
        auto it = first;

        while (it != last && *it < first_key) {
            f(it, 0);
            ++it;
        }

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

    /** Returns the number of output buckets. (Compressed version of number of input objects) */
    [[nodiscard]] size_t size() const {
        return one_segment ? one_segment_data.n : std::get<4>(metadata());
    }

    /** Returns the number of segments in the index. */
    [[nodiscard]] size_t segments_count() const {
        return one_segment ? 1 : std::get<5>(metadata());
    }

    /**  Returns the size of the index in bytes. */
    [[nodiscard]] size_t size_in_bytes() const {
        if (one_segment)
            return sizeof(OneSegmentData);
        auto [_, key_bits, size_bits, intercept_bits, _1, n_segments, _2] = metadata();
        return bytes_needed(key_bits, size_bits, intercept_bits, n_segments);
    }

    void copyTo(char *ptr) {
        if (one_segment) {
            *((OneSegmentData *) ptr) = one_segment_data;
        } else {
            memcpy(ptr, (void *)data(), size_in_bytes());
            assert(!((OneSegmentData *) ptr)->_one_segment);
        }
    }

    explicit UnalignedPGMIndex(const char *ptr) {
        createdFromRawPtr = true;
        if (((OneSegmentData *) ptr)->_one_segment) {
            one_segment = true;
            one_segment_data = *((OneSegmentData *) ptr);
        } else {
            one_segment = false;
            raw_ptr = (uint64_t *) ptr;
        }
    }
};

}