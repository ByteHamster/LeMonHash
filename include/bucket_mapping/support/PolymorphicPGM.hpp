#pragma once

#include "OptimalPiecewiseLinearModel.hpp"
#include "SuccinctPGM.hpp"
#include "util.hpp"
#include <sdsl/bits.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pgm {

struct Segment {
    using K = uint64_t;
    using Floating = float;
    K key;             ///< The first key that the segment indexes.
    Floating slope;    ///< The slope of the segment.
    int32_t intercept; ///< The intercept of the segment.

    Segment() = default;

    Segment(K key, Floating slope, int32_t intercept) : key(key), slope(slope), intercept(intercept) {};

    explicit Segment(size_t n) : key(std::numeric_limits<K>::max()), slope(), intercept(n) {};

    explicit Segment(const typename internal::OptimalPiecewiseLinearModel<K, size_t>::CanonicalSegment &cs)
        : key(cs.first.x) {
        auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(key);
        if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
            throw std::overflow_error("Change the type of Segment::intercept to int64");
        slope = cs_slope;
        intercept = cs_intercept;
    }

    friend inline bool operator<(const K &k, const Segment &s) { return k < s.key; }

    /**
     * Returns the approximate position of the specified key.
     * @param k the key whose position must be approximated
     * @return the approximate position of the specified key
     */
    inline size_t operator()(const K &k) const {
        auto pos = int64_t(slope * (k - key)) + intercept;
        return pos > 0 ? size_t(pos) : 0ull;
    }
};

/**
 * Variant of the PGM that either stores its data uncompressed or compressed, depending on whether compression is
 * beneficial.
 * Also adjusts "rank space" just before a segment start that no key is mapped to.
 * So this does not actually return rank estimates.
 */
class PolymorphicPGMIndex {
    static constexpr auto align_val = std::align_val_t(std::max<size_t>(__STDCPP_DEFAULT_NEW_ALIGNMENT__, 4));
    static constexpr uint8_t first_key_bits = 64;
    static constexpr uint8_t bit_width_bits = 6;
    static constexpr uint8_t slope_bits = 32;

    static constexpr uint8_t one_segment_size_bits = 10;
    static constexpr uint8_t one_segment_intercept_bits = 20;
    static constexpr int64_t max_one_segment_intercept = int64_t(1) << (one_segment_intercept_bits - 1);
    static_assert(slope_bits + one_segment_size_bits + one_segment_intercept_bits == 62);

    using K = uint64_t;

    bool one_segment : 1;
    bool compressed : 1;
    uint64_t raw_ptr : 62;

    // if !one_segment and !compressed, then uncompressed_data() points to a memory area containing:
    //   first_key in first_key_bits
    //   key_bits-1 in bit_width_bits
    //   size_bits-1 in bit_width_bits
    //   size-1 in size_bits
    //   n_segments-1 in size_bits
    //   first segment (slope, intercept) in (slope_bits, size_bits)
    //   remaining segments, where ith segment: (key-first_key-i, slope, intercept) in (key_bits, slope_bits, size_bits)

    [[nodiscard]] uint64_t *uncompressed_data() const {
        assert(!one_segment && !compressed);
        return reinterpret_cast<uint64_t *>(raw_ptr << 2);
    }

    [[nodiscard]] SuccinctPGMIndex *compressed_data() const {
        assert(!one_segment && compressed);
        return reinterpret_cast<SuccinctPGMIndex *>(raw_ptr << 2);
    }

    [[nodiscard]] uint64_t *one_segment_data() const {
        assert(one_segment);
        return reinterpret_cast<uint64_t *>(raw_ptr << 2);
    }

    void free_data() {
        if (!one_segment) {
            if (compressed)
                compressed_data()->~SuccinctPGMIndex();
            else
                operator delete[] (uncompressed_data(), align_val);
        }
    }

    [[nodiscard]] std::tuple<uint32_t, float, int32_t> get_one_segment() const {
        assert(one_segment);
        auto n = uint32_t(raw_ptr >> (one_segment_intercept_bits + slope_bits));
        auto slope = as_float(uint32_t((raw_ptr >> one_segment_intercept_bits) & 0xFFFFFFFF));
        auto intercept = int32_t(raw_ptr & ((1ull << one_segment_intercept_bits) - 1)) - max_one_segment_intercept;
        return {n, slope, intercept};
    }

    [[nodiscard]] std::tuple<uint64_t, uint8_t, uint8_t, size_t, size_t, size_t> metadata() const {
        assert(!one_segment && !compressed);
        const uint64_t *ptr = uncompressed_data() + 1;
        uint8_t offset = 0;
        auto key_bits = sdsl::bits::read_int_and_move(ptr, offset, bit_width_bits) + 1;
        auto size_bits = sdsl::bits::read_int_and_move(ptr, offset, bit_width_bits) + 1;
        auto size = sdsl::bits::read_int_and_move(ptr, offset, size_bits) + 1;
        auto n_segments = sdsl::bits::read_int_and_move(ptr, offset, size_bits) + 1;
        return {uncompressed_data()[0], key_bits, size_bits, size, n_segments, 64 * (ptr - uncompressed_data()) + offset};
    }

    static size_t words_needed(uint8_t key_bits, uint8_t size_bits, size_t n_segments) {
        auto bits = first_key_bits + 2 * bit_width_bits + 2 * size_bits
            + (key_bits + slope_bits + size_bits) * n_segments - key_bits;
        return (bits + 63) / 64;
    }

    [[nodiscard]] uint64_t segment_key_delta(size_t i, size_t segments_offset,
                                             uint8_t segment_bits, uint8_t key_bits) const {
        assert(i > 0 && !one_segment && !compressed);
        segments_offset += i * segment_bits - key_bits;
        return i + sdsl::bits::read_int(uncompressed_data() + segments_offset / 64, segments_offset % 64, key_bits);
    }

    [[nodiscard]] std::tuple<uint64_t, float, int64_t, int64_t> segment(size_t i) const {
        assert(!one_segment && !compressed);
        auto [first_key, key_bits, size_bits, size, n_segments, segments_offset] = metadata();
        segments_offset += i * (key_bits + slope_bits + size_bits) - (i > 0 ? key_bits : 0);
        const uint64_t *ptr = uncompressed_data() + segments_offset / 64;
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

    PolymorphicPGMIndex() : one_segment(false), compressed(false), raw_ptr(0) {}
    PolymorphicPGMIndex(const PolymorphicPGMIndex &other) = delete;
    PolymorphicPGMIndex(PolymorphicPGMIndex &&other) {
        one_segment = other.one_segment;
        compressed = other.compressed;
        raw_ptr = other.raw_ptr;
        other.one_segment = false;
        other.compressed = false;
        other.raw_ptr = 0;
    }

    PolymorphicPGMIndex& operator=(const PolymorphicPGMIndex &other) = delete;

    PolymorphicPGMIndex& operator=(PolymorphicPGMIndex &&other) {
        if (this != &other) {
            free_data();
            one_segment = other.one_segment;
            compressed = other.compressed;
            raw_ptr = other.raw_ptr;
            other.one_segment = false;
            other.compressed = false;
            other.raw_ptr = 0;
        }
        return *this;
    };

    ~PolymorphicPGMIndex() {
        free_data();
    }

    template<typename RandomIt>
    PolymorphicPGMIndex(RandomIt first, RandomIt last, uint8_t epsilon) {
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
                raw_ptr = n << (one_segment_intercept_bits + slope_bits);
                raw_ptr |= uint64_t(as_uint32(slope)) << one_segment_intercept_bits;
                raw_ptr |= int64_t(intercept) + max_one_segment_intercept;
                assert(get_one_segment() == std::make_tuple((uint32_t) n, (float) slope, (int32_t) intercept));
                return;
            }
        }

        auto key_bits = BIT_WIDTH(std::max<uint64_t>(1, segments.back().key - segments.front().key - (segments.size() - 1)));
        auto size_bits = BIT_WIDTH(n - 1);
        auto words = words_needed(key_bits, size_bits, segments.size());
        if (words > 4 * sizeof(SuccinctPGMIndex)) {
            auto *succinct = new SuccinctPGMIndex(first, last, epsilon);
            if (succinct->size_in_bytes() < words * sizeof(uint64_t)) {
                one_segment = false;
                compressed = true;
                raw_ptr = reinterpret_cast<uint64_t>(succinct) >> 2;
                assert(compressed_data() == succinct);
                return;
            } else {
                delete succinct;
            }
        }

        auto ptr = new(align_val) uint64_t[words];
        one_segment = false;
        compressed = false;
        raw_ptr = reinterpret_cast<uint64_t>(ptr) >> 2;
        assert(uncompressed_data() == ptr);

        uint8_t offset = 0;
        sdsl::bits::write_int_and_move(ptr, segments.front().key, offset, first_key_bits);
        sdsl::bits::write_int_and_move(ptr, key_bits - 1, offset, bit_width_bits);
        sdsl::bits::write_int_and_move(ptr, size_bits - 1, offset, bit_width_bits);
        sdsl::bits::write_int_and_move(ptr, std::min(n, (size_t) after_last_point_of_previous_segment) - 1, offset, size_bits);
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
            auto [n, slope, intercept] = get_one_segment();
            auto p = int64_t(std::round(key * slope)) + intercept;
            return std::min<size_t>(p > 0 ? size_t(p) : 0ull, n - 1);
        }
        if (compressed)
            return compressed_data()->approximate_rank(key);

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
            auto [n, slope, intercept] = get_one_segment();
            while (first != last) {
                auto p = int64_t(std::round(*first * slope)) + intercept;
                f(first, std::min<size_t>(p > 0 ? size_t(p) : 0ull, n - 1));
                ++first;
            }
            return;
        }
        if (compressed) {
            compressed_data()->for_each(first, last, f);
            return;
        }

        auto [first_key, key_bits, size_bits, size, n_segments, segments_offset] = metadata();
        auto segment_bits = key_bits + slope_bits + size_bits;
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
        if (one_segment)
            return std::get<0>(get_one_segment());
        if (compressed)
            return compressed_data()->size();
        return std::get<3>(metadata());
    }

    /** Returns the number of segments in the index. */
    [[nodiscard]] size_t segments_count() const {
        if (one_segment)
            return 1;
        if (compressed)
            return compressed_data()->segments_count();
        return std::get<4>(metadata());
    }

    /**  Returns the size of the index in bytes. */
    [[nodiscard]] size_t size_in_bytes() const {
        if (one_segment)
            return sizeof(*this);
        if (compressed)
            return sizeof(*this) + compressed_data()->size_in_bytes();
        auto [_, key_bits, size_bits, _1, n_segments, _2] = metadata();
        return sizeof(*this) + words_needed(key_bits, size_bits, n_segments) * sizeof(uint64_t);
    }
};

}