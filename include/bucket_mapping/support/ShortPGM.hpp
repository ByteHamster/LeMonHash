#pragma once

#include "OptimalPiecewiseLinearModel.hpp"
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

#pragma pack(push, 1)

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

template<size_t Capacity=0>
class ShortPGMIndex {
    static constexpr size_t max_local_capacity = 15;

    struct Local {
        Segment segments[Capacity];

        Local(size_t) {}
        auto begin() const { return segments; }
        auto end() const { return segments + Capacity; }
        auto size() const { return Capacity; }
        auto capacity() const { return Capacity; }
        auto size_in_bytes() const { return 0; }
    };

    struct Heap {
        static constexpr uint8_t segments_size_bits = 24;
        size_t segments_size : segments_size_bits;
        Segment *segments;

        Heap(size_t size) : segments_size(size), segments(new Segment[size]) {}
        ~Heap() { delete[] segments; }
        auto begin() const { return segments; }
        auto end() const { return segments + segments_size; }
        auto size() const { return segments_size; }
        auto capacity() const { return (1ull << segments_size_bits) - 1; }
        auto size_in_bytes() const { return sizeof(Segment) * segments_size; }
    };

    template<size_t... Is>
    static auto to_tuple_helper(std::integer_sequence<size_t, Is...>) -> std::tuple<ShortPGMIndex<Is>...>;

    using K = uint64_t;
    using Storage = std::conditional_t<(Capacity == 0 || Capacity > max_local_capacity), Heap, Local>;
    Storage segments;
    uint8_t epsilon;

public:

    static constexpr size_t capacity = Capacity;

    using capacities_sequence = std::make_index_sequence<max_local_capacity + 1>;
    using all_types = decltype(to_tuple_helper(capacities_sequence{}));

    static std::pair<void *, size_t> make(auto first, auto last, uint8_t epsilon) {
        return make_helper(first, last, epsilon, capacities_sequence{});
    }

    template <typename F>
    static void visit(F f, void *ptr, size_t type) {
        switch (type) {
            case 0:  f(static_cast<std::tuple_element_t<0,  all_types>*>(ptr)); break;
            case 1:  f(static_cast<std::tuple_element_t<1,  all_types>*>(ptr)); break;
            case 2:  f(static_cast<std::tuple_element_t<2,  all_types>*>(ptr)); break;
            case 3:  f(static_cast<std::tuple_element_t<3,  all_types>*>(ptr)); break;
            case 4:  f(static_cast<std::tuple_element_t<4,  all_types>*>(ptr)); break;
            case 5:  f(static_cast<std::tuple_element_t<5,  all_types>*>(ptr)); break;
            case 6:  f(static_cast<std::tuple_element_t<6,  all_types>*>(ptr)); break;
            case 7:  f(static_cast<std::tuple_element_t<7,  all_types>*>(ptr)); break;
            case 8:  f(static_cast<std::tuple_element_t<8,  all_types>*>(ptr)); break;
            case 9:  f(static_cast<std::tuple_element_t<9,  all_types>*>(ptr)); break;
            case 10: f(static_cast<std::tuple_element_t<10, all_types>*>(ptr)); break;
            case 11: f(static_cast<std::tuple_element_t<11, all_types>*>(ptr)); break;
            case 12: f(static_cast<std::tuple_element_t<12, all_types>*>(ptr)); break;
            case 13: f(static_cast<std::tuple_element_t<13, all_types>*>(ptr)); break;
            case 14: f(static_cast<std::tuple_element_t<14, all_types>*>(ptr)); break;
            case 15: f(static_cast<std::tuple_element_t<15, all_types>*>(ptr)); break;
            default: f(static_cast<ShortPGMIndex<>*>(ptr));
        }
    }

    ShortPGMIndex() = default;

    ShortPGMIndex(const std::vector<Segment> &s, uint8_t epsilon) : segments(s.size()), epsilon(epsilon) {
        assert(s.size() <= segments.capacity());
        std::copy(s.begin(), s.end(), &segments.segments[0]);
    }

    /** Returns the approximate rank of @p key. */
    [[nodiscard]] size_t approximate_rank(const K &key, size_t n) const {
        auto k = std::max(segments.begin()->key, key);
        decltype(segments.begin()) it;
        if constexpr (Capacity > 0 && Capacity <= 7) {
            it = std::find_if(segments.begin(), segments.begin() + segments.size(),
                              [k](const auto &s) { return s.key > k; });
        } else {
            it = std::upper_bound(segments.begin(), segments.end(), k);
        }
        auto limit = it == segments.end() ? n : it->intercept;
        return std::min<size_t>((*std::prev(it))(k), limit - 1);
    }

    /** Execute a given function for each key in the sorted range [first, last). The function takes as the argument
     * an iterator to the current key and the corresponding approximate rank computed by the index. */
    template<typename ForwardIt, typename F>
    void for_each(ForwardIt first, ForwardIt last, F f, size_t n) const {
        auto segment_it = segments.begin();
        auto limit = segments.size() == 1 ? n : std::next(segment_it)->intercept;
        auto it = first;
        while (it != last) {
            if (std::next(segment_it) != segments.end() && *it >= std::next(segment_it)->key) {
                ++segment_it;
                limit = std::next(segment_it) == segments.end() ? n : std::next(segment_it)->intercept;
            }
            auto pos = std::min<size_t>((*segment_it)(*it), limit - 1);
            assert(std::abs(int64_t(pos) - int64_t(std::distance(first, it))) <= epsilon + 1);
            f(it, pos);
            ++it;
        }
    }

    [[nodiscard]] size_t epsilon_value() const { return epsilon; }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    [[nodiscard]] size_t segments_count() const { return segments.size(); }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    [[nodiscard]] size_t size_in_bytes() const { return sizeof(*this) + segments.size_in_bytes(); }

private:

    template<size_t ...Capacities>
    static std::pair<void *, size_t> make_helper(auto first, auto last, uint8_t epsilon,
                                                 std::integer_sequence<size_t, Capacities...> = {}) {
        std::vector<Segment> segments;
        auto n = (size_t) std::distance(first, last);
        segments.reserve(n / std::max(2, epsilon * epsilon));

        auto in_fun = [&](auto i) { return std::pair<K, size_t>(first[i], i); };
        auto out_fun = [&](auto cs, auto) { segments.emplace_back(cs); };
        internal::make_segmentation_mod(n, epsilon, in_fun, out_fun, false);

        size_t i = 0;
        void *ptr = nullptr;
        auto test = [&](auto threshold) {
            if (segments.size() <= threshold) {
                ptr = new ShortPGMIndex<threshold>(segments, epsilon);
                return false;
            }
            ++i;
            return true;
        };
        (test(std::integral_constant<std::size_t, Capacities>()) && ...)
            && (i = 0, ptr = new ShortPGMIndex<0>(segments, epsilon), true);
        return std::make_pair(ptr, i);
    }
};

#pragma pack(pop)

}