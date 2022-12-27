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
    static constexpr uint8_t type_bit_width = BIT_WIDTH(max_local_capacity);
    static constexpr auto align_val = std::align_val_t(std::max<size_t>(__STDCPP_DEFAULT_NEW_ALIGNMENT__,
                                                                        1ull << type_bit_width));

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
        Heap(const Heap &) = delete;
        Heap(Heap&&) = delete;
        Heap& operator=(const Heap &) = delete;
        Heap& operator=(Heap &&) = delete;
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
    uint32_t n;

public:

    static constexpr size_t capacity = Capacity;

    using capacities_sequence = std::make_index_sequence<max_local_capacity + 1>;
    using all_types = decltype(to_tuple_helper(capacities_sequence{}));

    class Holder {
        uint64_t ptr : 64 - type_bit_width;
        size_t type : type_bit_width;

        template <typename F>
        void visit(F f) const {
            auto p = ptr << type_bit_width;
            switch (type) {
                case 0:  f(reinterpret_cast<std::tuple_element_t<0,  all_types>*>(p)); break;
                case 1:  f(reinterpret_cast<std::tuple_element_t<1,  all_types>*>(p)); break;
                case 2:  f(reinterpret_cast<std::tuple_element_t<2,  all_types>*>(p)); break;
                case 3:  f(reinterpret_cast<std::tuple_element_t<3,  all_types>*>(p)); break;
                case 4:  f(reinterpret_cast<std::tuple_element_t<4,  all_types>*>(p)); break;
                case 5:  f(reinterpret_cast<std::tuple_element_t<5,  all_types>*>(p)); break;
                case 6:  f(reinterpret_cast<std::tuple_element_t<6,  all_types>*>(p)); break;
                case 7:  f(reinterpret_cast<std::tuple_element_t<7,  all_types>*>(p)); break;
                case 8:  f(reinterpret_cast<std::tuple_element_t<8,  all_types>*>(p)); break;
                case 9:  f(reinterpret_cast<std::tuple_element_t<9,  all_types>*>(p)); break;
                case 10: f(reinterpret_cast<std::tuple_element_t<10, all_types>*>(p)); break;
                case 11: f(reinterpret_cast<std::tuple_element_t<11, all_types>*>(p)); break;
                case 12: f(reinterpret_cast<std::tuple_element_t<12, all_types>*>(p)); break;
                case 13: f(reinterpret_cast<std::tuple_element_t<13, all_types>*>(p)); break;
                case 14: f(reinterpret_cast<std::tuple_element_t<14, all_types>*>(p)); break;
                case 15: f(reinterpret_cast<std::tuple_element_t<15, all_types>*>(p)); break;
                default: f(reinterpret_cast<ShortPGMIndex<>*>(p));
            }
        }

        template <typename F>
        void visit(F f) { const_cast<const Holder*>(this)->visit(f); }

    public:

        Holder() : ptr(0), type(0) {}
        Holder(void *ptr, size_t type) : ptr(reinterpret_cast<uint64_t>(ptr) >> type_bit_width), type(type) {}
        Holder(const Holder&) = delete;
        Holder& operator=(const Holder &) = delete;
        Holder(Holder&&) = delete;

        Holder& operator=(Holder &&other) {
            if (this != &other) {
                if (ptr != 0 && type != 0)
                    visit([](auto index) { operator delete (index, align_val); });

                ptr = other.ptr;
                type = other.type;

                other.ptr = 0;
                other.type = 0;
            }
            return *this;
        };

        [[nodiscard]] size_t approximate_rank(const K &key) const {
            size_t rank = 0;
            visit([&](const auto index) { rank = index->approximate_rank(key); });
            return rank;
        }

        template<typename ForwardIt, typename F>
        void for_each(ForwardIt first, ForwardIt last, F f) const {
            visit([&](const auto index) { index->for_each(first, last, f); });
        }

        [[nodiscard]] size_t size() const {
            size_t size = 0;
            visit([&](const auto index) { size = index->size(); });
            return size;
        }

        [[nodiscard]] size_t size_in_bytes() const {
            size_t bytes = sizeof(*this);
            visit([&](const auto index) { bytes += index->size_in_bytes(); });
            return bytes;
        }

        ~Holder() {
            if (ptr != 0 && type != 0) {
                visit([](auto index) { operator delete (index, align_val); });
                ptr = 0;
                type = 0;
            }
        }
    };

    static Holder make(auto first, auto last, uint8_t epsilon) {
        return make_helper(first, last, epsilon, capacities_sequence{});
    }

    ShortPGMIndex() = default;

    ShortPGMIndex(const std::vector<Segment> &s, const size_t n) : segments(s.size()), n(n) {
        assert(s.size() <= segments.capacity());
        std::copy(s.begin(), s.end(), &segments.segments[0]);
    }

    /** Returns the approximate rank of @p key. */
    [[nodiscard]] size_t approximate_rank(const K &key) const {
        auto k = std::max(segments.begin()->key, key);
        decltype(segments.begin()) it;
        if constexpr (Capacity > 0 && Capacity <= 7) {
            it = std::find_if(segments.begin(), segments.begin() + segments.size(),
                              [k](const auto &s) { return s.key > k; });
        } else {
            it = std::upper_bound(segments.begin(), segments.end(), k);
        }
        auto limit = it == segments.end() ? size() : it->intercept;
        return std::min<size_t>((*std::prev(it))(k), limit - 1);
    }

    /** Execute a given function for each key in the sorted range [first, last). The function takes as the argument
     * an iterator to the current key and the corresponding approximate rank computed by the index. */
    template<typename ForwardIt, typename F>
    void for_each(ForwardIt first, ForwardIt last, F f) const {
        auto segment_it = segments.begin();
        auto limit = segments.size() == 1 ? n : std::next(segment_it)->intercept;
        auto it = first;
        while (it != last) {
            if (std::next(segment_it) != segments.end() && *it >= std::next(segment_it)->key) {
                ++segment_it;
                limit = std::next(segment_it) == segments.end() ? n : std::next(segment_it)->intercept;
            }
            auto pos = std::min<size_t>((*segment_it)(*it), limit - 1);
            //assert(std::abs(int64_t(pos) - int64_t(std::distance(first, it))) <= epsilon + 1);
            f(it, pos);
            ++it;
        }
    }

    /** Returns the number of elements the index was built on. */
    [[nodiscard]] size_t size() const { return n; }

    /**  Returns the size of the index in bytes. */
    [[nodiscard]] size_t size_in_bytes() const { return sizeof(*this) + segments.size_in_bytes(); }

private:

    template<size_t ...Capacities>
    static Holder make_helper(auto first, auto last, uint8_t epsilon,
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
                ptr = new(align_val) ShortPGMIndex<threshold>(segments, n);
                return false;
            }
            ++i;
            return true;
        };
        (test(std::integral_constant<std::size_t, Capacities>()) && ...)
            && (i = 0, ptr = new(align_val) ShortPGMIndex<0>(segments, n), true);

        assert((reinterpret_cast<uint64_t>(ptr) & (size_t(align_val) - 1)) == 0);
        return {ptr, i};
    }
};

#pragma pack(pop)

}