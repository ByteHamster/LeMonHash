#pragma once

#include <sdsl/int_vector.hpp>
#include <compact_elias_fano.hpp>
#include <strict_elias_fano.hpp>
#include "OptimalPiecewiseLinearModel.hpp"
#include "support/sequence/EliasFanoModified.hpp"
#include "PointsStorage.hpp"

namespace lemonhash::pgm {

#pragma pack(push, 1)
template<typename PointsStorage = EFPointsStorage>
class SuccinctPGMIndex {
    using K = uint64_t;
    PointsStorage points;

public:

    /**
     * Constructs an empty index.
     */
    SuccinctPGMIndex() = default;

    /**
     * Constructs the index on the sorted keys in the range [first, last).
     * @param first, last the range containing the sorted keys to be indexed
     */
    template<typename RandomIt>
    SuccinctPGMIndex(RandomIt first, RandomIt last, uint8_t epsilon) {
        auto n = (size_t) std::distance(first, last);
        if (n == 0)
            return;

        std::vector<K> xs_data = {0};
        std::vector<size_t> ys_data = {0};
        std::vector<uint32_t> yshifts_data;
        auto yshift_width = yshift_bit_width(epsilon);

        auto expected_segments = n / std::max(2, epsilon * epsilon);
        xs_data.reserve(expected_segments + 1);
        ys_data.reserve(expected_segments + 2);
        yshifts_data.reserve(expected_segments * 2 + 1);

        bool skip_first = n > 2; // the first key will always be mapped to rank 0, so we skip it in the segmentation
        auto first_key = first[skip_first];

        auto in_fun = [&](auto i) { return std::pair<K, size_t>(first[i + skip_first], i + skip_first); };
        auto out_fun = [&](auto cs, auto last_point) {
            auto x0 = cs.rectangle[1].x;
            auto y0 = cs.rectangle[1].y;
            auto x1 = cs.rectangle[3].x;
            auto y1 = cs.rectangle[3].y;
            auto eval1 = evaluate<false>(cs.first.x, x0, x1, y0, y1);
            auto eval2 = evaluate<false>(last_point.first, x0, x1, y0, y1);
            auto yshift1 = cs.first.y - eval1 + epsilon;
            auto yshift2 = last_point.second - eval2 + epsilon;
            if (std::bit_width(yshift1) > yshift_width || std::bit_width(yshift2) > yshift_width)
                throw std::runtime_error("Unexpected shift > epsilon");

            xs_data.emplace_back(last_point.first - first_key);
            ys_data.emplace_back(last_point.second);
            yshifts_data.push_back(yshift1);
            yshifts_data.push_back(yshift2);
        };
        internal::make_segmentation_mod(n - skip_first, epsilon, in_fun, out_fun, true);

        ys_data.push_back(n);
        yshifts_data.push_back(epsilon);
        points.load(first_key, xs_data, ys_data, yshifts_data, epsilon);
        assert(size() == n);
    }

    /** Execute a given function for each key in the sorted range [first, last). The function takes as the argument
     * an iterator to the current key and the corresponding approximate rank computed by the index. */
    template<typename ForwardIt, typename F>
    void for_each(ForwardIt first, ForwardIt last, F f) const {
        size_t segment_i = 0;
        auto segment_it = points.begin();
        auto [x0, x1, y0, y1, y2] = *segment_it;

        auto it = first;
        auto first_key = points.first_key();
        while (it != last && *it < first_key) {
            f(it, 0);
            ++it;
        }

        while (it != last) {
            assert(*it >= first_key);
            if (*it >= first_key + x1 && segment_i != segments_count() - 1) {
                ++segment_i;
                ++segment_it;
                std::tie(x0, x1, y0, y1, y2) = *segment_it;
            }
            auto eval = evaluate<true>(*it - first_key, x0, x1, y0, y1);
            auto pos = std::min<size_t>(eval > 0 ? size_t(eval) : 0ull, y2 - 1);
            assert(std::abs(int64_t(pos) - int64_t(std::distance(first, it))) <= points.epsilon_value() + 1);
            f(it, pos);
            ++it;
        }
    }

    [[nodiscard]] size_t approximate_rank(K key) const {
        auto first_key = points.first_key();
        if (key < first_key)
            return 0;
        auto [point_index, x0, x1, y0, y1, y2] = points.segment_for_key(key - first_key);
        auto eval = evaluate<true>(key - first_key, x0, x1, y0, y1);
        auto pos = std::min<size_t>(eval > 0 ? size_t(eval) : 0ull, y2 - 1);
        return pos;
    }

    /**
     * Returns the number of segments in the index.
     * @return the number of segments
     */
    [[nodiscard]] size_t segments_count() const { return points.segments_count(); }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    [[nodiscard]] size_t size_in_bytes() const { return sizeof(*this) + points.size_in_bytes(); }

    [[nodiscard]] size_t epsilon_value() const { return points.epsilon_value(); }

    /**
     * Returns the number of elements the index was built on.
     * @return the number of elements the index was built on
     */
    [[nodiscard]] size_t size() const { return points.size(); }

private:

    template<bool Query>
    [[nodiscard]] static inline size_t evaluate(K k, K x0, K x1, size_t y0, size_t y1) {
        assert(x0 <= x1 && y0 <= y1);
        if constexpr (Query) {
            assert(k >= x0);
            size_t mul;
            if (!__builtin_mul_overflow(k - x0, y1 - y0, &mul))
                return y0 + mul / (x1 - x0);
        } else {
            int64_t dk;
            int64_t mul;
            if (!__builtin_sub_overflow(k, x0, &dk) && !__builtin_mul_overflow(dk, int64_t(y1 - y0), &mul))
                return int64_t(y0) + mul / int64_t(x1 - x0);
        }
        return y0 + ((__int128(k) - x0) * (y1 - y0)) / (x1 - x0);
    }
};

#pragma pack(pop)

}