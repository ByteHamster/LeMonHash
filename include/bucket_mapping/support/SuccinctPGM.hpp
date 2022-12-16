#pragma once

#include <sdsl/int_vector.hpp>
#include "OptimalPiecewiseLinearModel.hpp"
#include "EliasFanoModified.hpp"

namespace pgm {

#pragma pack(push, 1)

template<typename K>
class SuccinctPGMIndex {
protected:
    using xs_type = util::EliasFanoM;
    using ys_type = util::EliasFanoM;

    uint16_t epsilon;
    K first_key;
    xs_type *xs;
    ys_type *ys;
    sdsl::int_vector<> yshifts;

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
    SuccinctPGMIndex(RandomIt first, RandomIt last, uint16_t epsilon)
        : epsilon(epsilon),
          first_key(std::distance(first, last) ? *first : K(0))  {
        auto n = (size_t) std::distance(first, last);
        if (n == 0)
            return;

        using namespace internal;

        auto ignore_last = *std::prev(last) == std::numeric_limits<K>::max(); // max() is the sentinel value
        auto last_n = n - ignore_last;
        last -= ignore_last;

        std::vector<K> xs_data = {0};
        std::vector<size_t> ys_data = {0};
        std::vector<size_t> yshifts_data;
        yshifts.width(sdsl::bits::hi(2 * epsilon) + 1);

        auto in_fun = [&](auto i) { return std::pair<K, size_t>(first[i], i); };

        typename OptimalPiecewiseLinearModel<K, size_t>::CanonicalSegment prev_cs;
        auto out_fun = [&](auto cs, auto last_point) {
            auto x0 = cs.rectangle[1].x;
            auto y0 = cs.rectangle[1].y;
            auto x1 = cs.rectangle[3].x;
            auto y1 = cs.rectangle[3].y;
            auto eval1 = evaluate(cs.first.x, x0, x1, y0, y1);
            auto eval2 = evaluate(last_point.first, x0, x1, y0, y1);
            auto yshift1 = cs.first.y - eval1 + epsilon;
            auto yshift2 = last_point.second - eval2 + epsilon;
            if (yshift1 < 0 || yshift2 < 0
                || sdsl::bits::hi(yshift1) + 1 > yshifts.width()
                || sdsl::bits::hi(yshift2) + 1 > yshifts.width())
                throw std::runtime_error("Unexpected shift > epsilon");

            xs_data.emplace_back(last_point.first - first_key);
            ys_data.emplace_back(last_point.second);
            yshifts_data.push_back(yshift1);
            yshifts_data.push_back(yshift2);
        };

        make_segmentation_mod(last_n, epsilon, in_fun, out_fun, true);

        ys_data.push_back(n);
        yshifts_data.push_back(epsilon);
        xs = new xs_type(xs_data.size(), xs_data.back() + 1);
        ys = new ys_type(ys_data.size(), ys_data.back() + 1);
        for (auto x: xs_data)
            xs->push_back(x);
        for (auto y: ys_data)
            ys->push_back(y);
        xs->buildRankSelect();
        ys->buildRankSelect();

        yshifts.resize(yshifts_data.size());
        std::copy(yshifts_data.begin(), yshifts_data.end(), yshifts.begin());
    }

    /** Execute a given function for each key in the sorted range [first, last). The function takes as the argument
     * an iterator to the current key and the corresponding approximate rank computed by the index. */
    template<typename ForwardIt, typename F>
    void for_each(ForwardIt first, ForwardIt last, F f) {
        size_t segment_i = 0;
        auto [x0, x1] = segment_at(0);
        auto [y0, y1, y2] = get_ys(0);
        auto it = first;
        while (it != last) {
            assert(*it >= first_key);
            if (*it >= first_key + x1 && segment_i != segments_count() - 2) {
                ++segment_i;
                std::tie(x0, x1) = segment_at(segment_i);
                std::tie(y0, y1, y2) = get_ys(segment_i);
            }
            auto eval = evaluate(*it - first_key, x0, x1, y0, y1);
            auto pos = std::min<size_t>(eval > 0 ? size_t(eval) : 0ull, y2 - 1);
            assert(std::abs(int64_t(pos) - int64_t(std::distance(first, it))) <= epsilon + 1);
            f(it, pos);
            ++it;
        }
    }

    [[nodiscard]] size_t approximate_rank(K key) const {
        auto k = std::max(first_key, key);
        auto [point_index, x0, x1] = segment_for_key(k - first_key);
        auto [y0, y1, y2] = get_ys(point_index);
        auto eval = evaluate(k - first_key, x0, x1, y0, y1);
        auto pos = std::min<size_t>(eval > 0 ? size_t(eval) : 0ull, y2 - 1);
        return pos;
    }

    /**
     * Returns the number of segments in the index.
     * @return the number of segments
     */
    [[nodiscard]] size_t segments_count() const { return xs->size(); }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    [[nodiscard]] size_t size_in_bytes() const {
        return sizeof(*this) + xs->space() + ys->space() + sdsl::size_in_bytes(yshifts);
    }

    [[nodiscard]] size_t epsilon_value() const { return epsilon; }

    ~SuccinctPGMIndex() {
        delete xs;
        delete ys;
    }

private:

    [[nodiscard]] size_t evaluate(K k, K x0, K x1, int64_t y0, int64_t y1) const {
        return y0 + ((__int128(k) - x0) * (y1 - y0)) / (x1 - x0);
    }

    [[nodiscard]] std::tuple<int64_t, int64_t, int64_t> get_ys(size_t point_index) const {
        auto p = ys->at(point_index);
        int64_t y0 = *p;
        int64_t y1 = *++p;

        auto shift1 = int64_t(yshifts[point_index * 2]) - int64_t(epsilon);
        auto shift2 = int64_t(yshifts[point_index * 2 + 1]) - int64_t(epsilon);
        auto shift3 = int64_t(yshifts[point_index * 2 + 2]) - int64_t(epsilon);

        return {y0 - shift1, y1 - shift2, y1 - shift3};
    }

    [[nodiscard]] std::tuple<size_t, K, K> segment_for_key(K key) const {
        auto p = key >= xs->universe_size() - 1 ? xs->at(segments_count() - 2) : xs->predecessorPosition(key);
        auto np = p;
        ++np;
        return {p.index(), *p, *np};
    }

    [[nodiscard]] std::pair<K, K> segment_at(size_t i) const {
        auto p = xs->at(i);
        auto np = p;
        ++np;
        return {*p, *np};
    }
};

#pragma pack(pop)

}