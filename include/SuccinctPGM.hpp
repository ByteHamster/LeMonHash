#pragma once

#include <sdsl/int_vector.hpp>
#include "OptimalPiecewiseLinearModel.hpp"
#include "EliasFanoModified.hpp"

namespace pgm {

template<typename K, size_t Epsilon = 64>
class SuccinctPGMIndex {
protected:
    static_assert(Epsilon > 0);

    using xs_type = util::EliasFanoM;
    using ys_type = util::EliasFanoM;

    size_t n;
    K first_key;
    xs_type *xs;
    ys_type *ys;
    sdsl::int_vector<> yshifts;

public:

    static constexpr size_t epsilon_value = Epsilon;

    /**
     * Constructs an empty index.
     */
    SuccinctPGMIndex() = default;

    /**
     * Constructs the index on the given sorted vector.
     * @param data the vector of keys to be indexed, must be sorted
     */
    explicit SuccinctPGMIndex(const std::vector<K> &data) : SuccinctPGMIndex(data.begin(), data.end()) {}

    /**
     * Constructs the index on the sorted keys in the range [first, last).
     * @param first, last the range containing the sorted keys to be indexed
     */
    template<typename RandomIt>
    SuccinctPGMIndex(RandomIt first, RandomIt last)
        : n(std::distance(first, last)),
          first_key(n ? *first : K(0))  {
        if (n == 0)
            return;

        using namespace internal2;

        auto ignore_last = *std::prev(last) == std::numeric_limits<K>::max(); // max() is the sentinel value
        auto last_n = n - ignore_last;
        last -= ignore_last;

        std::vector<K> xs_data = {0};
        std::vector<size_t> ys_data = {0};
        std::vector<size_t> yshifts_data;
        yshifts.width(sdsl::bits::hi(2 * epsilon_value) + 1);

        auto in_fun = [&](auto i) {
            auto x = first[i];
            // Here there is an adjustment for inputs with duplicate keys: at the end of a run of duplicate keys equal
            // to x=first[i] such that x+1!=first[i+1], we map the values x+1,...,first[i+1]-1 to their correct rank i
            auto flag = i > 0 && i + 1u < n && x == first[i - 1] && x != first[i + 1] && x + 1 != first[i + 1];
            return std::pair<K, size_t>(x + flag, i);
        };

        typename OptimalPiecewiseLinearModel<K, size_t>::CanonicalSegment prev_cs;
        auto out_fun = [&](auto cs, auto last_point) {
            auto x0 = cs.rectangle[1].x;
            auto y0 = cs.rectangle[1].y;
            auto x1 = cs.rectangle[3].x;
            auto y1 = cs.rectangle[3].y;
            auto eval1 = y0 + ((__int128(cs.first.x) - x0) * (y1 - y0)) / (x1 - x0);
            auto eval2 = y0 + ((__int128(last_point.first) - x0) * (y1 - y0)) / (x1 - x0);
            auto yshift1 = cs.first.y - eval1 + epsilon_value;
            auto yshift2 = last_point.second - eval2 + epsilon_value;
            if (yshift1 < 0 || yshift2 < 0
                || sdsl::bits::hi(yshift1) + 1 > yshifts.width()
                || sdsl::bits::hi(yshift2) + 1 > yshifts.width())
                throw std::runtime_error("Unexpected shift > epsilon");

            xs_data.emplace_back(last_point.first - first_key);
            ys_data.emplace_back(last_point.second);
            yshifts_data.push_back(yshift1);
            yshifts_data.push_back(yshift2);
        };

        make_segmentation_mod(last_n, epsilon_value, in_fun, out_fun, true);

        ys_data.push_back(n);
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

    [[nodiscard]] size_t approximate_rank(K key) const {
        auto k = std::max(first_key, key);
        auto [point_index, x0, x1] = segment_for_key(k - first_key);
        auto [y0, y1, y2] = get_ys(point_index);
        auto eval = y0 + (__int128(k - (x0 + first_key)) * (y1 - y0)) / (x1 - x0);
        auto pos = std::min<size_t>(eval > 0 ? size_t(eval) : 0ull, y2);
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
        return xs->space() + ys->space() + sdsl::size_in_bytes(yshifts);
    }

    ~SuccinctPGMIndex() {
        delete xs;
        delete ys;
    }

private:

    [[nodiscard]] std::tuple<int64_t, int64_t, int64_t> get_ys(size_t point_index) const {
        auto shift1 = int64_t(yshifts[point_index * 2]) - int64_t(epsilon_value);
        auto shift2 = int64_t(yshifts[point_index * 2 + 1]) - int64_t(epsilon_value);

        auto p = ys->at(point_index);
        int64_t y0 = *p;
        ++p;
        int64_t y1 = *p;
        ++p;
        int64_t y2 = *p;

        y0 -= shift1;
        y1 -= shift2;
        return {y0, y1, y2};
    }

    [[nodiscard]] std::tuple<size_t, K, K> segment_for_key(K key) const {
        auto p = key >= xs->universe_size() - 1 ? xs->at(xs->size() - 2) : xs->predecessorPosition(key);
        auto np = p;
        ++np;
        return {p.index(), *p, *np};
    }

};

}