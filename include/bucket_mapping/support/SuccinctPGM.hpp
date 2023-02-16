#pragma once

#include <sdsl/int_vector.hpp>
#include "OptimalPiecewiseLinearModel.hpp"
#include "EliasFanoModified.hpp"

namespace pgm {

#pragma pack(push, 1)

uint8_t yshift_bit_width(uint8_t epsilon) { return sdsl::bits::hi(2 * epsilon) + 1; }

/** Store points using Elias-Fano representations with support for random access. */
struct EFPointsStorage {
    using xs_type = util::EliasFanoM;
    using ys_type = util::EliasFanoM;
    xs_type *xs;
    ys_type *ys;
    sdsl::int_vector<> yshifts;
    uint8_t epsilon;

    EFPointsStorage() = default;

    void load(const auto &xs_data, const auto &ys_data, const auto &yshifts_data, uint8_t epsilon) {
        this->epsilon = epsilon;
        xs = new xs_type(xs_data.size(), xs_data.back() + 1);
        ys = new ys_type(ys_data.size(), ys_data.back() + 1);
        for (auto x: xs_data)
            xs->push_back(x);
        for (auto y: ys_data)
            ys->push_back(y);
        xs->buildRankSelect();
        ys->buildRankSelect();
        yshifts.width(yshift_bit_width(epsilon));
        yshifts.resize(yshifts_data.size());
        std::copy(yshifts_data.begin(), yshifts_data.end(), yshifts.begin());
    }

    [[nodiscard]] size_t size_in_bytes() const { return xs->space() + ys->space() + sdsl::size_in_bytes(yshifts); }

    [[nodiscard]] uint8_t epsilon_value() const { return epsilon; }

    [[nodiscard]] size_t size() const { return *ys->at(ys->size() - 1); }

    [[nodiscard]] size_t segments_count() const { return xs->size() - 1; }

    [[nodiscard]] std::tuple<size_t, uint64_t, uint64_t, int64_t, int64_t, int64_t> segment_for_key(auto key) const {
        auto p = key >= xs->universe_size() - 1 ? xs->at(segments_count() - 1) : xs->predecessorPosition(key);
        auto np = p;
        ++np;
        auto [y0, y1, y2] = get_ys(p.index());
        return {p.index(), *p, *np, y0, y1, y2};
    }

    [[nodiscard]] std::tuple<uint64_t, uint64_t, int64_t, int64_t, int64_t> segment_at(size_t i) const {
        auto p = xs->at(i);
        auto np = p;
        ++np;
        auto [y0, y1, y2] = get_ys(i);
        return {*p, *np, y0, y1, y2};
    }

    ~EFPointsStorage() {
        delete xs;
        delete ys;
    }

    class Iterator {
        xs_type::ElementPointer xs_it;
        ys_type::ElementPointer ys_it;
        EFPointsStorage *storage;

    public:

        Iterator(xs_type::ElementPointer xs_it, ys_type::ElementPointer ys_it, EFPointsStorage *storage)
            : xs_it(xs_it), ys_it(ys_it), storage(storage) {}

        std::tuple<uint64_t, uint64_t, int64_t, int64_t, int64_t> operator*() {
            auto tmp1 = xs_it;
            auto x0 = *tmp1;
            auto x1 = *++tmp1;
            auto tmp2 = ys_it;
            auto y0 = *tmp2;
            auto y1 = *++tmp2;
            auto point_index = xs_it.index();
            auto shift1 = int64_t(storage->yshifts[point_index * 2]) - int64_t(storage->epsilon);
            auto shift2 = int64_t(storage->yshifts[point_index * 2 + 1]) - int64_t(storage->epsilon);
            auto shift3 = int64_t(storage->yshifts[point_index * 2 + 2]) - int64_t(storage->epsilon);
            return {x0, x1, y0 - shift1, y1 - shift2, y1 - shift3};
        }

        void operator++() {
            ++xs_it;
            ++ys_it;
        }
    };

    Iterator begin() { return {xs->begin(), ys->begin(), this}; }

private:

    [[nodiscard]] std::tuple<int64_t, int64_t, int64_t> get_ys(size_t point_index) const {
        auto p = ys->at(point_index);
        int64_t y0 = *p;
        int64_t y1 = *++p;
        auto shift1 = int64_t(yshifts[point_index * 2]) - int64_t(epsilon);
        auto shift2 = int64_t(yshifts[point_index * 2 + 1]) - int64_t(epsilon);
        auto shift3 = int64_t(yshifts[point_index * 2 + 2]) - int64_t(epsilon);
        return {y0 - shift1, y1 - shift2, y1 - shift3};
    }
};

template<typename PointsStorage = EFPointsStorage>
class SuccinctPGMIndex {
    using K = uint64_t;
    K first_key;
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
    SuccinctPGMIndex(RandomIt first, RandomIt last, uint8_t epsilon)
        : first_key(std::distance(first, last) ? *first : 0ull)  {
        auto n = (size_t) std::distance(first, last);
        if (n == 0)
            return;

        std::vector<K> xs_data = {0};
        std::vector<size_t> ys_data = {0};
        std::vector<uint32_t> yshifts_data;
        auto yshift_width = yshift_bit_width(epsilon);

        auto in_fun = [&](auto i) { return std::pair<K, size_t>(first[i], i); };
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
                || sdsl::bits::hi(yshift1) + 1 > yshift_width
                || sdsl::bits::hi(yshift2) + 1 > yshift_width)
                throw std::runtime_error("Unexpected shift > epsilon");

            xs_data.emplace_back(last_point.first - first_key);
            ys_data.emplace_back(last_point.second);
            yshifts_data.push_back(yshift1);
            yshifts_data.push_back(yshift2);
        };

        auto ignore_last = *std::prev(last) == std::numeric_limits<K>::max();
        internal::make_segmentation_mod(n - ignore_last, epsilon, in_fun, out_fun, true);

        ys_data.push_back(n);
        yshifts_data.push_back(epsilon);
        points.load(xs_data, ys_data, yshifts_data, epsilon);
        assert(size() == n);
    }

    /** Execute a given function for each key in the sorted range [first, last). The function takes as the argument
     * an iterator to the current key and the corresponding approximate rank computed by the index. */
    template<typename ForwardIt, typename F>
    void for_each(ForwardIt first, ForwardIt last, F f) {
        size_t segment_i = 0;
        auto segment_it = points.begin();
        auto [x0, x1, y0, y1, y2] = *segment_it;

        auto it = first;
        while (it != last) {
            assert(*it >= first_key);
            if (*it >= first_key + x1 && segment_i != segments_count() - 1) {
                ++segment_i;
                ++segment_it;
                std::tie(x0, x1, y0, y1, y2) = *segment_it;
            }
            auto eval = evaluate(*it - first_key, x0, x1, y0, y1);
            auto pos = std::min<size_t>(eval > 0 ? size_t(eval) : 0ull, y2 - 1);
            assert(std::abs(int64_t(pos) - int64_t(std::distance(first, it))) <= points.epsilon_value() + 1);
            f(it, pos);
            ++it;
        }
    }

    [[nodiscard]] size_t approximate_rank(K key) const {
        auto k = std::max(first_key, key);
        auto [point_index, x0, x1, y0, y1, y2] = points.segment_for_key(k - first_key);
        auto eval = evaluate(k - first_key, x0, x1, y0, y1);
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

    [[nodiscard]] size_t evaluate(K k, K x0, K x1, int64_t y0, int64_t y1) const {
        return y0 + ((__int128(k) - x0) * (y1 - y0)) / (x1 - x0);
    }
};

#pragma pack(pop)

}