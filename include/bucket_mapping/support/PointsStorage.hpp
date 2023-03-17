#pragma once

#include <sdsl/int_vector.hpp>
#include "support/sequence/EliasFanoModified.hpp"
#include <compact_elias_fano.hpp>
#include <strict_elias_fano.hpp>
#include <partitioned_sequence.hpp>
#include <uniform_partitioned_sequence.hpp>

namespace pgm {

#pragma pack(push, 1)

uint8_t yshift_bit_width(uint8_t epsilon) { return sdsl::bits::hi(2 * epsilon) + 1; }

/** Store points using Elias-Fano representations with support for random access. */
struct EFPointsStorage {
    using xs_type = lemonhash::EliasFanoM;
    using ys_type = lemonhash::EliasFanoM;
    xs_type *xs;
    ys_type *ys;
    sdsl::int_vector<> yshifts;
    uint8_t epsilon;
    uint64_t first;

    EFPointsStorage() = default;

    EFPointsStorage(const EFPointsStorage &) = delete;

    EFPointsStorage(EFPointsStorage &&other) {
        xs = other.xs;
        ys = other.ys;
        yshifts = std::move(other.yshifts);
        epsilon = other.epsilon;
        first = other.first;
        other.xs = nullptr;
        other.ys = nullptr;
    }

    EFPointsStorage &operator=(const EFPointsStorage &) = delete;

    EFPointsStorage &operator=(EFPointsStorage &&other) {
        if (this != &other) {
            xs = other.xs;
            ys = other.ys;
            yshifts = std::move(other.yshifts);
            epsilon = other.epsilon;
            first = other.first;
            other.xs = nullptr;
            other.ys = nullptr;
        }
        return *this;
    }

    void load(uint64_t first, const auto &xs_data, const auto &ys_data, const auto &yshifts_data, uint8_t epsilon) {
        this->first = first;
        this->epsilon = epsilon;
        xs = new xs_type(xs_data.size(), xs_data.back() + 1);
        ys = new ys_type(ys_data.size(), ys_data.back() - ys_data.size() + 2);
        for (auto x: xs_data)
            xs->push_back(x);
        size_t i = 0;
        for (auto y: ys_data)
            ys->push_back(y - i++);
        xs->buildRankSelect();
        ys->buildRankSelect();
        yshifts.width(yshift_bit_width(epsilon));
        yshifts.resize(yshifts_data.size());
        std::copy(yshifts_data.begin(), yshifts_data.end(), yshifts.begin());
    }

    [[nodiscard]] size_t size_in_bytes() const { return xs->space() + ys->space() + sdsl::size_in_bytes(yshifts); }

    [[nodiscard]] uint8_t epsilon_value() const { return epsilon; }

    [[nodiscard]] size_t size() const { return *ys->at(ys->size() - 1) + ys->size() - 1; }

    [[nodiscard]] size_t segments_count() const { return xs->size() - 1; }

    [[nodiscard]] std::tuple<size_t, uint64_t, uint64_t, int64_t, int64_t, int64_t> segment_for_key(auto key) const {
        auto p = key >= xs->universe_size() - 1 ? xs->at(segments_count() - 1) : xs->predecessorPosition(key);
        auto np = p;
        ++np;
        auto [y0, y1, y2] = get_ys(p.index());
        return {p.index(), *p, *np, y0, y1, y2};
    }

    [[nodiscard]] uint64_t first_key() const { return first; }

    ~EFPointsStorage() {
        delete xs;
        delete ys;
    }

    class Iterator {
        xs_type::ElementPointer xs_it;
        ys_type::ElementPointer ys_it;
        const EFPointsStorage *storage;

    public:

        Iterator(const xs_type::ElementPointer xs_it, const ys_type::ElementPointer ys_it, const EFPointsStorage *storage)
            : xs_it(xs_it), ys_it(ys_it), storage(storage) {}

        std::tuple<uint64_t, uint64_t, int64_t, int64_t, int64_t> operator*() {
            auto tmp1 = xs_it;
            auto x0 = *tmp1;
            auto x1 = *++tmp1;
            auto tmp2 = ys_it;
            auto y0 = *tmp2 + tmp2.index();
            ++tmp2;
            auto y1 = *tmp2 + tmp2.index();
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

    Iterator begin() const { return {xs->begin(), ys->begin(), this}; }

private:

    [[nodiscard]] std::tuple<int64_t, int64_t, int64_t> get_ys(size_t point_index) const {
        auto p = ys->at(point_index);
        int64_t y0 = *p + point_index;
        int64_t y1 = *++p + point_index + 1;
        auto shift1 = int64_t(yshifts[point_index * 2]) - int64_t(epsilon);
        auto shift2 = int64_t(yshifts[point_index * 2 + 1]) - int64_t(epsilon);
        auto shift3 = int64_t(yshifts[point_index * 2 + 2]) - int64_t(epsilon);
        return {y0 - shift1, y1 - shift2, y1 - shift3};
    }
};

/** Store points using Elias-Fano representations with support for random access. */
class EFPointsStorageV2 {
    using xs_type = ds2i::compact_elias_fano; //ds2i::partitioned_sequence<ds2i::indexed_sequence>;
    using ys_type = ds2i::strict_elias_fano;
    succinct::bit_vector bv;

    static ds2i::global_parameters xs_params;
    static ds2i::global_parameters ys_params;

    // data layout:
    // first_key in key_bit_width bits
    // xs's universe size in key_bit_width bits
    // ys's universe size in size_bit_width bits
    // epsilon in epsilon_bit_width bits
    // nsegments in nsegments_bit_width bits
    // (nsegments*2+1) yshifts, each in yshift_bit_width(epsilon) bits
    // ys in Elias-Fano representation
    // xs in Elias-Fano representation

    static constexpr uint8_t key_bit_width = 64;
    static constexpr uint8_t size_bit_width = 32;
    static constexpr uint8_t epsilon_bit_width = 8;
    static constexpr uint8_t nsegments_bit_width = 24;

public:

    EFPointsStorageV2() {
        [[maybe_unused]] static bool execute_once = [](){
            ys_params.ef_log_sampling0 = 63; // rank is not needed on the ys
            return true;
        } ();
    }

    EFPointsStorageV2(const EFPointsStorageV2 &) = delete;

    EFPointsStorageV2(EFPointsStorageV2 &&other) {
        bv.swap(other.bv);
    }

    EFPointsStorageV2 &operator=(const EFPointsStorageV2 &) = delete;

    EFPointsStorageV2 &operator=(EFPointsStorageV2 &&other) {
        if (this != &other)
            bv.swap(other.bv);
        return *this;
    }

    void load(uint64_t first, const auto &xs_data, const auto &ys_data, const auto &yshifts_data, uint8_t epsilon) {
        succinct::bit_vector_builder xs_bvb;
        xs_type::write(xs_bvb, xs_data.begin(), xs_data.back() + 1, xs_data.size(), xs_params);

        auto xs_bit_size = xs_bvb.size();
        auto ys_bit_size = ys_type::bitsize(ys_params, ys_data.back() + 1, ys_data.size());

        auto yshift_width = yshift_bit_width(epsilon);
        auto total_bit_size = key_bit_width + key_bit_width + size_bit_width + epsilon_bit_width + nsegments_bit_width
            + yshifts_data.size() * yshift_width + ys_bit_size + xs_bit_size;

        succinct::bit_vector_builder bvb;
        bvb.reserve(total_bit_size);
        auto nsegments = xs_data.size() - 1;

        bvb.append_bits(first, key_bit_width);
        bvb.append_bits(xs_data.back() + 1, key_bit_width);
        bvb.append_bits(ys_data.back() + 1, size_bit_width);
        bvb.append_bits(epsilon, epsilon_bit_width);
        bvb.append_bits(nsegments, nsegments_bit_width);

        for (size_t i = 0; i < yshifts_data.size(); ++i)
            bvb.append_bits(yshifts_data[i], yshift_width);

        ys_type::write(bvb, ys_data.begin(), ys_data.back() + 1, ys_data.size(), ys_params);
        bvb.append(xs_bvb);

        succinct::bit_vector(&bvb).swap(bv);

        assert(first_key() == first);
        assert(epsilon_value() == epsilon);
        assert(size() == ys_data.back());
        assert(segments_count() == nsegments);
        assert(xs_size() == xs_data.size());
        assert(ys_size() == ys_data.size());
        assert(xs_universe() == xs_data.back() + 1);
        assert(ys_universe() == ys_data.back() + 1);
    }

    [[nodiscard]] std::tuple<size_t, uint64_t, uint64_t, int64_t, int64_t, int64_t> segment_for_key(auto key) const {
        auto xs_it = xs_enumerator();
        auto next_geq = xs_it.next_geq(key);
        uint64_t x0, x1, pos;
        if (next_geq.first >= xs_size() - 1) {
            std::tie(pos, x0) = xs_it.move(xs_size() - 2);
            x1 = xs_it.next().second;
        } else if (next_geq.second != key && next_geq.first > 0) {
            x0 = xs_it.prev_value();
            x1 = next_geq.second;
            pos = next_geq.first - 1;
        } else {
            x0 = next_geq.second;
            x1 = xs_it.next().second;
            pos = next_geq.first;
        }

        auto [y0, y1, y2] = get_ys(pos);
        return {pos, x0, x1, y0, y1, y2};
    }

    [[nodiscard]] size_t size() const {
        auto bit_offset = key_bit_width + key_bit_width;
        return bv.get_bits(bit_offset, size_bit_width) - 1;
    }

    [[nodiscard]] size_t size_in_bytes() const { return bv.size() / 8; }

    [[nodiscard]] uint8_t epsilon_value() const {
        auto bit_offset = key_bit_width + key_bit_width + size_bit_width;
        return bv.get_bits(bit_offset, epsilon_bit_width);
    }

    [[nodiscard]] size_t segments_count() const {
        auto bit_offset = key_bit_width + key_bit_width + size_bit_width + epsilon_bit_width;
        return bv.get_bits(bit_offset, nsegments_bit_width);
    }

    [[nodiscard]] uint64_t first_key() const {
        return bv.get_bits(0, key_bit_width);
    }

    class Iterator {
        xs_type::enumerator xs_it;
        ys_type::enumerator ys_it;
        size_t point_index;
        uint64_t x_value;
        uint64_t y_value;
        const EFPointsStorageV2 *storage;

    public:

        Iterator(xs_type::enumerator xs_it, ys_type::enumerator ys_it, const EFPointsStorageV2 *storage)
            : xs_it(xs_it), ys_it(ys_it), storage(storage) {
            point_index = 0;
            x_value = this->xs_it.move(0).second;
            y_value = this->ys_it.move(0).second;
        }

        std::tuple<uint64_t, uint64_t, int64_t, int64_t, int64_t> operator*() {
            auto tmp1 = xs_it;
            auto tmp2 = ys_it;
            auto x0 = x_value;
            auto x1 = tmp1.next().second;
            auto y0 = y_value;
            auto y1 = tmp2.next().second;
            int64_t epsilon = storage->epsilon_value();
            auto shift1 = int64_t(storage->get_yshift(point_index * 2)) - epsilon;
            auto shift2 = int64_t(storage->get_yshift(point_index * 2 + 1)) - epsilon;
            auto shift3 = int64_t(storage->get_yshift(point_index * 2 + 2)) - epsilon;
            return {x0, x1, y0 - shift1, y1 - shift2, y1 - shift3};
        }

        void operator++() {
            ++point_index;
            x_value = xs_it.next().second;
            y_value = ys_it.next().second;
        }
    };

    Iterator begin() const { return {xs_enumerator(), ys_enumerator(), this}; }

private:

    xs_type::enumerator xs_enumerator() const { return {bv, xs_offset(), xs_universe(), xs_size(), xs_params}; }

    ys_type::enumerator ys_enumerator() const { return {bv, ys_offset(), ys_universe(), ys_size(), ys_params}; }

    uint64_t xs_universe() const { return bv.get_bits(key_bit_width, key_bit_width); }

    uint64_t ys_universe() const { return size() + 1; }

    size_t xs_size() const { return segments_count() + 1; }

    size_t ys_size() const { return segments_count() + 2; }

    size_t yshifts_size() const { return segments_count() * 2 + 1; }

    size_t yshifts_offset() const {
        return key_bit_width + key_bit_width + size_bit_width + epsilon_bit_width + nsegments_bit_width;
    }

    size_t ys_offset() const { return yshifts_offset() + yshifts_size() * yshift_bit_width(epsilon_value()); }

    size_t xs_offset() const { return ys_offset() + ys_type::bitsize(ys_params, ys_universe(), ys_size()); }

    uint64_t get_yshift(size_t index) const {
        auto bit_offset = yshifts_offset() + index * yshift_bit_width(epsilon_value());
        return bv.get_bits(bit_offset, yshift_bit_width(epsilon_value()));
    }

    [[nodiscard]] std::tuple<int64_t, int64_t, int64_t> get_ys(size_t point_index) const {
        auto it = ys_enumerator();
        int64_t y0 = it.move(point_index).second;
        int64_t y1 = it.next().second;
        int64_t epsilon = epsilon_value();
        auto shift1 = int64_t(get_yshift(point_index * 2)) - epsilon;
        auto shift2 = int64_t(get_yshift(point_index * 2 + 1)) - epsilon;
        auto shift3 = int64_t(get_yshift(point_index * 2 + 2)) - epsilon;
        return {y0 - shift1, y1 - shift2, y1 - shift3};
    }
};

ds2i::global_parameters EFPointsStorageV2::xs_params = {};
ds2i::global_parameters EFPointsStorageV2::ys_params = {};

/** Store points using Elias-Fano representations without support for efficient random access. */
class EFSequentialPointsStorage {
    uint64_t *data;
    // data layout:
    // first_key in key_bit_width bits
    // epsilon in epsilon_bit_width bits
    // nsegments in nsegments_bit_width bits
    // encoding size of xs in bit_offset_width bits
    // encoding size of ys in bit_offset_width bits
    // (nsegments*2+1) yshifts, each in yshift_bit_width(epsilon) bits
    // xs in Elias-Fano representation
    // ys in Elias-Fano representation

    using EFIterator = lemonhash::SequentialEliasFano::Iterator;

    static constexpr uint8_t key_bit_width = 64;
    static constexpr uint8_t epsilon_bit_width = 8;
    static constexpr uint8_t nsegments_bit_width = 24;
    static constexpr uint8_t bit_offset_width = 32;

public:

    EFSequentialPointsStorage() = default;

    EFSequentialPointsStorage(const EFSequentialPointsStorage &) = delete;

    EFSequentialPointsStorage(EFSequentialPointsStorage &&other) {
        data = other.data;
        other.data = nullptr;
    }

    EFSequentialPointsStorage &operator=(const EFSequentialPointsStorage &) = delete;

    EFSequentialPointsStorage &operator=(EFSequentialPointsStorage &&other) {
        if (this != &other) {
            data = other.data;
            other.data = nullptr;
        }
        return *this;
    }

    void load(uint64_t first, const auto &xs_data, auto &ys_data, const auto &yshifts_data, uint8_t epsilon) {
        for (size_t i = 0; i < ys_data.size(); ++i)
            ys_data[i] -= i;
        auto xs_bit_size = lemonhash::SequentialEliasFano::encodingSize(xs_data);
        auto ys_bit_size = lemonhash::SequentialEliasFano::encodingSize(ys_data);

        auto yshift_width = yshift_bit_width(epsilon);
        auto total_bit_size = key_bit_width + epsilon_bit_width + nsegments_bit_width + bit_offset_width +
            bit_offset_width + xs_bit_size + ys_bit_size + yshifts_data.size() * yshift_width;
        data = new uint64_t[(total_bit_size + 63) / 64]();

        auto ptr = data;
        uint8_t word_offset = 0;
        auto nsegments = xs_data.size() - 1;
        sdsl::bits::write_int_and_move(ptr, first, word_offset, key_bit_width);
        sdsl::bits::write_int_and_move(ptr, epsilon, word_offset, epsilon_bit_width);
        sdsl::bits::write_int_and_move(ptr, nsegments, word_offset, nsegments_bit_width);
        sdsl::bits::write_int_and_move(ptr, xs_bit_size, word_offset, bit_offset_width);
        sdsl::bits::write_int_and_move(ptr, ys_bit_size, word_offset, bit_offset_width);
        assert(first_key() == first);
        assert(epsilon_value() == epsilon);
        assert(segments_count() == nsegments);
        assert(yshifts_offset() == size_t(ptr - data) * 64 + word_offset);

        for (size_t i = 0; i < yshifts_data.size(); ++i) {
            sdsl::bits::write_int_and_move(ptr, yshifts_data[i], word_offset, yshift_width);
            assert(get_yshift(i) == yshifts_data[i]);
        }

        assert(size_t(ptr - data) * 64 + word_offset == xs_offset());
        assert(size_t(ptr - data) * 64 + word_offset + xs_bit_size == ys_offset());
        lemonhash::SequentialEliasFano::write(data, (ptr - data) * 64 + word_offset, xs_data);
        lemonhash::SequentialEliasFano::write(data, (ptr - data) * 64 + word_offset + xs_bit_size, ys_data);
        assert(size_in_bytes() / sizeof(uint64_t) == (total_bit_size + 63) / 64);
    }

    [[nodiscard]] std::tuple<size_t, uint64_t, uint64_t, int64_t, int64_t, int64_t> segment_for_key(auto key) const {
        EFIterator it(data, xs_offset(), xs_size());
        do {
            ++it;
        } while (it.index() + 1 < xs_size() && *it <= key);

        auto x1 = *it;
        --it;
        auto x0 = *it;
        auto [y0, y1, y2] = get_ys(it.index());
        return {it.index(), x0, x1, y0, y1, y2};
    }

    [[nodiscard]] size_t size() const {
        return EFIterator(data, ys_offset(), ys_size()).at(ys_size() - 1) + ys_size() - 1;
    }

    [[nodiscard]] size_t size_in_bytes() const {
        auto bit_offset = key_bit_width + epsilon_bit_width + nsegments_bit_width + bit_offset_width;
        auto bit_size = ys_offset() + sdsl::bits::read_int(data + bit_offset / 64, bit_offset % 64, bit_offset_width);
        return ((bit_size + 63) / 64) * sizeof(uint64_t);
    }

    [[nodiscard]] uint8_t epsilon_value() const {
        auto bit_offset = key_bit_width;
        return sdsl::bits::read_int(data + bit_offset / 64, bit_offset % 64, epsilon_bit_width);
    }

    [[nodiscard]] size_t segments_count() const {
        auto bit_offset = key_bit_width + epsilon_bit_width;
        return sdsl::bits::read_int(data + bit_offset / 64, bit_offset % 64, nsegments_bit_width);
    }

    [[nodiscard]] uint64_t first_key() const {
        return sdsl::bits::read_int(data, 0, key_bit_width);
    }

    ~EFSequentialPointsStorage() { delete[] data; }

    class Iterator {
        EFIterator xs_it;
        EFIterator ys_it;
        const EFSequentialPointsStorage *storage;

    public:

        Iterator(EFIterator xs_it, EFIterator ys_it, const EFSequentialPointsStorage *storage)
            : xs_it(xs_it), ys_it(ys_it), storage(storage) {}

        std::tuple<uint64_t, uint64_t, int64_t, int64_t, int64_t> operator*() {
            auto tmp1 = xs_it;
            auto x0 = *tmp1;
            auto x1 = *++tmp1;
            auto tmp2 = ys_it;
            auto y0 = *tmp2 + tmp2.index();
            auto y1 = *++tmp2 + tmp2.index();
            auto point_index = xs_it.index();
            int64_t epsilon = storage->epsilon_value();
            auto shift1 = int64_t(storage->get_yshift(point_index * 2)) - epsilon;
            auto shift2 = int64_t(storage->get_yshift(point_index * 2 + 1)) - epsilon;
            auto shift3 = int64_t(storage->get_yshift(point_index * 2 + 2)) - epsilon;
            return {x0, x1, y0 - shift1, y1 - shift2, y1 - shift3};
        }

        void operator++() {
            ++xs_it;
            ++ys_it;
        }
    };

    Iterator begin() const {
        return {EFIterator(data, xs_offset(), xs_size()), EFIterator(data, ys_offset(), ys_size()), this};
    }

private:

    size_t xs_size() const { return segments_count() + 1; }

    size_t ys_size() const { return segments_count() + 2; }

    size_t yshifts_size() const { return segments_count() * 2 + 1; }

    size_t yshifts_offset() const {
        return key_bit_width + epsilon_bit_width + nsegments_bit_width + 2 * bit_offset_width;
    }

    size_t xs_offset() const { return yshifts_offset() + yshifts_size() * yshift_bit_width(epsilon_value()); }

    size_t ys_offset() const {
        auto bit_offset = key_bit_width + epsilon_bit_width + nsegments_bit_width;
        return xs_offset() + sdsl::bits::read_int(data + bit_offset / 64, bit_offset % 64, bit_offset_width);
    }

    uint64_t get_yshift(size_t index) const {
        auto bit_offset = yshifts_offset() + index * yshift_bit_width(epsilon_value());
        return sdsl::bits::read_int(data + bit_offset / 64, bit_offset % 64, yshift_bit_width(epsilon_value()));
    }

    [[nodiscard]] std::tuple<int64_t, int64_t, int64_t> get_ys(size_t point_index) const {
        EFIterator it(data, ys_offset(), ys_size());
        int64_t y0 = it.at(point_index) + point_index;
        int64_t y1 = *++it + point_index + 1;
        int64_t epsilon = epsilon_value();
        auto shift1 = int64_t(get_yshift(point_index * 2)) - epsilon;
        auto shift2 = int64_t(get_yshift(point_index * 2 + 1)) - epsilon;
        auto shift3 = int64_t(get_yshift(point_index * 2 + 2)) - epsilon;
        return {y0 - shift1, y1 - shift2, y1 - shift3};
    }
};

#pragma pack(pop)

}