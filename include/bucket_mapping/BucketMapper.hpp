#pragma once

#include <functional>

class BucketMapper {
    public:
        using Iterator = std::vector<uint64_t>::iterator;
        using Func = std::function<void(Iterator, size_t)>;

        virtual ~BucketMapper() = default;

        [[nodiscard]] virtual size_t bucketOf(uint64_t key) const = 0;

        virtual void bucketOf(Iterator first, Iterator last, Func f) const = 0;

        [[nodiscard]] virtual size_t size() const = 0;

        [[nodiscard]] virtual size_t numBuckets() const = 0;

};