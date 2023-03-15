#pragma once

#include <compact_elias_fano.hpp>
#include <vector>
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>
#include "support/sequence/EliasFanoModified.hpp"

/**
 * Stores a monotonic list of integers.
 * Filters out runs of repeated same values using rank/select.
 * Relays actual storage of the integers to the child data structure.
 */
template <typename ChildSequence, size_t FILTER_GRANULARITY = 3>
class DuplicateFilterRank {
    private:
        size_t n = 0;
        uint64_t lastValueOfPreviousArea = 0;
        std::vector<size_t> insertionBuffer;
        pasta::BitVector duplicateFilterIsAreaStored;
        pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES> duplicateFilterRank;
        ChildSequence child;
    public:
        DuplicateFilterRank() {
            duplicateFilterIsAreaStored.resize(100, false);
            duplicateFilterIsAreaStored[0] = true;
        }

        DuplicateFilterRank(DuplicateFilterRank &&other) noexcept {
            std::construct_at(&child, std::move(other.child));
            n = other.n;
            lastValueOfPreviousArea = other.lastValueOfPreviousArea;
            insertionBuffer = std::move(other.insertionBuffer);
            duplicateFilterIsAreaStored = std::move(other.duplicateFilterIsAreaStored);
            duplicateFilterRank = std::move(other.duplicateFilterRank);
        }

        void append(size_t x) {
            n++;
            if (n <= FILTER_GRANULARITY) {
                child.append(x);
                lastValueOfPreviousArea = x;
            } else {
                insertionBuffer.push_back(x);
                flushInsertionBuffer(false);
            }
        }

        void flushInsertionBuffer(bool forceFlush) {
            if (insertionBuffer.size() < FILTER_GRANULARITY && !forceFlush) {
                return;
            }
            if (insertionBuffer.empty()) {
                return;
            }
            bool allSame = true;
            for (size_t i : insertionBuffer) {
                allSame = allSame && i == lastValueOfPreviousArea;
            }
            size_t area = (n - 1) / FILTER_GRANULARITY;
            if (!allSame || forceFlush) {
                for (size_t i : insertionBuffer) {
                    child.append(i);
                    lastValueOfPreviousArea = i;
                }
                if (duplicateFilterIsAreaStored.size() < area + 1) {
                    duplicateFilterIsAreaStored.resize(area * 2, false);
                }
                duplicateFilterIsAreaStored[area] = true;
            }
            insertionBuffer.clear();
        }

        size_t size() const {
            return n;
        }

        void complete() {
            if (n == 0) {
                return;
            }
            flushInsertionBuffer(true);
            duplicateFilterIsAreaStored.resize(n / FILTER_GRANULARITY + 1, false);
            duplicateFilterRank = pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES>(duplicateFilterIsAreaStored);
            child.complete();
        }

        size_t at(size_t i) {
            assert(i < n);
            size_t area = i / FILTER_GRANULARITY;
            size_t areaLocation = FILTER_GRANULARITY * duplicateFilterRank.rank1(area);
            size_t index = areaLocation;
            bool isAreaStored = duplicateFilterIsAreaStored[area];
            if (isAreaStored) {
                index += i % FILTER_GRANULARITY;
            } else {
                // Take last key of previous area
                index--;
            }
            assert(index <= i);
            return child.at(index);
        }

        size_t bit_size() {
            return 8 * (sizeof(*this) - sizeof(ChildSequence))
                    + duplicateFilterIsAreaStored.size()
                    + 8 * duplicateFilterRank.space_usage()
                    + child.bit_size();
        }
};

/**
 * Stores the positions and offsets at which runs of duplicate values start.
 * Runs need a size of at least MIN_DUPLICATES.
 */
template <typename ChildSequence>
class DuplicateFilterRunlength {
    private:
        static constexpr size_t MIN_DUPLICATES = 8;
        size_t n = 0;
        size_t previousValue = -1;
        size_t duplicatesOfPreviousValue = 0;
        ChildSequence child;
        std::vector<size_t> duplicateRunsStart;
        util::EliasFanoM *duplicateRunsStartEf = nullptr;
        std::vector<size_t> duplicateRunsLengthPrefix;
        util::EliasFanoM *duplicateRunsLengthPrefixEf = nullptr;
        size_t duplicateRunsLengthPrefixTmp = 0;
    public:
        DuplicateFilterRunlength() = default;

        DuplicateFilterRunlength(DuplicateFilterRunlength &&other) noexcept {
            std::construct_at(&child, std::move(other.child));
            n = other.n;
            previousValue = other.previousValue;
            duplicatesOfPreviousValue = other.duplicatesOfPreviousValue;
            duplicateRunsStart = std::move(other.duplicateRunsStart);
            duplicateRunsStartEf = other.duplicateRunsStartEf;
            other.duplicateRunsStartEf = nullptr;
            duplicateRunsLengthPrefix = std::move(other.duplicateRunsLengthPrefix);
            duplicateRunsLengthPrefixEf = other.duplicateRunsLengthPrefixEf;
            other.duplicateRunsLengthPrefixEf = nullptr;
            duplicateRunsLengthPrefixTmp = other.duplicateRunsLengthPrefixTmp;
        }

        void flushDuplicates() {
            if (duplicatesOfPreviousValue < MIN_DUPLICATES) {
                // Simply store them
                for (size_t i = 0; i < duplicatesOfPreviousValue; i++) {
                    child.append(previousValue);
                }
            } else {
                duplicateRunsStart.push_back(n - duplicatesOfPreviousValue);
                duplicateRunsLengthPrefix.push_back(duplicateRunsLengthPrefixTmp);
                duplicateRunsLengthPrefixTmp += duplicatesOfPreviousValue;
            }
        }

        void append(size_t x) {
            if (x == previousValue) {
                duplicatesOfPreviousValue++;
            } else {
                if (n != 0) {
                    flushDuplicates();
                }
                child.append(x);
                previousValue = x;
                duplicatesOfPreviousValue = 0;
            }
            n++;
        }

        size_t size() const {
            return n;
        }

        void complete() {
            if (n == 0) {
                return;
            }
            flushDuplicates();
            duplicateRunsStart.push_back(n);
            duplicateRunsLengthPrefix.push_back(duplicateRunsLengthPrefixTmp);

            duplicateRunsStartEf = new util::EliasFanoM(duplicateRunsStart.size(),
                               duplicateRunsStart.back() + 1);
            for (size_t i = 0; i < duplicateRunsStart.size(); i++) {
                duplicateRunsStartEf->push_back(duplicateRunsStart.at(i));
            }
            duplicateRunsStartEf->buildRankSelect();
            duplicateRunsStart.resize(0);
            duplicateRunsStart.shrink_to_fit();

            duplicateRunsLengthPrefixEf = new util::EliasFanoM(duplicateRunsLengthPrefix.size(),
                                duplicateRunsLengthPrefix.back() - (duplicateRunsLengthPrefix.size() - 1) * MIN_DUPLICATES + 1);
            for (size_t i = 0; i < duplicateRunsLengthPrefix.size(); i++) {
                duplicateRunsLengthPrefixEf->push_back(duplicateRunsLengthPrefix.at(i) - i * MIN_DUPLICATES);
            }
            duplicateRunsLengthPrefixEf->buildRankSelect();
            duplicateRunsLengthPrefix.resize(0);
            duplicateRunsLengthPrefix.shrink_to_fit();

            child.complete();
        }

        size_t at(size_t i) {
            if (duplicateRunsStartEf->size() == 0 || i < *duplicateRunsStartEf->at(0)) {
                return child.at(i); // Before any run
            }
            auto predecessorPtr = duplicateRunsStartEf->predecessorPosition(i);
            size_t predecessorIndex = predecessorPtr.index();
            size_t runStart = *predecessorPtr;

            auto numEmptyBeforePtr = duplicateRunsLengthPrefixEf->at(predecessorIndex);
            size_t removedEmptyBefore = *numEmptyBeforePtr + predecessorIndex * MIN_DUPLICATES;
            ++numEmptyBeforePtr;
            size_t removedEmptyAfter = *numEmptyBeforePtr + (predecessorIndex + 1) * MIN_DUPLICATES;
            size_t runStorageStart = runStart - removedEmptyBefore - 1;
            size_t runSize = removedEmptyAfter - removedEmptyBefore;
            size_t posInRun = i - runStart;
            if (posInRun < runSize) {
                return child.at(runStorageStart);
            }
            return child.at(runStorageStart + (posInRun - runSize) + 1);
        }

        size_t bit_size() {
            return 8 * (sizeof(*this) - sizeof(ChildSequence))
                   + 8 * duplicateRunsStartEf->space()
                   + 8 * duplicateRunsLengthPrefixEf->space()
                   + child.bit_size();
        }
};
