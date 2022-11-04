#include <vector>
#include <algorithm>
#include <iostream>
#include <SimpleRibbon.h>

/**
 * Simple binary search on the raw input data.
 */
class BinarySearchMmphf {
    private:
        std::vector<uint64_t> data;
    public:
        explicit BinarySearchMmphf(std::vector<uint64_t> &data) : data(data) {
        }

        static std::string name() {
            return "BinarySearchMmphf";
        }

        size_t operator()(uint64_t key) {
            size_t a = 0;
            size_t b = data.size() - 1;
            while (a != b) {
                size_t p = (a + b + 1) / 2; // +1 ==> ceil
                if (data[p] > key) {
                    b = p - 1;
                } else {
                    a = p;
                }
            }
            return a;
        }

        size_t spaceBits() {
            return data.size() * sizeof(uint64_t) * 8 + sizeof(*this) * 8;
        }
};

/**
 * Explicitly store the ranks of all input objects using one large retrieval data structure.
 */
class DirectRetrievalMmphf {
    private:
        SimpleRibbon<32, 64, uint32_t> *retrieval;
    public:
        explicit DirectRetrievalMmphf(std::vector<uint64_t> &data) {
            std::vector<std::pair<uint64_t, uint32_t>> retrievalInput;
            for (size_t i = 0; i < data.size(); i++) {
                retrievalInput.emplace_back(std::make_pair(data.at(i), i));
            }
            retrieval = new SimpleRibbon<32, 64, uint32_t>(retrievalInput);
        }

        static std::string name() {
            return "DirectRetrievalMmphf";
        }

        size_t operator()(uint64_t key) {
            return retrieval->retrieve(key);
        }

        size_t spaceBits() {
            return retrieval->size() * 8 + sizeof(*this) * 8;
        }
};
