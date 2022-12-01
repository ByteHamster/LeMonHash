#pragma once

#include <SimpleRibbon.h>
#include "chunking_strategy/support/util.h"

/**
 * Container with multiple retrieval data structures of different sizes.
 * Can be indexed by specifying the *range* of the values to store/query.
 * An alternative method would be to specify the number of bits of each object.
 */
class MultiRetrievalDataStructure {
    private:
        SimpleRibbon<1> *retrieval1;
        SimpleRibbon<2> *retrieval2;
        SimpleRibbon<3> *retrieval3;
        SimpleRibbon<4> *retrieval4;
        SimpleRibbon<5> *retrieval5;
        SimpleRibbon<6> *retrieval6;
        SimpleRibbon<7> *retrieval7;
        SimpleRibbon<8> *retrieval8;
        SimpleRibbon<9, 64, uint16_t> *retrieval9;
        SimpleRibbon<10, 64, uint16_t> *retrieval10;
        SimpleRibbon<11, 64, uint16_t> *retrieval11;
        SimpleRibbon<12, 64, uint16_t> *retrieval12;
        SimpleRibbon<32, 64, uint32_t> *retrieval13Plus;

        std::vector<std::pair<uint64_t, uint8_t>> retrieval1Input;
        std::vector<std::pair<uint64_t, uint8_t>> retrieval2Input;
        std::vector<std::pair<uint64_t, uint8_t>> retrieval3Input;
        std::vector<std::pair<uint64_t, uint8_t>> retrieval4Input;
        std::vector<std::pair<uint64_t, uint8_t>> retrieval5Input;
        std::vector<std::pair<uint64_t, uint8_t>> retrieval6Input;
        std::vector<std::pair<uint64_t, uint8_t>> retrieval7Input;
        std::vector<std::pair<uint64_t, uint8_t>> retrieval8Input;
        std::vector<std::pair<uint64_t, uint16_t>> retrieval9Input;
        std::vector<std::pair<uint64_t, uint16_t>> retrieval10Input;
        std::vector<std::pair<uint64_t, uint16_t>> retrieval11Input;
        std::vector<std::pair<uint64_t, uint16_t>> retrieval12Input;
        std::vector<std::pair<uint64_t, uint32_t>> retrieval13PlusInput;
    public:
        void addInput(size_t range, uint64_t key, uint32_t value) {
            if (range <= (1 << 1)) {
                retrieval1Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 2)) {
                retrieval2Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 3)) {
                retrieval3Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 4)) {
                retrieval4Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 5)) {
                retrieval5Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 6)) {
                retrieval6Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 7)) {
                retrieval7Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 8)) {
                retrieval8Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 9)) {
                retrieval9Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 10)) {
                retrieval10Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 11)) {
                retrieval11Input.emplace_back(std::make_pair(key, value));
            } else if (range <= (1 << 12)) {
                retrieval12Input.emplace_back(std::make_pair(key, value));
            } else {
                retrieval13PlusInput.emplace_back(std::make_pair(key, value));
            }
        }

        void build() {
            retrieval1 = new SimpleRibbon<1>(retrieval1Input);
            retrieval2 = new SimpleRibbon<2>(retrieval2Input);
            retrieval3 = new SimpleRibbon<3>(retrieval3Input);
            retrieval4 = new SimpleRibbon<4>(retrieval4Input);
            retrieval5 = new SimpleRibbon<5>(retrieval5Input);
            retrieval6 = new SimpleRibbon<6>(retrieval6Input);
            retrieval7 = new SimpleRibbon<7>(retrieval7Input);
            retrieval8 = new SimpleRibbon<8>(retrieval8Input);
            retrieval9 = new SimpleRibbon<9, 64, uint16_t>(retrieval9Input);
            retrieval10 = new SimpleRibbon<10, 64, uint16_t>(retrieval10Input);
            retrieval11 = new SimpleRibbon<11, 64, uint16_t>(retrieval11Input);
            retrieval12 = new SimpleRibbon<12, 64, uint16_t>(retrieval12Input);
            retrieval13Plus = new SimpleRibbon<32, 64, uint32_t>(retrieval13PlusInput);

            retrieval1Input.clear();
            retrieval2Input.clear();
            retrieval3Input.clear();
            retrieval4Input.clear();
            retrieval5Input.clear();
            retrieval6Input.clear();
            retrieval7Input.clear();
            retrieval8Input.clear();
            retrieval9Input.clear();
            retrieval10Input.clear();
            retrieval11Input.clear();
            retrieval12Input.clear();
            retrieval13PlusInput.clear();

            retrieval1Input.shrink_to_fit();
            retrieval2Input.shrink_to_fit();
            retrieval3Input.shrink_to_fit();
            retrieval4Input.shrink_to_fit();
            retrieval5Input.shrink_to_fit();
            retrieval6Input.shrink_to_fit();
            retrieval7Input.shrink_to_fit();
            retrieval8Input.shrink_to_fit();
            retrieval9Input.shrink_to_fit();
            retrieval10Input.shrink_to_fit();
            retrieval11Input.shrink_to_fit();
            retrieval12Input.shrink_to_fit();
            retrieval13PlusInput.shrink_to_fit();
        }

        size_t query(size_t range, uint64_t key) {
            if (range <= (1 << 1)) {
                return retrieval1->retrieve(key);
            } else if (range <= (1 << 2)) {
                return retrieval2->retrieve(key);
            } else if (range <= (1 << 3)) {
                return retrieval3->retrieve(key);
            } else if (range <= (1 << 4)) {
                return retrieval4->retrieve(key);
            } else if (range <= (1 << 5)) {
                return retrieval5->retrieve(key);
            } else if (range <= (1 << 6)) {
                return retrieval6->retrieve(key);
            } else if (range <= (1 << 7)) {
                return retrieval7->retrieve(key);
            } else if (range <= (1 << 8)) {
                return retrieval8->retrieve(key);
            } else if (range <= (1 << 9)) {
                return retrieval9->retrieve(key);
            } else if (range <= (1 << 10)) {
                return retrieval10->retrieve(key);
            } else if (range <= (1 << 11)) {
                return retrieval11->retrieve(key);
            } else if (range <= (1 << 12)) {
                return retrieval12->retrieve(key);
            } else {
                return retrieval13Plus->retrieve(key);
            }
        }

        ~MultiRetrievalDataStructure() {
            delete retrieval1;
            delete retrieval2;
            delete retrieval3;
            delete retrieval4;
            delete retrieval5;
            delete retrieval6;
            delete retrieval7;
            delete retrieval8;
            delete retrieval9;
            delete retrieval10;
            delete retrieval11;
            delete retrieval12;
            delete retrieval13Plus;
        }

        /**
         * @param N If given, print size relative to N
         */
        size_t spaceBits(size_t N = 0) {
            if (N != 0) {
                std::cout<<"Space usage (bits/key) of the retrieval data structures:"<<std::endl;
                std::cout<<"1-bit: "<<(8.0*retrieval1->size()/N)<<", ";
                std::cout<<"2-bit: "<<(8.0*retrieval2->size()/N)<<", ";
                std::cout<<"3-bit: "<<(8.0*retrieval3->size()/N)<<", ";
                std::cout<<"4-bit: "<<(8.0*retrieval4->size()/N)<<std::endl;
                std::cout<<"5-bit: "<<(8.0*retrieval5->size()/N)<<", ";
                std::cout<<"6-bit: "<<(8.0*retrieval6->size()/N)<<", ";
                std::cout<<"7-bit: "<<(8.0*retrieval7->size()/N)<<", ";
                std::cout<<"8-bit: "<<(8.0*retrieval8->size()/N)<<", ";
                std::cout<<"9-bit: "<<(8.0*retrieval9->size()/N)<<", ";
                std::cout<<"10-bit: "<<(8.0*retrieval10->size()/N)<<", ";
                std::cout<<"11-bit: "<<(8.0*retrieval11->size()/N)<<", ";
                std::cout<<"12-bit: "<<(8.0*retrieval12->size()/N)<<", ";
                std::cout<<"13+-bit: "<<(8.0*retrieval13Plus->size()/N)<<std::endl;
            }
            size_t bytes = retrieval1->size() + retrieval2->size() + retrieval3->size()
                           + retrieval4->size() + retrieval5->size()
                           + retrieval6->size() + retrieval7->size()
                           + retrieval8->size() + retrieval9->size()
                           + retrieval10->size() + retrieval11->size()
                           + retrieval12->size() + retrieval13Plus->size();
            return (bytes + sizeof(*this)) * 8;
        }
};