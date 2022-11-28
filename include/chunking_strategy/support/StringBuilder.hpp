#pragma once

#include <string>
#include <cassert>

class StringBuilder {
    private:
        std::string string;
        unsigned char nextChar = 0;
        size_t nextCharPosition = 0;
    public:
        /**
         * Appends an integer, reversing the bytes (little endian) to keep the lexicographical order.
         * If the number of bits is not a multiple of 8, this method shifts the number (and all later numbers)
         * to tightly pack the integers.
         *
         * appendInt(0b11100, 5); appendInt(0b1101; 4); ==> 11100110 10000000
         */
        void appendInt(uint64_t x, size_t numBits) {
            assert(numBits < 48);
            uint64_t shiftBits = 16 - nextCharPosition - numBits;
            shiftBits = shiftBits % 8;

            x = x << shiftBits;
            long nextContentByte = static_cast<long>(numBits + shiftBits - 1) / 8;
            nextChar |= ((char*) &x)[nextContentByte];
            nextContentByte--;

            while (nextContentByte >= 0) {
                flush();
                nextChar = ((char*) &x)[nextContentByte];
                nextContentByte--;
            }
            nextCharPosition = 8 - shiftBits;
            if (nextCharPosition == 8) {
                flush();
            }
        }

        std::string toString() {
            if (nextCharPosition != 0) {
                flush();
            }
            return string;
        }

    private:
        void flush() {
            string += nextChar;
            nextChar = 0;
            nextCharPosition = 0;
        }
};