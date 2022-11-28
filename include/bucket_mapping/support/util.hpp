#pragma once

/** Computes the number of bits needed to store x, that is, 0 if x is 0, 1 + floor(log2(x)) otherwise. */
#define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))
