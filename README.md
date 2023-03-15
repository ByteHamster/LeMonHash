# LeMonHash: Learned Monotone Minimal Perfect Hashing

Monotone Minimal Perfect Hashing through combining
the [PGM-Index](https://github.com/gvinciguerra/PGM-index) for space-efficient ranking
and [BuRR](https://github.com/lorenzhs/BuRR) for low-overhead retrieval.

A variant for variable-length strings provides significantly faster queries than competitors.

## Requirements

- GCC 11 or later
- [libxxhash v0.8.0](https://github.com/Cyan4973/xxHash/releases/tag/v0.8.0) or later