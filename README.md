# Learned Monotone Minimal Perfect Hashing

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="img/lemon_wordmark_dark.png">
  <img src="img/lemon_wordmark.png" width="400" alt="Logo">
</picture>

LeMonHash is a Monotone Minimal Perfect Hash function that combines
the [PGM-Index](https://github.com/gvinciguerra/PGM-index) for space-efficient ranking
and [BuRR](https://github.com/lorenzhs/BuRR) for low-overhead retrieval.

A variant for variable-length strings provides significantly faster queries than competitors.

[<img src="img/plots.png" alt="Screenshot of measurements in paper">](https://arxiv.org/pdf/2304.11012)

### Requirements

- GCC 11 or later
- [libxxhash v0.8.0](https://github.com/Cyan4973/xxHash/releases/tag/v0.8.0) or later

### Usage

Clone the repository (as a submodule) and add the following to your `CMakeLists.txt`.

```cmake
add_subdirectory(path/to/LeMonHash)
target_link_libraries(YourTarget PRIVATE LeMonHash)
```

Then you can use the straight-forward interface of LeMonHash:

```cpp
std::vector<uint64_t> inputData {0, 1, 7, 15, 23, 42, 250};
lemonhash::LeMonHash<> hashFunc(inputData);
for (uint64_t x : inputData) {
    std::cout << x << ": \t" << hashFunc(x) << std::endl;
}
```

### License

This code is licensed under the [GPLv3](/LICENSE).
If you use the project in an academic context or publication, please cite our [paper](https://arxiv.org/pdf/2304.11012):

```bibtex
@article{ferragina2023learned,
    title   = {Learned Monotone Minimal Perfect Hashing},
    author  = {Paolo Ferragina and Hans-Peter Lehmann and Peter Sanders and Giorgio Vinciguerra},
    journal = {arXiv preprint arXiv:2304.11012},
    doi     = {10.48550/ARXIV.2304.11012},
    year    = {2023},
}
```

The code of the experiments comparing LeMonHash to competitors from the literature is available [here](https://github.com/ByteHamster/MMPHF-Experiments).
