## PLEX: Towards Practical Learned Indexing

[PLEX](https://arxiv.org/abs/2108.05117) only has a single hyperparameter ϵ (maximum prediction error) and offers a better trade-off between build and lookup time than state-of-the-art approaches.

## Build

```
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./example
```

## Examples

Using ``ts::Builder`` to index sorted data:

```c++
// Create random keys.
std::vector<uint64_t> keys(1e6);
generate(keys.begin(), keys.end(), rand);
keys.push_back(424242);
std::sort(keys.begin(), keys.end());

// Build PLEX
uint64_t min = keys.front();
uint64_t max = keys.back();
ts::Builder tsb(min, max, keys.size());

for (const auto& key : keys) tsb.AddKey(key);
auto ts = tsb.Finalize();

// Search using PLEX
ts::SearchBound bound = ts.GetSearchBound(424242);
std::cout << "The search key is in the range: ["
			<< bound.begin << ", " << bound.end << ")" << std::endl;
auto start = std::begin(keys) + bound.begin, last = std::begin(keys) + bound.end;
auto pos = std::lower_bound(start, last, 424242) - begin(keys);
assert(keys[pos] == 424242);
std::cout << "The key is at position: " << pos << std::endl;
```

## Cite

Please cite our [AIDB@VLDB 2021 paper](https://arxiv.org/abs/2108.05117) if you use this code in your own work.
