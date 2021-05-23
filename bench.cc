#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <queue>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "include/ts/builder.h"
#include "include/ts/common.h"

using namespace std;
using namespace std::chrono;

static std::string hdd = "/media/mihail/seagate/";

auto randomNumberBetween = [](uint64_t low, uint64_t high) {
  auto randomFunc =
      [distribution_ = std::uniform_int_distribution<uint64_t>(low, high),
       random_engine_ = std::mt19937{std::random_device{}()}]() mutable {
        return distribution_(random_engine_);
      };
  return randomFunc;
};

static std::vector<uint64_t> generator(unsigned size, uint64_t maxValue,
                                       bool withSort = false) {
  vector<uint64_t> numbers;
  std::generate_n(std::back_inserter(numbers), size / 4,
                  randomNumberBetween(0, maxValue / 4));
  std::generate_n(std::back_inserter(numbers), size / 2,
                  randomNumberBetween(maxValue / 4, maxValue / 2));
  std::generate_n(std::back_inserter(numbers), size / 4,
                  randomNumberBetween(maxValue / 2, maxValue));
  if (withSort) {
    sort(std::begin(numbers), std::end(numbers));
#if 0
		numbers.erase(std::unique(numbers.begin(), numbers.end()), numbers.end());
    assert(std::is_sorted(numbers.begin(), numbers.end()));
#endif
  }
  return numbers;
}

void createFile(std::string name) {
  // TODO: by now only unique keys!
  static constexpr uint64_t numKeys = 5e7;
  static constexpr uint64_t maxValue = std::numeric_limits<uint64_t>::max();
  vector<uint64_t> keys = generator(numKeys, maxValue, true);
  std::ofstream outputKeys(hdd + name + ".k");
  outputKeys.write(reinterpret_cast<char*>(keys.data()),
                   keys.size() * sizeof(uint64_t));

  std::cerr << "Generated: " << keys.size() << " keys" << std::endl;

  vector<uint64_t> queries = generator(1e6, maxValue);
  std::ofstream outputQueries(hdd + name + ".q");
  outputQueries.write(reinterpret_cast<char*>(queries.data()),
                      queries.size() * sizeof(uint64_t));
}

void read(std::vector<uint64_t>& v, std::string fileName) {
  ifstream in(hdd + fileName);
  auto pos = in.tellg();
  in.seekg(0, ios::end);
  auto size = in.tellg() - pos;
  in.seekg(0, ios::beg);
  unsigned elements = size / sizeof(uint64_t);
  v.resize(elements);
  in.read(reinterpret_cast<char*>(v.data()), elements * sizeof(uint64_t));
}

void test(std::string name) {
  vector<uint64_t> keys, queries;
#if 1
  read(keys, name + ".k");
#else
  keys = {0,  0,  0,  0,  0,  0,  0,  2,  2,  3,  3,  3,  3,  4,  4,  4,  5,
          5,  6,  6,  6,  7,  7,  7,  7,  9,  9,  11, 13, 13, 14, 14, 14, 14,
          14, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 18, 19, 22, 23, 25};
#endif
  read(queries, name + ".q");
  assert(!keys.empty());
  assert(!queries.empty());
  std::cerr << "keys size=" << keys.size() << std::endl;
  std::cerr << "queries size=" << queries.size() << std::endl;

#if 0
	const auto debugKeys = [&]() {
		std::cerr << "Keys" << std::endl;
		std::cerr << "size=" << keys.size() << std::endl;
		
		auto index = 0;
		for (auto key : keys) {
			auto key_length = std::to_string(key).size();
			auto index_length = std::to_string(index).size();
			std::cerr << index;
			std::cerr << "|";
			++index;
		}
		std::cerr << std::endl;
		index = 0;
		for (auto key : keys) {
			auto key_length = std::to_string(key).size();
			auto index_length = std::to_string(index).size();
			//std::cerr << "key=" << key << " index=" << index << " key_length=" << key_length << " index_length=" << index_length << std::endl;
			for (unsigned i = 0, limit = (index_length > key_length) ? (index_length - key_length) : 0; i != limit; ++i)
				std::cerr << " ";
			std::cerr << key << " ";
			++index;
		}
		std::cerr << std::endl;
	};
	
	debugKeys();
#endif

  auto min = keys.front();
  auto max = keys.back();

  size_t spline_max_error_ = 32, num_bins_ = 64, tree_max_error_ = 16;
  ts::Builder<uint64_t> tsbSinglePassWithoutCache(
      min, max, spline_max_error_, num_bins_, tree_max_error_, true, false);
  ts::Builder<uint64_t> tsbMultiplePassWithoutCache(
      min, max, spline_max_error_, num_bins_, tree_max_error_, false, false);
  ts::Builder<uint64_t> tsbMultiplePassWithCache(
      min, max, spline_max_error_, num_bins_, tree_max_error_, false, true);

  for (auto key : keys) {
    tsbSinglePassWithoutCache.AddKey(key);
    tsbMultiplePassWithoutCache.AddKey(key);
    tsbMultiplePassWithCache.AddKey(key);
  }
  ts::TrieSpline tsSinglePassWithoutCache =
      tsbSinglePassWithoutCache.Finalize();
  ts::TrieSpline tsMultiplePassWithoutCache =
      tsbMultiplePassWithoutCache.Finalize();
  ts::TrieSpline tsMultiplePassWithCache = tsbMultiplePassWithCache.Finalize();

  const auto checkResult = [&](ts::SearchBound searchBound,
                               uint64_t query) -> bool {
    auto pos = std::lower_bound(keys.begin(), keys.end(), query) - keys.begin();
    auto rangePos = std::lower_bound(keys.begin() + searchBound.begin,
                                     keys.begin() + searchBound.end, query) -
                    keys.begin();
    assert(pos == rangePos);
    return ((searchBound.begin <= pos) && (pos <= searchBound.end));
  };

  const auto measureTime = [&](std::string type) -> void {
    auto start = high_resolution_clock::now();
    if (type == "ts-10") {
      for (uint64_t query : queries) {
        auto searchBound = tsSinglePassWithoutCache.GetSearchBound(query);
#if 0
				assert(checkResult(searchBound, query));
#endif
      }
    } else if (type == "ts-00") {
      for (uint64_t query : queries) {
        auto searchBound = tsMultiplePassWithoutCache.GetSearchBound(query);
#if 0
				assert(checkResult(searchBound, query));
#endif
      }
    } else if (type == "ts-01") {
      for (uint64_t query : queries) {
        auto searchBound = tsMultiplePassWithCache.GetSearchBound(query);
#if 0
				assert(checkResult(searchBound, query));
#endif
      }
    }

    auto stop = high_resolution_clock::now();
    auto answer =
        duration_cast<nanoseconds>(stop - start).count() / queries.size();
    std::cout << type << ": " << answer << " ns" << std::endl;
  };

  // ts-{`single_pass`}-{`use_cache`}
  measureTime("ts-10");
  measureTime("ts-00");
  measureTime("ts-01");
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <option>{0: createFile(s), 1: test}"
              << std::endl;
    exit(-1);
  }

  auto option = atoi(argv[1]);
  if (option == 0) {
    createFile("test");
  } else {
    test("test");
  }
  return 0;
}
