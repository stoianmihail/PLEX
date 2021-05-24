#pragma once

#include <cstddef>
#include <cstdint>

namespace ts {

// A CDF coordinate.
template <class KeyType>
struct Coord {
  KeyType x;
  double y;
};

// A radix config.
struct RadixConfig {
	unsigned shiftBits;
	unsigned prevPrefix;
	unsigned prevSplineIndex;
	double cost;
};

struct SearchBound {
  size_t begin;
  size_t end;  // Exclusive.
};

}  // namespace ts
