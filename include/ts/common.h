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

struct SearchBound {
  size_t begin;
  size_t end;  // Exclusive.
};

struct TreeStats {
  TreeStats() {}
  
  TreeStats(unsigned treeMaxError, double cost, unsigned numNodes = 0)
    : treeMaxError(treeMaxError),
      cost(cost),
      numNodes(numNodes) {}

  unsigned treeMaxError;
  double cost;
  unsigned numNodes;
};

}  // namespace ts
