#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "cht/cht.h"
#include "common.h"

namespace ts {
enum class Mode : bool {
  RadixTable, CHT
};

// Approximates a cumulative distribution function (CDF) using spline
// interpolation.
template <class KeyType>
class TrieSpline {
 public:
  TrieSpline() = default;

  TrieSpline(ts::Mode mode, KeyType min_key, KeyType max_key, size_t num_keys,
             size_t num_radix_bits, size_t num_shift_bits,
             size_t spline_max_error, size_t cht_num_bins, size_t cht_max_error,
             std::vector<uint32_t> radix_table,
             cht::CompactHistTree<KeyType> cht,
             std::vector<ts::Coord<KeyType>> spline_points)
      : mode(mode),
        min_key_(min_key),
        max_key_(max_key),
        num_keys_(num_keys),
        num_radix_bits_(num_radix_bits),
        num_shift_bits_(num_shift_bits),
        spline_max_error_(spline_max_error),
        cht_num_bins_(cht_num_bins),
        cht_max_error_(cht_max_error),
        radix_table_(std::move(radix_table)),
        cht_(std::move(cht)),
        spline_points_(std::move(spline_points)) {
          std::cerr << "min_key=" << min_key << std::endl;
          std::cerr << "shift=" << num_shift_bits << std::endl;
        }

  // Returns the estimated position of `key`.
  double GetEstimatedPosition(const KeyType key) const {
    // Truncate to data boundaries.
    if (key <= min_key_) return 0;
    if (key >= max_key_) return num_keys_ - 1;

    // Find spline segment with `key` ∈ (spline[index - 1], spline[index]].
    const size_t index = GetSplineSegment(key);
    const Coord<KeyType> down = spline_points_[index - 1];
    const Coord<KeyType> up = spline_points_[index];

    // Compute slope.
    const double x_diff = up.x - down.x;
    const double y_diff = up.y - down.y;
    const double slope = y_diff / x_diff;

    // Interpolate.
    const double key_diff = key - down.x;
    return std::fma(key_diff, slope, down.y);
  }

  // Returns a search bound [begin, end) around the estimated position.
  ts::SearchBound GetSearchBound(const KeyType key) const {
    const size_t estimate = GetEstimatedPosition(key);
    const size_t begin =
        (estimate < spline_max_error_) ? 0 : (estimate - spline_max_error_);
    // `end` is exclusive.
    const size_t end = (estimate + spline_max_error_ + 2 > num_keys_)
                           ? num_keys_
                           : (estimate + spline_max_error_ + 2);
    return ts::SearchBound{begin, end};
  }

  // Returns the size in bytes.
  size_t GetSize() const {
    return sizeof(*this) + cht_.GetSize() +
           spline_points_.size() * sizeof(Coord<KeyType>);
  }

 private:
  // Returns the index of the spline point that marks the end of the spline
  // segment that contains the `key`: `key` ∈ (spline[index - 1], spline[index]]
  size_t GetSplineSegment(const KeyType key) const {
    return (mode == ts::Mode::RadixTable) ? RadixTableLookup(key) : CompactHistTreeLookup(key);
  }

  size_t CompactHistTreeLookup(const KeyType key) const {
    // Narrow search range.
    const auto range = cht_.GetSearchBound(key);
  
    if (cht_max_error_ < 32) {
      // Do linear search over narrowed range.
      uint32_t current = range.begin;
      while (spline_points_[current].x < key) ++current;
      return current;
    }

    // Do binary search over narrowed range.
    const auto lb =
        std::lower_bound(spline_points_.begin() + range.begin,
                         spline_points_.begin() + range.end, key,
                         [](const Coord<KeyType>& coord, const KeyType key) {
                           return coord.x < key;
                         });
    return std::distance(spline_points_.begin(), lb);
  }

  size_t RadixTableLookup(const KeyType key) const {
    // Narrow search range using radix table.
    const KeyType prefix = (key - min_key_) >> num_shift_bits_;
    assert(prefix + 1 < radix_table_.size());
    const uint32_t begin = radix_table_[prefix];
    const uint32_t end = radix_table_[prefix + 1];

    if (end - begin < 32) {
      // Do linear search over narrowed range.
      uint32_t current = begin;
      while (spline_points_[current].x < key) ++current;
      return current;
    }

    // Do binary search over narrowed range.
    const auto lb = std::lower_bound(
        spline_points_.begin() + begin, spline_points_.begin() + end, key,
        [](const Coord<KeyType>& coord, const KeyType key) {
          return coord.x < key;
        });
    return std::distance(spline_points_.begin(), lb);
  }

  ts::Mode mode;
  KeyType min_key_;
  KeyType max_key_;
  size_t num_keys_;
  size_t num_radix_bits_;
  size_t num_shift_bits_;
  size_t spline_max_error_;
  size_t cht_num_bins_;
  size_t cht_max_error_;

  std::vector<uint32_t> radix_table_;
  cht::CompactHistTree<KeyType> cht_;  
  std::vector<ts::Coord<KeyType>> spline_points_;
};

}  // namespace ts
