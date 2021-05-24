#pragma once

#include <cassert>
#include <cmath>
#include <limits>
#include <map>
#include <optional>

#include "cht/builder.h"
#include "cht/cht.h"
#include "common.h"
#include "ts.h"

namespace ts {

// Allows building an adaptive `TrieSpline`.
template <class KeyType>
class Builder {
 public:
  Builder(KeyType min_key, KeyType max_key, size_t num_radix_bits, size_t spline_max_error,
                                            size_t cht_num_bins = 1024, size_t cht_max_error = 16)
      : min_key_(min_key),
        max_key_(max_key),
        num_radix_bits_(num_radix_bits),
        num_shift_bits_(GetNumShiftBits(max_key - min_key, num_radix_bits)),
        spline_max_error_(spline_max_error),
        cht_num_bins_(cht_num_bins),
        cht_max_error_(cht_max_error),
        curr_num_keys_(0),
        curr_num_distinct_keys_(0),
        prev_key_(min_key),
        prev_position_(0),
        prev_prefix_(0),
        chtb_(min_key, max_key, cht_num_bins, cht_max_error) {
          // Initialize radix configurations.
          radix_analyzer_.resize(1 + maxNumRadixBits);
          for (unsigned index = 1; index <= maxNumRadixBits; ++index) {
            radix_analyzer_[index].shiftBits = GetNumShiftBits(max_key_ - min_key_, index);
            radix_analyzer_[index].prevPrefix = 0;
            radix_analyzer_[index].prevSplineIndex = 0;
            radix_analyzer_[index].cost = 0;       
          }

          // And initialize radix table, needs to contain all prefixes up to the largest
          // key + 1.
          const uint32_t max_prefix = (max_key - min_key) >> num_shift_bits_;
          std::cerr << "max_prefix=" << max_prefix << std::endl;
          radix_table_.resize(max_prefix + 2, 0);
        }

  // Adds a key. Assumes that keys are stored in a dense array.
  void AddKey(KeyType key) {
    if (curr_num_keys_ == 0) {
      AddKey(key, /*position=*/0);
      return;
    }
    AddKey(key, prev_position_ + 1);
  }

  // Finalizes the construction and returns a read-only `TrieSpline`.
  TrieSpline<KeyType> Finalize() {
    // Last key needs to be equal to `max_key_`.
    assert(curr_num_keys_ == 0 || prev_key_ == max_key_);

    // Ensure that `prev_key_` (== `max_key_`) is last key on spline.
    if (curr_num_keys_ > 0 && spline_points_.back().x != prev_key_)
      AddKeyToSpline(prev_key_, prev_position_);

    // Finalyze the radix configurations.
    FinalizeRadixConfigurations();

    std::cerr << "slope=" << ComputeSlope() << std::endl;
    std::cerr << "theo=" << ComputeTheoreticalSlope() << std::endl;

    // Should we keep the radix table?
    if (ComputeSlope() < ComputeTheoreticalSlope() + precision) {
      std::cerr << "radix table!" << std::endl;
      // Then finalize the radix table.
      FinalizeRadixTable();

      std::cerr << "finish radix table!" << std::endl;

      std::cerr << "start finalizing cht!" << std::endl;
      // Dummy CHT instance.  
      auto cht_ = chtb_.Finalize(false);
    

      std::cerr << "finish finalizing cht!" << std::endl;
      std::cerr << "Min_key=" << min_key_ << std::endl;
      std::cerr << "shift=" << num_shift_bits_ << std::endl;

      // And return.
      return TrieSpline<KeyType>(ts::Mode::RadixTable,
        min_key_, max_key_, curr_num_keys_, num_radix_bits_, num_shift_bits_,
        spline_max_error_, cht_num_bins_, cht_max_error_, std::move(radix_table_), std::move(cht_), std::move(spline_points_));
    }
    
    std::cerr << "chose cht!" << std::endl;

    // Otherwise, build CHT.
    for (const auto [key, _]: spline_points_)
      chtb_.AddKey(key);
    
    // Finalize CHT
    auto cht_ = chtb_.Finalize();

    // And return the read-only instance
    return TrieSpline<KeyType>(ts::Mode::CHT, min_key_, max_key_, curr_num_keys_,
                               num_radix_bits_, num_shift_bits_, spline_max_error_,
                               cht_num_bins_, cht_max_error_,
                               std::move(radix_table_), std::move(cht_),
                               std::move(spline_points_));
  }

 private:
  // The maximum number of radix bits allowed.
  static constexpr unsigned maxNumRadixBits = 30;

  // Returns the number of shift bits based on the `diff` between the largest
  // and the smallest key. KeyType == uint32_t.
  static size_t GetNumShiftBits(uint32_t diff, size_t num_radix_bits) {
    const uint32_t clz = __builtin_clz(diff);
    if ((32 - clz) < num_radix_bits) return 0;
    return 32 - num_radix_bits - clz;
  }
  
  // KeyType == uint64_t.
  static size_t GetNumShiftBits(uint64_t diff, size_t num_radix_bits) {
    const uint32_t clzl = __builtin_clzl(diff);
    if ((64 - clzl) < num_radix_bits) return 0;
    return 64 - num_radix_bits - clzl;
  }

  // `int(ceil(log_2(distance)))`. I
  static size_t ComputeCost(uint32_t value) {
		assert(value);
    if (value == 1) return 1;
    return 31 - __builtin_clz(value) + ((value & (value - 1)) != 0);
  }


  void AddKey(KeyType key, size_t position) {
    assert(key >= min_key_ && key <= max_key_);
    // Keys need to be monotonically increasing.
    assert(key >= prev_key_);
    // Positions need to be strictly monotonically increasing.
    assert(position == 0 || position > prev_position_);

    PossiblyAddKeyToSpline(key, position);

    ++curr_num_keys_;
    prev_key_ = key;
    prev_position_ = position;
  }

  void AddKeyToSpline(KeyType key, double position) {
    spline_points_.push_back({key, position});
    PossiblyAddKeyToRadixTable(key);
    UpdateRadixConfigurations(spline_points_.size() - 1);
  }

  enum Orientation { Collinear, CW, CCW };
  static constexpr double precision = std::numeric_limits<double>::epsilon();

  static Orientation ComputeOrientation(const double dx1, const double dy1,
                                        const double dx2, const double dy2) {
    const double expr = std::fma(dy1, dx2, -std::fma(dy2, dx1, 0));
    if (expr > precision)
      return Orientation::CW;
    else if (expr < -precision)
      return Orientation::CCW;
    return Orientation::Collinear;
  };

  void SetUpperLimit(KeyType key, double position) {
    upper_limit_ = {key, position};
  }
  void SetLowerLimit(KeyType key, double position) {
    lower_limit_ = {key, position};
  }
  void RememberPreviousCDFPoint(KeyType key, double position) {
    prev_point_ = {key, position};
  }

  // Implementation is based on `GreedySplineCorridor` from:
  // T. Neumann and S. Michel. Smooth interpolating histograms with error
  // guarantees. [BNCOD'08]
  void PossiblyAddKeyToSpline(KeyType key, double position) {
    if (curr_num_keys_ == 0) {
      // Add first CDF point to spline.
      AddKeyToSpline(key, position);
      ++curr_num_distinct_keys_;
      RememberPreviousCDFPoint(key, position);
      return;
    }

    if (key == prev_key_) {
      // No new CDF point if the key didn't change.
      return;
    }

    // New CDF point.
    ++curr_num_distinct_keys_;

    if (curr_num_distinct_keys_ == 2) {
      // Initialize `upper_limit_` and `lower_limit_` using the second CDF
      // point.
      SetUpperLimit(key, position + spline_max_error_);
      SetLowerLimit(key, (position < spline_max_error_)
                             ? 0
                             : position - spline_max_error_);
      RememberPreviousCDFPoint(key, position);
      return;
    }

    // `B` in algorithm.
    const Coord<KeyType>& last = spline_points_.back();

    // Compute current `upper_y` and `lower_y`.
    const double upper_y = position + spline_max_error_;
    const double lower_y =
        (position < spline_max_error_) ? 0 : position - spline_max_error_;

    // Compute differences.
    assert(upper_limit_.x >= last.x);
    assert(lower_limit_.x >= last.x);
    assert(key >= last.x);
    const double upper_limit_x_diff = upper_limit_.x - last.x;
    const double lower_limit_x_diff = lower_limit_.x - last.x;
    const double x_diff = key - last.x;

    assert(upper_limit_.y >= last.y);
    assert(position >= last.y);
    const double upper_limit_y_diff = upper_limit_.y - last.y;
    const double lower_limit_y_diff = lower_limit_.y - last.y;
    const double y_diff = position - last.y;

    // `prev_point_` is the previous point on the CDF and the next candidate to
    // be added to the spline. Hence, it should be different from the `last`
    // point on the spline.
    assert(prev_point_.x != last.x);

    // Do we cut the error corridor?
    if ((ComputeOrientation(upper_limit_x_diff, upper_limit_y_diff, x_diff,
                            y_diff) != Orientation::CW) ||
        (ComputeOrientation(lower_limit_x_diff, lower_limit_y_diff, x_diff,
                            y_diff) != Orientation::CCW)) {
      // Add previous CDF point to spline.
      AddKeyToSpline(prev_point_.x, prev_point_.y);

      // Update limits.
      SetUpperLimit(key, upper_y);
      SetLowerLimit(key, lower_y);
    } else {
      assert(upper_y >= last.y);
      const double upper_y_diff = upper_y - last.y;
      if (ComputeOrientation(upper_limit_x_diff, upper_limit_y_diff, x_diff,
                             upper_y_diff) == Orientation::CW) {
        SetUpperLimit(key, upper_y);
      }

      const double lower_y_diff = lower_y - last.y;
      if (ComputeOrientation(lower_limit_x_diff, lower_limit_y_diff, x_diff,
                             lower_y_diff) == Orientation::CCW) {
        SetLowerLimit(key, lower_y);
      }
    }

    RememberPreviousCDFPoint(key, position);
  }

  void PossiblyAddKeyToRadixTable(KeyType key) {
    const KeyType curr_prefix = (key - min_key_) >> num_shift_bits_;
    if (curr_prefix != prev_prefix_) {
      const uint32_t curr_index = spline_points_.size() - 1;
      for (KeyType prefix = prev_prefix_ + 1; prefix <= curr_prefix; ++prefix)
        radix_table_[prefix] = curr_index;
      prev_prefix_ = curr_prefix;
    }
  }

  void FinalizeRadixTable() {
    ++prev_prefix_;
    const uint32_t num_spline_points = spline_points_.size();
    for (; prev_prefix_ < radix_table_.size(); ++prev_prefix_)
      radix_table_[prev_prefix_] = num_spline_points;
  }

  void UpdateRadixConfigurations(unsigned splineIndex) {
    // Ignore the first spline knot.
    if (!splineIndex) return;

    // And update the configurations for all possible radix tables.
    for (unsigned radix = 1; radix <= 30; ++radix) {
      const KeyType currPrefix = (spline_points_[splineIndex].x - min_key_) >> radix_analyzer_[radix].shiftBits;
      
      // New prefix?
      if (currPrefix != radix_analyzer_[radix].prevPrefix) {
        // Then compute statistics.
        assert(splineIndex);
        const unsigned prevSplineIndex = radix_analyzer_[radix].prevSplineIndex;
        const size_t numDataKeys = spline_points_[splineIndex].y - spline_points_[prevSplineIndex].y;
        const size_t numSplineKeys = splineIndex - prevSplineIndex;
        assert(numSplineKeys);

        // Update cost.
        radix_analyzer_[radix].cost += numDataKeys * ComputeCost(numSplineKeys);
        
        // Update the parameters.
        radix_analyzer_[radix].prevPrefix = currPrefix;
        radix_analyzer_[radix].prevSplineIndex = splineIndex;
      }
    }
  }

  void FinalizeRadixConfigurations() {
    for (unsigned radix = 1; radix <= 30; ++radix) {
      // Compute statistics.
      const unsigned prevSplineIndex = radix_analyzer_[radix].prevSplineIndex;
      const size_t numDataKeys = spline_points_.back().y - spline_points_[prevSplineIndex].y;
      const size_t numSplineKeys = spline_points_.size() - prevSplineIndex;
      assert(numSplineKeys);
      
      // Update the cost.
      radix_analyzer_[radix].cost += numDataKeys * ComputeCost(numSplineKeys);

      // Normalize the cost.
      radix_analyzer_[radix].cost /= spline_points_.back().y;
    }
  }

  // Compute the slope of radix configurations.
  // Least squares method.
  double ComputeSlope() {
    // Formula: slope := (n * sum(xy) - sum(x)sum(y)) / (nsum(x^2) - sum(x)^2). 
    unsigned n = 0;
    double sum_x_y = 0, sum_x = 0, sum_y = 0, sum_x_2 = 0;
    for (unsigned radix = 1; radix <= 30; ++radix) {
      // Close enough to 1.0 (optimal cost)?
      if (std::fabs(radix_analyzer_[radix].cost - 1) < precision)
        break;
      
      // And compute.
      sum_x_y += radix * radix_analyzer_[radix].cost;
      sum_x += radix;
      sum_y += radix_analyzer_[radix].cost;
      sum_x_2 += radix * radix;
      ++n;
    }
    return (n * sum_x_y - sum_x * sum_y) / (n * sum_x_2 - sum_x * sum_x);
  }

  double ComputeTheoreticalSlope() {
    // Build function f(x) = ax + b, s.t.:
    // f(1) = log_2(#spline-knots) - 1
    // f(numMaxRadixBits) = 1.0 (optimal cost)
    // Return the slope `a`.
    return -((1.0 * ComputeCost(spline_points_.size()) - 1) - 1) / (maxNumRadixBits - 1);

  }

  const KeyType min_key_;
  const KeyType max_key_;
  const size_t spline_max_error_;
  const size_t num_radix_bits_;
  const size_t num_shift_bits_;
  const size_t cht_num_bins_;
  const size_t cht_max_error_;
  
  std::vector<uint32_t> radix_table_;
  std::vector<Coord<KeyType>> spline_points_;

  size_t curr_num_keys_;
  size_t curr_num_distinct_keys_;
  KeyType prev_key_;
  size_t prev_position_;
  KeyType prev_prefix_;
  cht::Builder<KeyType> chtb_;

  // Current upper and lower limits on the error corridor of the spline.
  Coord<KeyType> upper_limit_;
  Coord<KeyType> lower_limit_;

  // Previous CDF point.
  Coord<KeyType> prev_point_;

  // Analyzer.
  std::vector<RadixConfig> radix_analyzer_;
};

}  // namespace ts
