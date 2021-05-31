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

// Allows building a `TrieSpline` in a single pass over sorted data
template <class KeyType>
class Builder {
 public:
  Builder(KeyType min_key, KeyType max_key, size_t spline_max_error,
          size_t num_bins, size_t tree_max_error, bool single_pass = false,
          bool use_cache = false)
      : min_key_(min_key),
        max_key_(max_key),
        spline_max_error_(spline_max_error),
        tree_max_error_(tree_max_error),
        curr_num_keys_(0),
        curr_num_distinct_keys_(0),
        prev_key_(min_key),
        prev_position_(0),
        chtb_(min_key, max_key, num_bins, tree_max_error, single_pass,
              use_cache) {}

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

    ComputeStatistics();

    // Finalize CHT
    auto cht_ = chtb_.Finalize();

    std::cerr << "finalize cht!" << std::endl;

    // And return the read-only instance
    return TrieSpline<KeyType>(min_key_, max_key_, curr_num_keys_,
                               spline_max_error_, tree_max_error_,
                               std::move(cht_), std::move(spline_points_));
  }

 private:
  using Interval = std::pair<unsigned, unsigned>;

  static unsigned computeLog(uint32_t n, bool round = false) {
    assert(n);
    return 31 - __builtin_clz(n) + (round ? ((n & (n - 1)) != 0) : 0);
  }

  static unsigned computeLog(uint64_t n, bool round = false) {
    assert(n);
    return 63 - __builtin_clzl(n) + (round ? ((n & (n - 1)) != 0) : 0);
  }

  static unsigned computeLcp(uint32_t x, uint32_t y) {
    return __builtin_clz(x ^ y);
  }

  static unsigned computeLcp(uint64_t x, uint64_t y) {
    return __builtin_clzl(x ^ y);
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
    AddKeyToCHT(key);
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

  void AddKeyToCHT(KeyType key) { chtb_.AddKey(key); }

#define DEBUG 0
  void ComputeStatistics() {
    unsigned maxBitLevel = 64;
    unsigned lg = computeLog(max_key_ - min_key_, true);
    unsigned alreadyCommon = (sizeof(KeyType) << 3) - lg;

    std::cerr << "lg=" << lg << std::endl;

    // TODO: use lg as maxBitLevel?
    // TODO: also check when min_key != 0!!!!

    // Compute the longest-common-prefix.
    const auto ExtractLCP = [&](unsigned index) -> unsigned {
      return computeLcp(spline_points_[index].x - min_key_, spline_points_[index - 1].x - min_key_) - alreadyCommon;// __builtin_clzl((spline_points_[index].x - min_key_) ^ (spline_points_[index - 1].x - min_key_)) - alreadyCommon;
    };

    // Fill the lcp-array.
    std::vector<unsigned> lcp(spline_points_.size());
    std::vector<unsigned> counters(maxBitLevel + 1);
    lcp[0] = std::numeric_limits<unsigned>::max();
    for (unsigned index = 1, limit = spline_points_.size(); index != limit; ++index) {
      lcp[index] = ExtractLCP(index);
      //std::cerr << "index=" << index << " x1=" << spline_points_[index-1].x << " x2=" << spline_points_[index].x << " lcp=" << lcp[index] << std::endl;
      counters[lcp[index]]++;
    }

    // TODO: don't forget that we have [1, size()[
    std::vector<unsigned> offsets(maxBitLevel + 2);
    for (unsigned bit = 1; bit <= maxBitLevel; ++bit) {
      offsets[bit] = offsets[bit - 1] + counters[bit - 1];
    }

    std::vector<unsigned> sorted(spline_points_.size() - 1);
    for (unsigned index = 1, limit = spline_points_.size(); index != limit; ++index) {
      sorted[offsets[lcp[index]]++] = index;
    }
    //std::cerr << "offsets[lg]=" << offsets[lg] << " vs " << spline_points_.size() << std::endl;
    //assert(offsets[lg] == spline_points_.size());

    std::pair<unsigned, std::vector<Interval>> histogram[2];
    histogram[0].first = histogram[1].first = 0;
    histogram[0].second.resize(spline_points_.size());
    histogram[1].second.resize(spline_points_.size());
    unsigned side = 0;

    const auto AddNewInterval = [&](Interval interval) -> void {
      // Empty interval?
      if (interval.first == interval.second)
        return;
      histogram[side].second[histogram[side].first++] = interval;
    };


    // Init the first level of the histogram and move already to the next level.
    AddNewInterval({1, spline_points_.size()});
    unsigned ptrInSorted = 0;

    const auto print = [&](Interval interval) -> std::string {
      std::string ret = "[";
      ret += std::to_string(interval.first);
      ret += ", ";
      ret += std::to_string(interval.second);
      ret += "[";
      return ret;
    };

    const auto AnalyzeInterval = [&](unsigned level, Interval interval) -> void {
      const auto isInside = [&](unsigned pos) -> bool {
        return (interval.first <= pos) && (pos < interval.second);
      };
 #if DEBUG
      std::cerr << "[analyze interval] level=" << level << " interval=" << print(interval) << std::endl;
#endif

      //std::cerr << "ptrInSorted=" << ptrInSorted << " vs sorted.size()=" << sorted.size() << std::endl;
      // TODO: is this fine?
      if (ptrInSorted == sorted.size()) {
        std::cerr << "already here!" << std::endl;
        return;
      }

      // Could the next `lcp` be an interval-breaker?
      auto leftSide = interval.first;
      while ((ptrInSorted != sorted.size()) && (lcp[sorted[ptrInSorted]] < level)) {
        // Check if its position is inside our interval.
        // TODO: make sure that it's correct! (small example)
        
        // interval := [first, second[
        if (!isInside(sorted[ptrInSorted]))
          break;
        
        AddNewInterval({leftSide, sorted[ptrInSorted]});
        leftSide = sorted[ptrInSorted] + 1;
        ++ptrInSorted;
      }
      // TODO: check it's correct.
      // It can be the same interval!!!!!!
      AddNewInterval({leftSide, interval.second});
    };

    // TODO: mark numBins as unsigned! To avoid the same compilatoion error on mac!
    static constexpr unsigned maxNumPossibleBins = 20;
    unsigned numPossibleBins = std::min(maxNumPossibleBins, lg);
    static constexpr unsigned maxPossibleTreeError = 1u << 10;
    std::vector<unsigned> possibleNumBins(numPossibleBins);
    std::vector<std::vector<uint64_t>> matrix(numPossibleBins);
    for (unsigned index = 0; index != numPossibleBins; ++index) {
      possibleNumBins[index] = (1u << (index + 1));
      matrix[index].assign(1 + maxPossibleTreeError, 0);
    }

    const auto ConsumeLevel = [&](unsigned level) -> void {
      for (unsigned index = 0; index != numPossibleBins; ++index) {
        auto currNumBins = possibleNumBins[index];
        

        // Does this number of bins benefit from this level?
        // TODO: optimize this condition!!! (preprocess the number of bins for each level)
        //std::cerr << "currNumBins=" << currNumBins << " log=" << computeLog(currNumBins) << std::endl;
        if (level % computeLog(currNumBins) == 0) {
#if DEBUG
          std::cerr << "Consume level=" << level << std::endl;
          std::cerr << "currNumbins=" << currNumBins << std::endl;
#endif
          // Consume all intervals.
          for (unsigned ptr = 0, limit = histogram[side].first; ptr != limit; ++ptr) {
            // [first, second[ also takes into consideration the `first-1`th element.
            // This is due to `lcp`-array, which takes the previous element into consideration.
            // That's why `second` - `first` + 1.
            assert(histogram[side].second[ptr].second > histogram[side].second[ptr].first);
            auto intervalSize = histogram[side].second[ptr].second - histogram[side].second[ptr].first + 1;
            matrix[index][std::min(intervalSize - 1, maxPossibleTreeError)] += intervalSize;
#if DEBUG
            std::cerr << "take interval=" << print(histogram[side].second[ptr]) << " intervalSize=" << intervalSize << std::endl;
#endif
          }
#if DEBUG
          std::cerr << "--------------------------------" << std::endl;
#endif
        }
      }
    };

    // TODO: until which bit? log2(num_bits) possible. Consider max - min?
    // TODO: really max lcp???
    for (unsigned level = 1; level <= lg; ++level) {
      std::cerr << "level=" << level << std::endl;
      side = 1 - side;
      histogram[side].first = 0;
      for (unsigned index = 0, limit = histogram[1 - side].first; index != limit; ++index) {
        AnalyzeInterval(level, histogram[1 - side].second[index]);
      }
      ConsumeLevel(level);
    }

#if 0
    for (unsigned index = 0; index != numPossibleBins; ++index) {
      std::cerr << "num_bins=" << possibleNumBins[index] << std::endl << "\t";
      for (unsigned ptr = 0; ptr <= maxPossibleTreeError; ++ptr) {
        std::cerr << matrix[index][ptr] << ", "; 
      }
      std::cerr << std::endl;
    }
#endif


    // Build the prefix sums.
    std::vector<std::pair<unsigned, double>> bestConfigurations(numPossibleBins);
    for (unsigned index = 0; index != numPossibleBins; ++index) {
      unsigned localBest = maxPossibleTreeError;
      uint64_t localMinCost = matrix[index][maxPossibleTreeError] + computeLog(maxPossibleTreeError, true);
      for (unsigned backPtr = maxPossibleTreeError; backPtr; --backPtr) {
        matrix[index][backPtr - 1] += matrix[index][backPtr];
        if (matrix[index][backPtr - 1] + computeLog(backPtr - 1, true) < localMinCost) {
          localMinCost = matrix[index][backPtr - 1] + computeLog(backPtr - 1, true);
          localBest = backPtr - 1;
        }
      }
      double localAvg = 1.0 * localMinCost / spline_points_.size();
      bestConfigurations[index] = {localBest, localAvg};
      std::cerr << "numBins=" << possibleNumBins[index] << " localAvg=" << localAvg << " localBest=" << localBest << std::endl;
    }

#if 0
    for (unsigned index = 0; index != numPossibleBins; ++index) {
      std::cerr << "num_bins=" << possibleNumBins[index] << std::endl << "\t";
      for (unsigned ptr = 0; ptr <= maxPossibleTreeError; ++ptr) {
        std::cerr << matrix[index][ptr] << ", "; 
      }
      std::cerr << std::endl;
    }
#endif

    std::cerr << "finished statistic!" << std::endl;
  }

  const KeyType min_key_;
  const KeyType max_key_;
  const size_t spline_max_error_;
  const size_t tree_max_error_;

  std::vector<Coord<KeyType>> spline_points_;

  size_t curr_num_keys_;
  size_t curr_num_distinct_keys_;
  KeyType prev_key_;
  size_t prev_position_;
  cht::Builder<KeyType> chtb_;

  // Current upper and lower limits on the error corridor of the spline.
  Coord<KeyType> upper_limit_;
  Coord<KeyType> lower_limit_;

  // Previous CDF point.
  Coord<KeyType> prev_point_;
};

}  // namespace ts
