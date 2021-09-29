#include <algorithm>
#include <iostream>
#include <memory>

#include "cetlib_except/exception.h"

#include "Offline/Mu2eUtilities/inc/MedianCalculator.hh"

namespace mu2e {

//================================================================

void MedianCalculator::clear() {
  _vec.clear();
  _needsSorting = true;
  _goodWM = false;
  _goodUWM = false;
  _weightedMedian = 0;
  _unweightedMedian = 0;
  _totalWeight = 0;
}

//================================================================

float MedianCalculator::median(bool useWeights) {

  size_t v_size = _vec.size();
  if (v_size == 0) {
    throw cet::exception("MEDIANCALCULATOR_ZERO")
        << "No entries in the vector: median undefined" << std::endl;
  }

  if (v_size == 1) {
    return _vec[0].val;
  }

  if (_needsSorting) {
    std::sort(_vec.begin(), _vec.end(), lessByValue());
    _needsSorting = false;
    _goodWM = false;
    _goodUWM = false;
  }

  if (useWeights) {

    if (!_goodWM) {

      // integratng low to high, id will the first entry past the 50% point
      // the entry at id then satisfies the definition of weighted median:
      // the sum of weights above and below id (not including id)
      // is less than 50%.  id can be the first or last entry.

      size_t id = 0;
      float sum = _vec[0].wg;
      while (sum < 0.5 * _totalWeight) {
        ++id;
        sum += _vec[id].wg;
      }

      _weightedMedian = _vec[id].val;
      _goodWM = true;
    }

    return _weightedMedian;

  } else {

    if (!_goodUWM) {

      // v_size is >=2, so id >=1
      size_t id = v_size / 2;
      if (v_size % 2 == 0) {
        _unweightedMedian = (_vec[id - 1].val + _vec[id].val) / 2.0;
      } else {
        _unweightedMedian = _vec[id].val;
      }
      _goodUWM = true;
    }
    return _unweightedMedian;
  }
}

} // namespace mu2e
