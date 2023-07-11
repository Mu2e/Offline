#ifndef Mu2eUtilities_MedianCalculator_hh
#define Mu2eUtilities_MedianCalculator_hh
//
// this class is intended to be used for evaluating the median
// from a set of elements that are stored internally in a vector
//
// - intermediate answers can be retrieved before the full
// set of elements have been pushed
// - after clear(), it is ready for a new set of elements
//

#include <functional>
#include <vector>

namespace mu2e {
class MedianCalculator {

  struct MedianData {
    MedianData(float Val, float Wg) : val(Val), wg(Wg) {}
    float val;
    float wg;
  };
  struct lessByValue : public std::binary_function<MedianData, MedianData, bool> {
    bool operator()(MedianData const& p1, MedianData const& p2) { return p1.val < p2.val; }
  };

public:
  MedianCalculator(size_t nToReserve = 0) {
    clear();
    _vec.reserve(nToReserve);
  }

  float weightedMedian() { return median(true); }
  float unweightedMedian() { return median(false); }

  inline void push(float value, float weight = 1) {
    _vec.emplace_back(value, weight);
    _needsSorting = true;
    _totalWeight += weight;
  }

  inline size_t size() { return _vec.size(); }
  void clear();

private:
  float median(bool useWeights);

  std::vector<MedianData> _vec;
  bool _needsSorting;
  bool _goodWM;
  bool _goodUWM;
  float _weightedMedian;
  float _unweightedMedian;
  float _totalWeight;
};

} // namespace mu2e

#endif /* Mu2eUtilities_MedianCalculator_hh */
