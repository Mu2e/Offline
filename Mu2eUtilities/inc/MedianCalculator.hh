#ifndef Mu2eUtilities_MedianCalculator_hh
#define Mu2eUtilities_MedianCalculator_hh
//
// Original author G. Pezzullo
//
// this class is intended to be used for evaluaitng the median 
// from a set of elements that are stored internally in a vector
//

#include <stddef.h>
#include <functional>
//#include <utility>
#include <numeric>
#include <vector>

namespace mu2e { 
  class MedianCalculator{
    struct  MedianData {
      MedianData(float Val, float Wg): val(Val), wg(Wg){}
      float    val;
      float    wg;
    };
    struct MedianDatacomp : public std::binary_function<MedianData,MedianData,bool> {
      bool operator()(MedianData const& p1, MedianData const& p2) { return p1.val < p2.val; }
    };
   
  public:
    MedianCalculator(size_t nToReserve=0){
      _vec.reserve(nToReserve);
    }
    
    float  weightedMedian();
    float  unweightedMedian();
    
    inline void     push(float  value, float   weight=1){
      _vec.emplace_back(MedianData(value, weight));
      _needsSorting = true;
      _totalWeight  += weight;
    }
    
    inline size_t   size(){ return _vec.size(); }
  private:
    
    std::vector<MedianData>  _vec; 
    bool                     _needsSorting     = true;
    float                    _weightedMedian   = 0;
    float                    _unweightedMedian = 0;
    float                    _totalWeight      = 0;
  };
} // namespace mu2e

#endif /* Mu2eUtilities_MedianCalculator_hh */
