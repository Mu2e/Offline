#ifndef Mu2eUtilities_MedianCalculator_hh
#define Mu2eUtilities_MedianCalculator_hh
//
// Original author G. Pezzullo
//
// this class is intended to be used for evaluaitng the median 
// from a set of elements that are stored internally in a vector
//

#include <vector>
//#include <utility>
#include <numeric>
#include <functional>

namespace mu2e { 
  struct  MedianData {
    MedianData(float Val, float Wg){
      val = Val;
      wg  = Wg; 
    }
    float    val;
    float    wg;
  };
  struct MedianDatacomp : public std::binary_function<mu2e::MedianData,mu2e::MedianData,bool> {
    bool operator()(MedianData const& p1, MedianData const& p2) { return p1.val < p2.val; }
  };

  class MedianCalculator{
    
  public:
    MedianCalculator():_needsSorting(true){}
    
    float  weightedMedian();
    float  unweightedMedian();
    
    void   push(float  value, float   weight=1.);

  private:
    
    std::vector<MedianData>  _vec; 
    bool                     _needsSorting;
    float                    _weightedMedian;
    float                    _unweightedMedian;
    
  };
} // namespace mu2e

#endif /* Mu2eUtilities_MedianCalculator_hh */
