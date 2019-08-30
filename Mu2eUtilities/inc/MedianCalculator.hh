#ifndef Mu2eUtilities_MedianCalculator_hh
#define Mu2eUtilities_MedianCalculator_hh
//
//
// $Id: $
// $Author: $
// $Date: $
//
// Original author G. Pezzullo
//

#include <vector>
#include <utility>
#include <math.h>
#include <cmath>
#include <numeric>

namespace mu2e {
  struct  MedianData {
      MedianData(float Val, float Wg){
	val = Val;
	wg  = Wg; 
      }
      float    val;
      float    wg;
    };
    struct MedianDatacomp : public std::binary_function<MedianData,MedianData,bool> {
      bool operator()(MedianData const& p1, MedianData const& p2) { return p1.val < p2.val; }
    };
    
  class MedianCalculator{
    
  public:
    MedianCalculator(){}
    
    float  extractMedian(std::vector<MedianData> &v){
      std::sort(v.begin(), v.end(), MedianDatacomp());
      //now, we need to loop over it and evaluate the median
      size_t   v_size = v.size();
      if (v_size ==0) {
	return -1e5;
      }
      if (v_size == 1) return v[0].val;

      float   totWg(std::accumulate(v.begin(), v.end(), 0, [](float sum, const MedianData& curr){return sum + curr.wg;}));
      float   sum(0);
      size_t  id(0);

      sum = totWg - v[0].wg;
      while (sum > 0.5*totWg){
	++id;
	sum -= v[id].wg;
      }

      float   over((sum)/totWg);
      float   interpolation(0);
      if (v_size %2 == 0) {
	interpolation =  v[id].val * over + v[id+1].val * (1.-over);
      }else {
	float  w2     = (sum)/totWg;
	float  w1     = (sum + v[id].wg )/totWg;
	float  val1   = v[id-1].val*w1 + v[id].val*(1.-w1);
	float  val2   = v[id].val*w2 + v[id+1].val*(1.-w2);
	interpolation = 0.5*(val1 + val2);
      }
      return interpolation;
    }
  };
} // namespace mu2e

#endif /* Mu2eUtilities_MedianCalculator_hh */
