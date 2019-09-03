#include "Mu2eUtilities/inc/MedianCalculator.hh"
#include <iostream>
#include "cetlib_except/exception.h"

namespace mu2e {

  float    MedianCalculator::weightedMedian(){
    //now, we need to loop over it and evaluate the median
    size_t   v_size = _vec.size();
    if (v_size ==0) {
      throw cet::exception("MATH")<<"No entries in the vector: median undefined" << std::endl;
    }
    if (v_size == 1){
      return _vec[0].val;
    }
    
    if (_needsSorting){
      std::sort(_vec.begin(), _vec.end(), MedianDatacomp());
      _needsSorting = false;
    }else {
      return   _weightedMedian;
    }
  
    float   totWg(std::accumulate(_vec.begin(), _vec.end(), 0, [](float sum, const MedianData& curr){return sum + curr.wg;}));
    float   sum(0);
    size_t  id(0);

    sum = totWg - _vec[0].wg;
    while (sum > 0.5*totWg){
      ++id;
      sum -= _vec[id].wg;
    }

    float   over((sum)/totWg);
    float   interpolation(0);
    if (v_size %2 == 0) {
      interpolation =  _vec[id].val * over + _vec[id+1].val * (1.-over);
    }else {
      float  w2     = (sum)/totWg;
      float  w1     = (sum + _vec[id].wg )/totWg;
      float  val1   = _vec[id-1].val*w1 + _vec[id].val*(1.-w1);
      float  val2   = _vec[id].val*w2 + _vec[id+1].val*(1.-w2);
      interpolation = 0.5*(val1 + val2);
    }

    //cache the result
    _weightedMedian = interpolation;

    return interpolation;
  }
	
  float    MedianCalculator::unweightedMedian(){
    //now, we need to loop over it and evaluate the median
    size_t   v_size = _vec.size();
    if (v_size ==0) {
      throw cet::exception("MATH")<<"No entries in the vector: median undefined" << std::endl;
    }
    if (v_size == 1){
      return _vec[0].val;
    }

    if (_needsSorting){
      std::sort(_vec.begin(), _vec.end(), MedianDatacomp());
      _needsSorting = false;
    }else {
      return   _unweightedMedian;
    }
 
    float   totWg(_vec.size());
    float   sum(0);
    size_t  id(0);

    sum = totWg - 1.;//_vec[0].wg;
    while (sum > 0.5*totWg){
      ++id;
      sum -= 1.;//_vec[id].wg;
    }

    float   over((sum)/totWg);
    float   interpolation(0);
    if (v_size %2 == 0) {
      interpolation =  _vec[id].val * over + _vec[id+1].val * (1.-over);
    }else {
      float  w2     = (sum)/totWg;
      float  w1     = (sum + 1.)/totWg;
      float  val1   = _vec[id-1].val*w1 + _vec[id].val*(1.-w1);
      float  val2   = _vec[id].val*w2 + _vec[id+1].val*(1.-w2);
      interpolation = 0.5*(val1 + val2);
    }

    //cache the result
    _weightedMedian = interpolation;

    return interpolation;
  }
  
  void     MedianCalculator::push(float  value, float   weight){
    _vec.push_back(MedianData(value, weight));
    _needsSorting = true;
  }

}
