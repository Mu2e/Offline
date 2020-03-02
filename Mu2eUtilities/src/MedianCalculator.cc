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
  
    float   sum(0);
    size_t  id(0);

    sum = _totalWeight - _vec[0].wg;
    while (sum > 0.5*_totalWeight){
      ++id;
      sum -= _vec[id].wg;
    }

    float   over((sum)/_totalWeight);
    float   interpolation(0);
    if (v_size %2 == 0) {
      interpolation =  _vec[id].val * over + _vec[id+1].val * (1.-over);
    }else {
      float  w2     = (sum)/_totalWeight;
      float  w1     = (sum + _vec[id].wg )/_totalWeight;
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
    size_t  id(0);

    float   interpolation(0);
    if (v_size %2 == 0) {
      id = v_size/2 - 1;
      interpolation =  _vec[id].val * 0.5 + _vec[id+1].val * 0.5;
    }else {
      id = v_size/2;
      float  sum(id);      
      float  w2     = (sum)/totWg;
      float  w1     = (sum + 1.)/totWg;
      float  val1   = _vec[id-1].val*w1 + _vec[id].val  *(1.-w1);
      float  val2   = _vec[id].val  *w2 + _vec[id+1].val*(1.-w2);
      interpolation = 0.5*(val1 + val2);
    }

    //cache the result
    _weightedMedian = interpolation;

    return interpolation;
  }
  
}
