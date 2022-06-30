#include "Offline/Mu2eKinKal/inc/WHSIterator.hh"
namespace mu2e {
  bool WHSIterator::increment() {
    size_t ihit(0);
    while(ihit < nhits_){
      ++indices_[ihit];
      if(indices_[ihit] < allowed_.size()){
        break;
      } else {
        // flip to a higher hit, and zero all the lower indices
        for(size_t jhit = 0; jhit <= ihit; ++jhit)indices_[jhit] = 0;
        ++ihit;
      }
    }
    if(ihit < nhits_){
      for(size_t jhit=0; jhit < nhits_; ++jhit){
        // reverse the order, so we iterate higher-index hits first
        size_t khit = nhits_-jhit-1;
        current_[khit] = allowed_[indices_[jhit]];
      }
      return true;
    } else {
      return false;
    }
  }
  void WHSIterator::reset() {
    for(auto& index : indices_)index = 0;
  }
}
