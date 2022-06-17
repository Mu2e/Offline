#include "Offline/Mu2eKinKal/inc/WHSIterator.hh"
namespace mu2e {
  bool WHSIterator::increment() {
    while(ihit_ < nhits_){
      ++indices_[ihit_];
      if(indices_[ihit_] < allowed_.size()){
        break;
      } else {
        for(size_t jhit = 0; jhit <= ihit_; ++jhit)indices_[jhit] = 0;
        ++ihit_;
      }
    }
    if(ihit_ < nhits_){
      for(size_t jhit=0; jhit < nhits_; ++jhit){
        current_[jhit] = allowed_[indices_[jhit]];
      }
      return true;
    } else {
      return false;
    }
  }
  void WHSIterator::reset() {
    ihit_ = 0;
    for(size_t ihit =0; ihit < nhits_; ++ihit)indices_[ihit] = 0;
  }
}
