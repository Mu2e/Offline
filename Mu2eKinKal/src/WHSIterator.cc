#include "Offline/Mu2eKinKal/inc/WHSIterator.hh"
namespace mu2e {
  bool WHSIterator::increment() {
    size_t ihit(0);
    while(ihit < nhits_){
      ++indices_[ihit];
      if(indices_[ihit] < allowed_.size()){
        break;
      } else {
        for(size_t jhit = 0; jhit <= ihit; ++jhit)indices_[jhit] = 0;
        ++ihit;
      }
    }
    if(ihit < nhits_){
      for(size_t jhit=0; jhit < nhits_; ++jhit){
        current_[jhit] = allowed_[indices_[jhit]];
      }
      return true;
    } else {
      return false;
    }
  }
  void WHSIterator::reset() {
    for(size_t ihit =0; ihit < nhits_; ++ihit)indices_[ihit] = 0;
  }
}
