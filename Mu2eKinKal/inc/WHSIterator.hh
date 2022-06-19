#ifndef Mu2eKinKal_WHSIteratorator_hh
#define Mu2eKinKal_WHSIteratorator_hh
#include "KinKal/Detector/WireHitStructs.hh"
#include <vector>
#include <cmath>
//
//  Iterator across all allowed states of a cluster of StrawHits
//
namespace mu2e {
  using KinKal::WireHitState;
  using WHSCOL = std::vector<WireHitState>;
  class WHSIterator {
    public:
      WHSIterator(size_t nhits, WHSCOL const& allowed) : nhits_(nhits),allowed_(allowed), current_(nhits, allowed_.front()),
      indices_(4,0){}
      size_t nCombo() const { return static_cast<size_t>(std::rint(std::pow(allowed_.size(),nhits_))); }
      bool increment();
      void reset();
      auto const& current() const { return current_; }
    private:
      size_t nhits_; // number of hits
      WHSCOL const& allowed_; // allowed states for each hit
      WHSCOL current_; // current state of all hits
      std::vector<size_t> indices_; // current indices into allowed states for each hit
  };
}
#endif
