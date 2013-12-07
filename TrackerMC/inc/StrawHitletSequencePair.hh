#ifndef TrackerMC_StrawHitletSequencePair_hh
#define TrackerMC_StrawHitletSequencePair_hh
// Define HitletSequencePair, ie the hitlet sequence from both ends of a straw
//
// $Id: StrawHitletSequencePair.hh,v 1.1 2013/12/07 19:51:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:51:42 $
//
// Original author David Brown, LBNL
#include "TrackerMC/inc/StrawHitletSequence.hh"
#include <utility>
namespace mu2e {
  class StrawHitletSequencePair{
    public:
      typedef StrawHitlet StrawHitletPair[2];
      StrawHitletSequencePair();
      StrawHitletSequencePair(StrawIndex index);
      StrawHitletSequencePair(StrawHitletSequencePair const& other);
      StrawHitletSequencePair& operator =(StrawHitletSequencePair const& other);
      StrawHitletSequence& hitletSequence(StrawEnd end) { if(end==StrawEnd::plus) return _plus; else return _minus; }
      StrawHitletSequence const& hitletSequence(StrawEnd end) const { if(end==StrawEnd::plus) return _plus; else return _minus; }
      void insert(StrawHitletPair const& hpair);
      StrawIndex strawIndex() const { return _plus.strawIndex(); }
    private:
      StrawHitletSequence _plus, _minus; // plus and minus end sequences
  };
}
#endif
