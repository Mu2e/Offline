#ifndef TrackerMC_StrawHitletSequence_hh
#define TrackerMC_StrawHitletSequence_hh
//
// StrawHitletSequence is a time-ordered sequence of StrawHitlets
//
// $Id: StrawHitletSequence.hh,v 1.1 2013/12/07 19:51:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:51:42 $
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <list>
// Mu2e includes
#include "TrackerMC/inc/StrawHitlet.hh"
#include "DataProducts/inc/StrawIndex.hh"

namespace mu2e {
  typedef std::list<StrawHitlet> HitletList;
  class StrawHitletSequence {
    public:
// constructors
      StrawHitletSequence();
      StrawHitletSequence(StrawHitlet const& hitlet);
      StrawHitletSequence(StrawIndex const& index, StrawEnd end);
      StrawHitletSequence(StrawHitletSequence const& other);
      StrawHitletSequence& operator =(StrawHitletSequence const& other);
      // accessors: just hand over the list!
      HitletList const& hitletList() const { return _hlist; }
      // insert a new hitlet, in time order.
      HitletList::iterator insert(StrawHitlet const& hitlet);
      StrawIndex strawIndex() const { return _strawIndex; }
      StrawEnd strawEnd() const { return _end; }
    private:
      StrawIndex _strawIndex;
      StrawEnd _end;
      HitletList _hlist; // time-ordered sequence of hitlets
  };
}
#endif


