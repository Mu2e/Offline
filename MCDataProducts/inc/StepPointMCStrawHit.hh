#ifndef MCDataProducts_StepPointMCStrawHit_hh
#define MCDataProducts_StepPointMCStrawHit_hh
//
// $Id: StepPointMCStrawHit.hh,v 1.1 2012/08/22 22:19:42 genser Exp $
// $Author: genser $
// $Date: 2012/08/22 22:19:42 $
//
// Original author KLG based on Rob's MakeStrawHit_module StrawHit

// Mu2e includes
#include "MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {

  // Utility class (structure) to hold calculated drift time for G4 hits

  class StepPointMCStrawHit {

    // name it StepPointMCStrawHit?

  public:

    art::Ptr<StepPointMC> _ptr;
    double _edep;
    double _dca;
    double _driftTimeNonSm;
    double _driftTime;
    double _distanceToMid;
    double _t1;
    double _t2;

    StepPointMCStrawHit() {}
    StepPointMCStrawHit(art::Ptr<StepPointMC> const& ptr, 
            double edep, 
            double dca, 
            double driftTNonSm, 
            double driftT, 
            double toMid, 
            double t1, 
            double t2):
      _ptr(ptr), 
      _edep(edep), 
      _dca(dca), 
      _driftTimeNonSm(driftTNonSm),
      _driftTime(driftT),
      _distanceToMid(toMid), 
      _t1(t1), 
      _t2(t2) {}

    // This operator is overloaded in order to time-sort the hits
    bool operator <(const StepPointMCStrawHit& b) const { return (_t1 < b._t1); }

  };

} // namespace mu2e

#endif /* MCDataProducts_StepPointMCStrawHit_hh */
