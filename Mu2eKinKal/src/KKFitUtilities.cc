#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
namespace mu2e {
  namespace Mu2eKinKal {
    KinKal::Line hitLine(ComboHit const& ch, Straw const& straw,StrawResponse const& strawresponse) {
      double sprop = 2*strawresponse.halfPropV(ch.strawId(),ch.energyDep());
      // construct a kinematic line trajectory from this straw. the measurement point is the signal end
      KinKal::VEC3 vp0(straw.wireEnd(ch.driftEnd()));
      KinKal::VEC3 vp1(straw.wireEnd(StrawEnd(ch.driftEnd().otherEnd())));
      return KinKal::Line(vp0,vp1,ch.time(),sprop);
    }
  }
}
