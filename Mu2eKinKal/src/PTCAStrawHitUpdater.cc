#include "Offline/Mu2eKinKal/inc/PTCAStrawHitUpdater.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  WireHitState PTCAStrawHitUpdater::wireHitState(ClosestApproachData const& tpdata, Straw const& straw ) const {
    WireHitState whstate(WireHitState::inactive);
    double doca = tpdata.doca();
    double absdoca = fabs(doca);
    if( absdoca < maxdoca_ && absdoca > mindoca_
        && tpdata.deltaT() < maxdt_ && tpdata.deltaT() > mindt_){  // in the sweet spot: use the DOCA to sign the ambiguity
      whstate = doca > 0.0 ? WireHitState::right : WireHitState::left;
    } else if(absdoca < mindoca_ || tpdata.deltaT() < mindt_) {
      whstate = WireHitState::null;
    }
    return whstate;
  }

  NullHitInfo PTCAStrawHitUpdater::nullHitInfo(StrawResponse const& sresponse, Straw const& straw) const {
    NullHitInfo nhinfo;
    // compute time and distance parameters used for null ambiguity (wire constraint)
    double vdriftinst = sresponse.driftInstantSpeed(straw.id(),mindoca_,0.0,true);
    nhinfo.toff_ = 0.5*mindoca_/vdriftinst;
    nhinfo.tvar_ = 0.25*dvar_/(vdriftinst*vdriftinst);
    nhinfo.dvar_ = dvar_;
    nhinfo.usetime_ = nulltime_;
    return nhinfo;
  }

  bool PTCAStrawHitUpdater::insideStraw(KinKal::ClosestApproachData const& ca,Straw const& straw) const {
    static const double ubuffer(10.0); // should be a parameter FIXME
    // compute the position along the wire and compare to the 1/2 length
    // have to translate from CLHEP, should be native to Straw FIXME
    double upos = VEC3(straw.wireDirection()).Dot((ca.sensorPoca().Vect() - VEC3(straw.origin())));
    return fabs(upos) < straw.halfLength() + ubuffer;
  }

}
