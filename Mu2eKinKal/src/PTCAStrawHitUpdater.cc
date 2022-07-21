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
    if(tvar_ <= 0.0){
      // defer evaluation as StrawResponse is needed.  This should be a prodition FIXME
      // scan the drift info and take the average
      static const unsigned nscan(10);
      static const double dstep = mindoca_/(nscan-1);
      toff_ = tvar_ = 0.0;
      for(unsigned iscan=0;iscan<nscan;++iscan){
        double dist=iscan*dstep;
        auto dinfo = sresponse.driftInfoAtDistance(straw.id(),dist,0.0);
        toff_ += dinfo.time;
        tvar_ += dinfo.variance;
      }
      toff_ /= nscan; tvar_ /= nscan;
    }
    nhinfo.toff_ = toff_;
    nhinfo.tvar_ = tvar_;
    nhinfo.dvar_ = dvar_;
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
