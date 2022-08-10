#include "Offline/Mu2eKinKal/inc/PTCAStrawHitUpdater.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  WireHitState PTCAStrawHitUpdater::wireHitState(ClosestApproachData const& tpdata, Straw const& straw, StrawResponse const& sresponse ) const {
    WireHitState whstate(WireHitState::inactive,StrawHitUpdaters::PTCA);
    double doca = tpdata.doca();
    double absdoca = fabs(doca);
    if( absdoca < maxdoca_ && absdoca > mindoca_
        && tpdata.deltaT() < maxdt_ && tpdata.deltaT() > mindt_){  // in the sweet spot: use the DOCA to sign the ambiguity
      whstate.state_ = doca > 0.0 ? WireHitState::right : WireHitState::left;
    } else if(absdoca < mindoca_ || tpdata.deltaT() < mindt_) {
      whstate.state_ = WireHitState::null;
    }
    // compute time and distance parameters used for null ambiguity (wire constraint)
    auto& nhinfo = whstate.nhinfo_;
    if(nulltime_ == usecombo){
      nhinfo.usetime_ = true;
      nhinfo.toff_ =  0.0;
      nhinfo.useComboDriftTime_ = true;
      nhinfo.tvar_ = 50.0; // should come from ComboHit FIXME
      nhinfo.dvar_ = 2.1;
    } else if(nulltime_ == usedoca){
      double vdrift = sresponse.driftConstantSpeed(); // use average speed for now
      nhinfo.usetime_ = true;
      nhinfo.toff_ = 0.5*mindoca_/vdrift; // this calculation is unreliable currently
      nhinfo.useComboDriftTime_ = false;
      nhinfo.dvar_ = dvar_;
      nhinfo.tvar_ = dvar_/(vdrift*vdrift);
    } else {
      nhinfo.usetime_ = false;
    }
    return whstate;
  }
}
