#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  std::vector<std::string> WireHitState::StateNames_ = {"Unusable", "Inactive", "Left", "Right", "Null"};
  std::vector<std::string> WireHitState::TOTUseNames_ = { "Unused", "NullOnly", "DriftOny", "All"};
  std::vector<std::string> WireHitState::NullDistVarNames_ = { "StrawRadius", "DriftRadius" };
  std::ostream& operator <<(std::ostream& ost, WireHitState const& whs) {
    ost << "WireHitState ";
    if(whs.frozen()) ost << " Frozen ";
    ost << WireHitState::StateNames_[whs.state_+3];
    ost << " TOTUse "  << WireHitState::TOTUseNames_[whs.totuse_];
    ost << " Null Dist Var " << WireHitState::NullDistVarNames_[whs.nulldvar_];
    return ost;
  }

  bool WireHitState::isIn(WHSMask const& whsmask) const {
    switch (state_) {
      case WireHitState::inactive:
        return whsmask.hasAnyProperty(WHSMask::inactive);
      case WireHitState::left: case WireHitState::right:
        return whsmask.hasAnyProperty(WHSMask::drift);
      case WireHitState::null:
        return whsmask.hasAnyProperty(WHSMask::null);
      default:
        return false;
    }
  }

  WireHitState::TOTUse WireHitState::totUse(std::string const& totuse){
    for (size_t index=0;index < TOTUseNames_.size(); ++index){
      if(TOTUseNames_[index].compare(totuse)){
        return static_cast<TOTUse>(index);
      }
    }
    throw cet::exception("RECO")<<"mu2e::WireHitState: unknown TOTUse " << totuse << std::endl;
  }

  WireHitState::NullDistVar WireHitState::nullDistVar(std::string const& ndvar){
    for (size_t index=0;index < NullDistVarNames_.size(); ++index){
      if(NullDistVarNames_[index].compare(ndvar)){
        return static_cast<NullDistVar>(index);
      }
    }
    throw cet::exception("RECO")<<"mu2e::WireHitState: unknown NullDistVar" << ndvar << std::endl;
  }
}
