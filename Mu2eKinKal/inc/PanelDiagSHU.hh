#ifndef Mu2eKinKal_PanelDiagSHU_hh
#define Mu2eKinKal_PanelDiagSHU_hh
//
// Diagnostic updater of StrawHits that can force panels to be inactive
//
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include <tuple>
#include <string>
#include <iostream>

namespace mu2e {
  // Update based just on PTCA to the wire
  class PanelDiagSHU {
    public:
      using Config = std::tuple<std::string,std::string,int>;
      PanelDiagSHU(Config const& config);
      static std::string const& configDescription(); // description of the variables
      // set the state based on the current PTCA value
      WireHitState wireHitState(WireHitState const& input, StrawId const& id) const;
    private:
      StrawIdMask mask_;
      StrawId sid_;
      int diag_ =0; // diag print level
  };
}
#endif
