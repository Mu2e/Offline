#ifndef TrackerConditions_StrawResponse_hh
#define TrackerConditions_StrawResponse_hh
//
// StrawResponse collects the net response features of straws used in reconstruction 
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <vector>
// Mu2e includes
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {
  class StrawHit;
  class StrawResponse : virtual public ConditionsEntity {
    public:
      // construct from parameters
      explicit StrawResponse(fhicl::ParameterSet const& pset);
      virtual ~StrawResponse();
      bool wireDistance(StrawHit const& strawhit, float shlen, float& wdist, float& wderr) const;
      void print(std::ostream& os) const;
    private:
      float halfPropV(float kedep) const;
      float wpRes(float kedep) const;
      // parametric data for calibration functions
      // TD reconstruction uses 1/2 the propagation velocity and depends on the
      // Dependence on position and straw length still needed FIXME!
      // (reconstructed) energy deposit
      std::vector<float> _pvParams;
      // TD resolution, vs energy deposit.
      std::vector<float> _wpresParams;
      float _wbuf; // buffer at the edge of the straws
  };
}
#endif

