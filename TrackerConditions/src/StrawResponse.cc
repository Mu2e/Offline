//
// Reconstruction functions for straws
// Original author David Brown, LBNL
//
#include "TrackerConditions/inc/StrawResponse.hh"
// data products
#include "RecoDataProducts/inc/StrawHit.hh"
#include <math.h>
#include <algorithm>

using namespace std;
namespace mu2e {
  StrawResponse::StrawResponse(fhicl::ParameterSet const& pset) :
    _pvParams(pset.get<vector<float> >("PropagationVelocityParams",vector<float>{3.228,127.8,16.88})),
    _wpresParams(pset.get<vector<float> >("WirePositionResolutionParams",vector<float>{2.804,21.23,-24.91})),
    _wbuf(pset.get<float>("WireLengthBuffer",-5.0)) //mm
  {}
  StrawResponse::~StrawResponse(){}
  bool StrawResponse::wireDistance(StrawHit const& strawhit, float shlen, float& wdist, float& wderr) const {
    bool retval(true);
  // convert edep from Mev to KeV (should be standardized, FIXME!)
    float kedep = 1000.0*strawhit.energyDep();
    wdist = halfPropV(kedep)*(strawhit.dt());
    wderr = wpRes(kedep);
    // truncate positions that exceed the length of the straw: these come from missing a cluster on one end
    if(fabs(wdist) > shlen+_wbuf){
      wdist = copysign(shlen+_wbuf,wdist);
      retval = false;
      // does the error change in this case: FIXME!!
    }
    return retval;
  }

  float StrawResponse::halfPropV(float kedep) const {
    float hvp = min(_pvParams[1],_pvParams[1]+(kedep-_pvParams[0])*_pvParams[2]);
    return hvp;
  }

  float StrawResponse::wpRes(float kedep) const {
    float wpr = max(_wpresParams[1],_wpresParams[1]+(kedep-_wpresParams[0])*_wpresParams[2]);
    return wpr;
  }

}
