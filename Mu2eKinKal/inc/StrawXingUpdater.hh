#ifndef Mu2eKinKal_StrawXingUpdater_hh
#define Mu2eKinKal_StrawXingUpdater_hh
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include <tuple>
namespace mu2e {
  // simple struct to hold crossing calculation configuration parameters
  struct StrawXingUpdater {
    using SXUConfig = std::tuple<float,float,bool,int>;
    static std::string const& configDescription(); // description of the variables
    double maxdoca_ =0; // maximum DOCA to activate straw materials without an associated hits. <0 means don't change the state
    double nsig_ =0; // Number of doca_sigma around doca value to use when averageing
    bool scalevar_ =false; // scale variance or not
    int diag_ =0; // diag print level
    // default constructor is functional but will always use the impact-parameter averaged material
    StrawXingUpdater(SXUConfig const& sxusetting);
    StrawXingUpdater(){}
  };
}
#endif
