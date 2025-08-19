#ifndef Mu2eKinKal_StrawXingUpdater_hh
#define Mu2eKinKal_StrawXingUpdater_hh
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include <tuple>
namespace mu2e {
  // simple struct to hold crossing calculation configuration parameters
  struct StrawXingUpdater {
    using SXUConfig = std::tuple<float,float,float,bool,int>;
    static std::string const& configDescription(); // description of the variables
    double maxdoca_ =0; // maximum DOCA to include this straw's material
    double maxddoca_ =0; // maximum DOCA to use non-averaged value
    double maxdocasig_ =0; // maximum doca error to use non-averaged value
    bool scalevar_ =false; // scale variance or not
    int diag_ =0; // diag print level
    // default constructor is functional but will always use the impact-parameter averaged material
    StrawXingUpdater(SXUConfig const& sxusetting);
    StrawXingUpdater(){}
  };
}
#endif
