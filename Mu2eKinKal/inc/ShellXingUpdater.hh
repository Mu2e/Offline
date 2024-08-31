#ifndef Mu2eKinKal_ShellXingUpdater_hh
#define Mu2eKinKal_ShellXingUpdater_hh
#include <tuple>
namespace mu2e {
  // simple struct to hold smalle material Xing calculation configuration parameters
  struct ShellXingUpdater {
    using SXUConfig = std::tuple<bool,int>;
    static std::string const& configDescription(); // description of the variables
    bool scalevar_ =false; // scale variance or not
    int diag_ =0; // diag print level
    // default constructor is functional but will always use the impact-parameter averaged material
    ShellXingUpdater(SXUConfig const& sxusetting);
  };
}
#endif
