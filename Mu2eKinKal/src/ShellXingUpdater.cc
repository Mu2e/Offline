#include "Offline/Mu2eKinKal/inc/ShellXingUpdater.hh"
#include <iostream>
namespace mu2e {
  ShellXingUpdater::ShellXingUpdater(SXUConfig const& sxusetting) {
    scalevar_ = std::get<0>(sxusetting);
    diag_ = std::get<1>(sxusetting);
    if(diag_ > 0) std::cout << "ShellXingUpdater ScaleVariance?" << scalevar_ << std::endl;
  }
  std::string const& ShellXingUpdater::configDescription() {
    static std::string descrip( "scale variance with annealing temp?; Diag level");
    return descrip;
  }
}
