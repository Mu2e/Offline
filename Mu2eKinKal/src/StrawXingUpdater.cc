#include "Offline/Mu2eKinKal/inc/StrawXingUpdater.hh"
#include <iostream>
namespace mu2e {
  StrawXingUpdater::StrawXingUpdater(SXUConfig const& sxusetting) {
    maxdoca_ = std::get<0>(sxusetting);
    nsig_ = std::get<1>(sxusetting);
    scalevar_ = std::get<2>(sxusetting);
    diag_ = std::get<3>(sxusetting);
    if(diag_ > 0) std::cout << "StrawXingUpdater maxdoca " << maxdoca_ << " averaging factor" << nsig_ << " scale variance? " << scalevar_ << std::endl;
  }
  std::string const& StrawXingUpdater::configDescription() {
    static std::string descrip( "Maximum DOCA to use straw material, Averaging factor (WRT DOCA error), scale variance with annealing temp?, Diag level");
    return descrip;
  }
}
