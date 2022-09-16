#include "Offline/Mu2eKinKal/inc/StrawXingUpdater.hh"
#include <iostream>
namespace mu2e {
  StrawXingUpdater::StrawXingUpdater(SXUConfig const& sxusetting) {
    maxdoca_ = std::get<0>(sxusetting);
    maxddoca_ = std::get<1>(sxusetting);
    maxdocasig_ = std::get<2>(sxusetting);
    scalevar_ = std::get<3>(sxusetting);
    diag_ = std::get<4>(sxusetting);
    if(diag_ > 0) std::cout << "StrawXingUpdater maxdoca " << maxdoca_ << " max unaveraged doca, doca error " << maxddoca_ << ", " << maxdocasig_ << " scale variance? " << scalevar_ << std::endl;
  }
  std::string const& StrawXingUpdater::configDescription() {
    static std::string descrip( "Maximum DOCA to use straw material, Maximum DOCA to use unaveraged material,Maximum DOCA error to use unaveraged material, scale variance with annealing temp?");
    return descrip;
  }
}
