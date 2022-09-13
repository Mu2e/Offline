#include "Offline/Mu2eKinKal/inc/StrawXingUpdater.hh"
namespace mu2e {
  StrawXingUpdater::StrawXingUpdater(SXUConfig const& sxusetting) {
      maxdoca_ = std::get<0>(sxusetting);
      maxddoca_ = std::get<1>(sxusetting);
      maxdocasig_ = std::get<2>(sxusetting);
      scalevar_ = std::get<3>(sxusetting);
  }
  std::string const& StrawXingUpdater::configDescription() {
    static std::string descrip( "Maximum DOCA to use straw material, Maximum DOCA to use unaveraged material,Maximum DOCA error to use unaveraged material, scale variance with annealing temp?");
    return descrip;
  }
}
