//
// Instantiation of RegrowKalSeed for LoopHelix fits
//
// Original author: D. Brown (LBNL) 4/18/2025
#include "Offline/Mu2eKinKal/inc/RegrowKalSeed_module.hh"
namespace mu2e {
  using RegrowLoopHelix = RegrowKalSeed<KinKal::LoopHelix>;
}
DEFINE_ART_MODULE(mu2e::RegrowLoopHelix)
