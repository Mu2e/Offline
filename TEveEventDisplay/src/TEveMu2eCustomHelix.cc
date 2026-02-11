#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCustomHelix.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

using namespace mu2e;
namespace mu2e{

  TEveMu2eCustomHelix::TEveMu2eCustomHelix(){}

  /*------------Function to build Infor after contruction:-------------*/
  void TEveMu2eCustomHelix::SetSeedInfo(KalSeed seed) {
    fKalSeed_ = seed;
    auto const& fseg = fKalSeed_.segments().front();
    this->Momentum_ = fseg.mom();
    this->PDGcode_ = fKalSeed_.particle();
    auto const& ptable = mu2e::GlobalConstantsHandle<mu2e::ParticleDataList>();
    this->Charge_ = ptable->particle(fKalSeed_.particle()).charge();
    this->Mass_ = ptable->particle(fKalSeed_.particle()).mass();
    this->Radius_ = fabs(1.0/fseg.helix().omega());
    this->Time_ = fKalSeed_.t0().t0();
  }

}
