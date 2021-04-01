#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCustomHelix.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

using namespace mu2e;
namespace mu2e{

  TEveMu2eCustomHelix::TEveMu2eCustomHelix(){};
 
  TEveMu2eCustomHelix::TEveMu2eCustomHelix(HelixSeed hseed){fHelixSeed = hseed;};
  TEveMu2eCustomHelix::TEveMu2eCustomHelix(KalSeed kseed){
    fKalSeed = kseed;
    this->Momentum = fKalSeed.helix()->helix().momentum();
    this->PDGcode = fKalSeed.particle();
    auto const& ptable = mu2e::GlobalConstantsHandle<mu2e::ParticleDataTable>();
    this->Charge = ptable->particle(fKalSeed.particle()).ref().charge();
    this->Mass = ptable->particle(fKalSeed.particle()).ref().mass().value();
    this->Radius = fKalSeed.helix()->helix().radius();
    this->Time = fKalSeed.t0().t0();
  };

  void TEveMu2eCustomHelix::SetSeedInfo(KalSeed seed) { 
    fKalSeed = seed;
    this->Momentum = fKalSeed.helix()->helix().momentum();
    this->PDGcode = fKalSeed.particle();
    auto const& ptable = mu2e::GlobalConstantsHandle<mu2e::ParticleDataTable>();
    this->Charge = ptable->particle(fKalSeed.particle()).ref().charge();
    this->Mass = ptable->particle(fKalSeed.particle()).ref().mass().value();
    this->Time = fKalSeed.t0().t0();
    this->Radius = fKalSeed.helix()->helix().radius();
  }

  void TEveMu2eCustomHelix::SetPostionAndDirectionFromHelixSeed(double zpos){
    fHelixSeed.helix().position(Position);
    fHelixSeed.helix().direction(zpos, Direction);
  }

  void TEveMu2eCustomHelix::SetPostionAndDirectionFromKalRep(double zpos){
    fKalSeed.helix()->helix().position(Position);
    fKalSeed.helix()->helix().direction(zpos, Direction);
  }
  
  void TEveMu2eCustomHelix::SetMomentumExt(){
    this->Momentum = fTrkExtTraj.front().momentum().mag();
  }

  void TEveMu2eCustomHelix::SetParticleExt(){
    this->PDGcode = 11; 
 
  }
  
}
