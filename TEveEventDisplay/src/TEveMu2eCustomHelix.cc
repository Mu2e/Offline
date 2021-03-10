#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCustomHelix.h"

using namespace mu2e;
namespace mu2e{

  TEveMu2eCustomHelix::TEveMu2eCustomHelix(){};
 
  TEveMu2eCustomHelix::TEveMu2eCustomHelix(HelixSeed hseed){fHelixSeed = hseed;};
  TEveMu2eCustomHelix::TEveMu2eCustomHelix(KalSeed kseed){
    fKalSeed = kseed;
    this->Momentum = fKalSeed.helix()->helix().momentum();
    this->PDGcode = fKalSeed.particle().particleType();
    this->Charge = fKalSeed.particle().charge();
    this->Mass = fKalSeed.particle().mass();
  };

  void TEveMu2eCustomHelix::SetSeedInfo(KalSeed seed) { 
    fKalSeed = seed;
    this->Momentum = fKalSeed.helix()->helix().momentum();
    this->PDGcode = fKalSeed.particle().particleType();
    this->Charge = fKalSeed.particle().charge();
    this->Mass = fKalSeed.particle().mass();
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

  std::string TEveMu2eCustomHelix::DataTitle(const std::string &pstr, int n){
        std::string dstr=" hit#" + std::to_string(n) + "\nLayer: ";
        std::string strlab=pstr+dstr;
        return (strlab);
  }
  
  const std::string TEveMu2eCustomHelix::Title(){
    const std::string title = "Track PDG " + std::to_string(this->PDGcode) +" Momentum = " + std::to_string(this->Momentum) + " Charge "+ std::to_string(this->Charge);
    return title;
  }
  
}
