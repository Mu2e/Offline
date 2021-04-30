#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCustomHelix.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

using namespace mu2e;
namespace mu2e{

  TEveMu2eCustomHelix::TEveMu2eCustomHelix(){};
 
  TEveMu2eCustomHelix::TEveMu2eCustomHelix(HelixSeed hseed){fHelixSeed = hseed;};
  TEveMu2eCustomHelix::TEveMu2eCustomHelix(KalSeed kseed){
    fKalSeed = kseed;  // why is a deep copy being made?  FIXME
    auto const& fseg = fKalSeed.segments().front(); // this is the segment at the earliest point, impossible to know if it's what the user wants FIXME!
    this->Momentum = fseg.mom();
    this->PDGcode = fKalSeed.particle();
    auto const& ptable = mu2e::GlobalConstantsHandle<mu2e::ParticleDataTable>();
    this->Charge = ptable->particle(fKalSeed.particle()).ref().charge();
    this->Mass = ptable->particle(fKalSeed.particle()).ref().mass().value();
    // the following uses deprecated legacy functions and should be refactored away FIXME
    this->Radius = fabs(1.0/fseg.helix().omega()); // radius only makes sense if there is a magnetic field, this implementation restricts the display to helices FIXME
    this->Time = fKalSeed.t0().t0();
  };

  void TEveMu2eCustomHelix::SetSeedInfo(KalSeed seed) { 
  // this function is identical to the constructor above, it should be consolidated FIXME
    fKalSeed = seed; // class variables have no underscore, violates Mu2e coding guidelines FIXME
    auto const& fseg = fKalSeed.segments().front();
    this->Momentum = fseg.mom();
    this->PDGcode = fKalSeed.particle();
    auto const& ptable = mu2e::GlobalConstantsHandle<mu2e::ParticleDataTable>();
    this->Charge = ptable->particle(fKalSeed.particle()).ref().charge();
    this->Mass = ptable->particle(fKalSeed.particle()).ref().mass().value();
    this->Radius = fabs(1.0/fseg.helix().omega());
    this->Time = fKalSeed.t0().t0();
  }

  void TEveMu2eCustomHelix::SetPostionAndDirectionFromHelixSeed(double zpos){// how is this function used?  Does the user really want helix seed info?  FIXME!
    fHelixSeed.helix().position(Position);
    fHelixSeed.helix().direction(zpos, Direction);
  }

  void TEveMu2eCustomHelix::SetPostionAndDirectionFromKalRep(double zpos){
    auto const& fseg = fKalSeed.segments().front();  // find the segment nearest zpos.  This should be iterative FIXME
    auto pos = fseg.position3();
    auto vel = fseg.velocity();
    double tz = fseg.tref() + (zpos-pos.Z())/vel.Z();
    auto zseg = fKalSeed.nearestSeg(tz);
    pos = zseg->position3();
    vel = zseg->velocity();
    tz = zseg->tref() + (zpos-pos.Z())/vel.Z();
    this->Momentum = zseg->mom();
    // these next assume a helix FIXME
    auto hel = zseg->centralHelix();
    Position = hel.position3(tz); // class variables have no underscore, violates Mu2e coding guidelines FIXME
    Direction = hel.direction(tz); // class variables have no underscore, violates Mu2e coding guidelines FIXME

  }
  
  void TEveMu2eCustomHelix::SetMomentumExt(){
    this->Momentum = fTrkExtTraj.front().momentum().mag(); // Not sure what this is supposed to do FIXME
  }

  void TEveMu2eCustomHelix::SetParticleExt(){
    this->PDGcode = 11; 
 
  }
  
}
