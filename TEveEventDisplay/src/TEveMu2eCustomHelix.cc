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

  /*------------Function tobuild position and direction based on Kal output:-------------*/
  void TEveMu2eCustomHelix::SetPostionAndDirectionFromKalRep(double zpos){
    auto const& fseg = fKalSeed_.segments().front();  // find the segment nearest zpos.  This should be iterative FIXME
    auto pos = fseg.position3();
    auto vel = fseg.velocity();
    double tz = fseg.tref() + (zpos-pos.Z())/vel.Z();
    auto zseg = fKalSeed_.nearestSegment(tz);
    pos = zseg->position3();
    vel = zseg->velocity();
    tz = zseg->tref() + (zpos-pos.Z())/vel.Z();
    this->Momentum_ = zseg->mom();
    // these next assume a helix FIXME
    auto hel = zseg->centralHelix();
    Position_ = hel.position3(tz);
    Direction_ = hel.direction(tz);

  }

  /*void TEveMu2eCustomHelix::SetMomentumExt(){
    this->Momentum_ = fTrkExtTraj.front().momentum().mag(); // Not sure what this is supposed to do FIXME
  }

  void TEveMu2eCustomHelix::SetParticleExt(){
    this->PDGcode_ = 11;

  }*/

}
