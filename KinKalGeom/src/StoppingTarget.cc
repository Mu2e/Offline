#include "Offline/KinKalGeom/inc/StoppingTarget.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Cylinder;
    using KinKal::Disk;
    using KinKal::Annulus;
    // currently use hard-coded geometry
    StoppingTarget::StoppingTarget() :
      outer_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-4300),75,400.0)},
      inner_{ std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-4300),21.5,400.0)},
      front_{ std::make_shared<Disk>(outer_->frontDisk())},
      back_{ std::make_shared<Disk>(outer_->backDisk())} {
        double startz = -4700;
        double endz = -3900;
        double dz = (endz-startz)/36.0;
        for(int ifoil=0;ifoil < 37; ++ifoil){
          double zpos = startz + ifoil*dz;
          foils_.push_back(std::make_shared<Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zpos),21.5,75));
        }
      }
  }
}
