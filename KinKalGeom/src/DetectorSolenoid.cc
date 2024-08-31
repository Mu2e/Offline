#include "Offline/KinKalGeom/inc/DetectorSolenoid.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Cylinder;
    using KinKal::Disk;
    using KinKal::Frustrum;
    using KinKal::Annulus;

    // currently use hard-coded geometry.  So far define only surfaces, volumes and materials will come later
    DetectorSolenoid::DetectorSolenoid() :
      outer_{ std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-1482),950,5450)},
      inner_{ std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-1482),1328,5450)}, // bounding surfaces
      front_{ std::make_shared<Disk>(outer_->frontDisk())},
      back_{ std::make_shared<Disk>(outer_->backDisk())},
      ipa_{ std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-2770),300.0,500.0)},
      ipa_front_{ std::make_shared<Disk>(ipa_->frontDisk())},
      ipa_back_{ std::make_shared<Disk>(ipa_->backDisk())},
      opa_{ std::make_shared<Frustrum>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-3766),454.0,728.4,2125.0)}, // inner surface
      tsda_{ std::make_shared<Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-5967),235.0,525.0)} // back surface
    { }
  }
}
