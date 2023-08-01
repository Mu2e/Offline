#include "Offline/KinKalGeom/inc/DetectorSolenoid.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    // currently use hard-coded geometry.  So far define only surfaces, volumes and materials will come later
    DetectorSolenoid::DetectorSolenoid() :
      outercyl_(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-1482),950,5450),
      innercyl_(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-1482),1328,5450), // bounding surfaces
      ipa_(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-2770),300.0,500.0), // inner surface
      opa_(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-3766),454.0,728.4,2125.0), // inner surface
      tsda_(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-5967),235.0,525.0) // 25.4mm half thickness
    { }
  }
}
