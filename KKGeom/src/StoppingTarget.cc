#include "Offline/KKGeom/inc/StoppingTarget.hh"
namespace mu2e {
  namespace KKGeom {
    using KinKal::VEC3;
    using KinKal::Annulus;
    // currently use hard-coded geometry
    StoppingTarget::StoppingTarget() :
      outercyl_(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-4300),75,400.0),
      innercyl_(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-4300),21.5,400.0) {
        double startz = -4700;
        double endz = -3900;
        double dz = (endz-startz)/36.0;
        for(int ifoil=0;ifoil < 37; ++ifoil){
          double zpos = startz + ifoil*dz;
          foils_.push_back(Annulus(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zpos),21.5,75));
        }
      }
  }
}
