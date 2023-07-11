#include "Offline/DataProducts/inc/GenVector.hh"
namespace mu2e {
  namespace GenVector {
    CLHEP::Hep3Vector Hep3Vec(XYZVectorF const& rvec) { return
      CLHEP::Hep3Vector((double)rvec.x(),(double)rvec.y(),(double)rvec.z()); }
    CLHEP::HepLorentzVector HepLorentzVec(XYZTVectorF const& rvec) { return
      CLHEP::HepLorentzVector((double)rvec.x(),(double)rvec.y(),(double)rvec.z(),(double)rvec.t()); }
    CLHEP::Hep3Vector Hep3Vec(XYZVectorD const& rvec) { return
      CLHEP::Hep3Vector((double)rvec.x(),(double)rvec.y(),(double)rvec.z()); }
    CLHEP::HepLorentzVector HepLorentzVec(XYZTVectorD const& rvec) { return
      CLHEP::HepLorentzVector((double)rvec.x(),(double)rvec.y(),(double)rvec.z(),(double)rvec.t()); }
    // z direction definition; this is missing from GenVector
    XYZVectorF const& ZDir() {
      static XYZVectorF _zdir(0.0,0.0,1.0);
      return _zdir;
    }
  }
}
