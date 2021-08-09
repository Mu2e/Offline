#include "Offline/DataProducts/inc/XYZVec.hh"
namespace Geom {
  XYZVec toXYZVec(CLHEP::Hep3Vector const& cvec) { return XYZVec(cvec); }
  CLHEP::Hep3Vector Hep3Vec(XYZVec const& rvec) { return
    CLHEP::Hep3Vector((double)rvec.x(),(double)rvec.y(),(double)rvec.z()); }
  CLHEP::HepLorentzVector HepLorentzVec(XYZTVec const& rvec) { return
    CLHEP::HepLorentzVector((double)rvec.x(),(double)rvec.y(),(double)rvec.z(),(double)rvec.t()); }
  // z direction definition; this is missing from GenVector
  XYZVec const& ZDir() {
    static XYZVec _zdir(0.0,0.0,1.0);
    return _zdir;
  }
  std::string XYZnames(const char* vname) {
    std::string svname(vname);
    static std::string leaves; leaves = svname + std::string("x/F:") +
      svname + std::string("y/F:") + svname + std::string("z/F");
    return leaves;
  }
}
