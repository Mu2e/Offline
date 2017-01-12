// Andrei Gaponenko, 2012

#include "ProductionSolenoidGeom/inc/PSShield.hh"

#include <iterator>
#include <algorithm>

namespace mu2e {

  // genreflex persistency requires default ctr
  PSShield::Groove::Groove() : theta_(), phi_(), r_(), halfLength_() {}

  PSShield::Groove::Groove(const CLHEP::Hep3Vector& ref, double theta, double phi, double r, double hl)
    //: placement_(CLHEP::HepRotationZ(phi)*CLHEP::HepRotationY(-theta), ref)
    : placement_(HepGeom::Translate3D(ref)*HepGeom::RotateZ3D(phi)*HepGeom::RotateY3D(-theta))
    , theta_(theta)
    , phi_(phi)
    , r_(r)
    , halfLength_(hl)
  {}

  // genreflex persistency requires default ctr
  PSShield::PSShield() : shells_(), version_(), endRings_() {}

  std::ostream& operator<<(std::ostream& os, const PSShield::Groove& g) {
    return os<<"Groove("
             <<"ref="<<g.refPoint()
             <<", theta="<<g.theta()
             <<", phi="<<g.phi()
             <<", r="<<g.r()
             <<", halfLength="<<g.halfLength()
             <<")";
  }


  std::ostream& operator<<(std::ostream& os, const PSShield& shield) {
    os<<"PSShield( shells={ ";
    std::copy(shield.shells().begin(), shield.shells().end(),
              std::ostream_iterator<Polycone>(os, ", "));
    os<<" } , grooves={ ";
    std::copy(shield.grooves().begin(), shield.grooves().end(),
              std::ostream_iterator<PSShield::Groove>(os, ", "));
    os<<" } )";
    return os;
  }

} // namespace mu2e
