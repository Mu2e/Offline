#ifndef ProductionSolenoidGeom_PSShield_hh
#define ProductionSolenoidGeom_PSShield_hh

//
// $Id: PSShield.hh,v 1.3 2012/04/05 18:43:39 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/05 18:43:39 $
//
// Original author Andrei Gaponenko
//

#include <ostream>

#include "Mu2eInterfaces/inc/Detector.hh"
#include "GeomPrimitives/inc/Polycone.hh"

#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {

  class PSShieldMaker;

  class PSShield : virtual public Detector {

  public:

    const Polycone& bulk() const { return bulk_; }
    const CLHEP::Hep3Vector& originInMu2e() const { return bulk_.originInMu2e(); }

    const CLHEP::Hep3Vector& cutoutRefPoint() const { return cutoutRefPoint_; }

    double cutoutR() const { return cutoutR_; }
    double cutoutHalfLength() const { return cutoutHalfLength_; }
    double cutoutRotY() const { return cutoutRotY_; }

  private:
    friend class PSShieldMaker;

    // Private ctr: the class should only be constructed via PSShield::PSShieldMaker.
    PSShield(const Polycone& bulk,
             const CLHEP::Hep3Vector& cutoutRefPoint,
             double cutoutR,
             double cutoutHalfLength,
             double cutoutRotY)
      : bulk_(bulk)
      , cutoutRefPoint_(cutoutRefPoint)
      , cutoutR_(cutoutR)
      , cutoutHalfLength_(cutoutHalfLength)
      , cutoutRotY_(cutoutRotY)
    {};

    // Or read back from persistent storage
    PSShield();
    template<class T> friend class art::Wrapper;

    Polycone bulk_;
    CLHEP::Hep3Vector cutoutRefPoint_;
    double   cutoutR_;
    double   cutoutHalfLength_;
    double   cutoutRotY_;
  };


  std::ostream& operator<<(std::ostream& os, const PSShield& shield);
}

#endif/*ProductionSolenoidGeom_PSShield_hh*/
