#ifndef ProductionSolenoidGeom_PSShield_hh
#define ProductionSolenoidGeom_PSShield_hh

//
// $Id: PSShield.hh,v 1.4 2012/04/27 05:37:32 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/27 05:37:32 $
//
// Original author Andrei Gaponenko
//

#include <ostream>
#include <vector>

#include <CLHEP/Geometry/Transform3D.h>

#include "Mu2eInterfaces/inc/Detector.hh"
#include "GeomPrimitives/inc/Polycone.hh"

#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {

  class PSShieldMaker;

  class PSShield : virtual public Detector {
  public:

    class Groove {
      HepGeom::Transform3D placement_;
      double theta_;
      double phi_;
      double r_;
      double halfLength_;

    public:
      CLHEP::Hep3Vector refPoint() const { return placement_.getTranslation(); }
      double theta() const { return theta_; }
      double phi() const { return phi_; }

      // combines all of the above info on position & rotation
      const HepGeom::Transform3D& placement() const { return placement_; }

      double r() const { return r_; }
      double halfLength() const { return halfLength_; }

      //
      Groove(const CLHEP::Hep3Vector& ref, double theta, double phi, double r, double hl);
      Groove(); // persistency for vector<Groove> requires default ctr.
    };

    const Polycone& bulk() const { return bulk_; }
    const CLHEP::Hep3Vector& originInMu2e() const { return bulk_.originInMu2e(); }

    const std::vector<Groove>& grooves() const { return grooves_; }
    unsigned nGrooves() const { return grooves_.size(); }

  private:
    friend class PSShieldMaker;

    // Private ctr: the class should only be constructed via PSShield::PSShieldMaker.
    PSShield(const Polycone& bulk) : bulk_(bulk) {};

    // Or read back from persistent storage
    PSShield();
    template<class T> friend class art::Wrapper;

    Polycone bulk_;
    std::vector<Groove> grooves_;
  };


  std::ostream& operator<<(std::ostream& os, const PSShield::Groove& groove);
  std::ostream& operator<<(std::ostream& os, const PSShield& shield);
}

#endif/*ProductionSolenoidGeom_PSShield_hh*/
