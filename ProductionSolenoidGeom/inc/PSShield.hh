#ifndef ProductionSolenoidGeom_PSShield_hh
#define ProductionSolenoidGeom_PSShield_hh

//
//
// Original author Andrei Gaponenko
//
// Added proton beam inlet to HRS.  David N. Brown (Louisville),
// March 2015
// Add upstream and downstream end rings to HRS, January 2017
#include <ostream>
#include <vector>

#include <CLHEP/Geometry/Transform3D.h>

#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/GeomPrimitives/inc/Polycone.hh"

#include "canvas/Persistency/Common/Wrapper.h"

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

    const int& version() const { return version_; } // Versioning implemented by DNB, Jan 2017

    const std::vector<Polycone>& shells() const { return shells_; }
    const std::vector<Polycone>& endRings() const { return endRings_;}

    const std::vector<Groove>& grooves() const { return grooves_; }
    unsigned nGrooves() const { return grooves_.size(); }


    const Tube beamInlet() const { return beamInlet_; }
    const double getBeamAngleY() const { return beamAngleY_; }
    const double getBeamAngleX() const { return beamAngleX_; }
    const CLHEP::Hep3Vector getBeamInletCenter() const { return beamInletCenter_; }
    const bool getCreateBeamPipe() const { return createBeamPipe_; };
  private:

    PSShield();

    // Private ctr: the class should only be constructed via PSShield::PSShieldMaker.
    friend class PSShieldMaker;

    // Or read back from persistent storage
    template<class T> friend class art::Wrapper;

    std::vector<Polycone> shells_;
    std::vector<Groove> grooves_;
    int version_;
    std::vector<Polycone> endRings_;
    Tube beamInlet_;
    double beamAngleY_;
    double beamAngleX_;
    CLHEP::Hep3Vector beamInletCenter_;
    bool createBeamPipe_;
  };

  std::ostream& operator<<(std::ostream& os, const PSShield::Groove& groove);
  std::ostream& operator<<(std::ostream& os, const PSShield& shield);
}

#endif/*ProductionSolenoidGeom_PSShield_hh*/
