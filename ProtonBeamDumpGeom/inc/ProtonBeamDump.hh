// Geometry of the proton beam dump.
//
// Andrei Gaponenko, 2011

#ifndef PROTONBEAMDUMP_HH
#define PROTONBEAMDUMP_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "canvas/Persistency/Common/Wrapper.h"

#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"

namespace mu2e {

  class ProtonBeamDumpMaker;

  class ProtonBeamDump : virtual public Detector {
  public:

    //----------------------------------------------------------------
    // Beam dump core placement and dimensions

    const CLHEP::Hep3Vector& coreCenterInMu2e() const { return _coreCenterInMu2e; }
    // absolute w.r.t to the Mu2e
    double coreRotY() const { return _coreRotY; }
    const CLHEP::HepRotation& coreRotationInMu2e() const { return _coreRotationInMu2e; }

    const std::vector<double>& coreHalfSize() const { return _coreHalfSize; }

    //----------------------------------------------------------------
    // The beam dump concrete poured around the core.  It is modeled
    // as a boolean subtraction volume.  The starting volume is as an
    // extruded solid parameterized here.  Some subtraction pieces are
    // defined below: the "mouth" and the "neutron cave", as well as
    // the ExtMon room cutout.  Another subtraction is based on the
    // Entrace collimator for the Extinction Monitor, its parameters
    // are not provided by this class.
    //
    // The parameters for the mouth and neutron cave subtraction
    // volumes here are nominal.  The margins to ensure proper
    // overlaps for G4 booleans are added by the code talking to G4.
    // OTOH extMonSubtractionOutline is tweaked here, as postponing
    // the extruded solid tweak would make the overall code too entangled.

    std::vector<CLHEP::Hep2Vector> dumpConcreteOutline() const { return _dumpConcreteOutline; }
    const CLHEP::Hep3Vector& dumpConcreteCenterInMu2e() const { return _dumpConcreteCenterInMu2e; }
    double dumpConcreteHalfHeight() const { return _dumpConcreteHalfHeight; }

    std::vector<CLHEP::Hep2Vector> extMonSubtractionOutline() const { return _extMonSubtractionOutline; }
    const CLHEP::Hep3Vector& extMonSubtractionCenterInMu2e() const { return _extMonSubtractionCenterInMu2e; }
    double extMonSubtractionHalfHeight() const { return _extMonSubtractionHalfHeight; }

    const std::vector<double>& mouthHalfSize() const { return _mouthHalfSize; }
    const CLHEP::Hep3Vector& mouthCenterInMu2e() const { return _mouthCenterInMu2e; }

    const std::vector<double>& neutronCaveHalfSize() const { return _neutronCaveHalfSize; }
    const CLHEP::Hep3Vector& neutronCaveCenterInMu2e() const { return _neutronCaveCenterInMu2e; }

    //----------------------------------------------------------------
    // More details of the beam dump construction

    // the gap between the core and the concrete for the cooling air
    const std::vector<double>& coreAirHalfSize() const { return _coreAirHalfSize; }
    const CLHEP::Hep3Vector& coreAirCenterInMu2e() const { return _coreAirCenterInMu2e; }

    // The shielding steel on top of the core
    const std::vector<double>& topSteelFlatHalfSize() const { return _topSteelFlatHalfSize; }
    const CLHEP::Hep3Vector& topSteelFlatCenterInMu2e() const { return _topSteelFlatCenterInMu2e; }

    //----------------------------------------------------------------
    // Transform to the "beam dump" coordinate system, which is centered
    // at the core center, and is rotated around the Y axis w.r.t the mu2e system

    CLHEP::Hep3Vector mu2eToBeamDump_position(const CLHEP::Hep3Vector& mu2epos) const;
    CLHEP::Hep3Vector mu2eToBeamDump_momentum(const CLHEP::Hep3Vector& mu2emom) const;

    CLHEP::Hep3Vector beamDumpToMu2e_position(const CLHEP::Hep3Vector& dumppos) const;
    CLHEP::Hep3Vector beamDumpToMu2e_momentum(const CLHEP::Hep3Vector& dumpmom) const;

    //----------------------------------------------------------------
  private:
    friend class ProtonBeamDumpMaker;
    // Private ctr: the class should be only obtained via ProtonBeamDumpMaker
    ProtonBeamDump();

    CLHEP::Hep3Vector _coreCenterInMu2e;
    double _coreRotY;
    CLHEP::HepRotation _coreRotationInMu2e;

    std::vector<double> _coreHalfSize;

    std::vector<CLHEP::Hep2Vector> _dumpConcreteOutline;
    double _dumpConcreteHalfHeight;
    CLHEP::Hep3Vector _dumpConcreteCenterInMu2e;

    std::vector<CLHEP::Hep2Vector> _extMonSubtractionOutline;
    double _extMonSubtractionHalfHeight;
    CLHEP::Hep3Vector _extMonSubtractionCenterInMu2e;

    std::vector<double> _mouthHalfSize;
    CLHEP::Hep3Vector _mouthCenterInMu2e;

    std::vector<double> _neutronCaveHalfSize;
    CLHEP::Hep3Vector _neutronCaveCenterInMu2e;

    std::vector<double> _coreAirHalfSize;
    CLHEP::Hep3Vector _coreAirCenterInMu2e;

    std::vector<double> _topSteelFlatHalfSize;
    CLHEP::Hep3Vector _topSteelFlatCenterInMu2e;
  };
}

#endif/*PROTONBEAMDUMP_HH*/
