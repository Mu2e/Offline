#ifndef GeometryService_DUSAFMu2eConverter_hh
#define GeometryService_DUSAFMu2eConverter_hh

//
// Convert coordinates between DUSAF and Mu2e coordinate systems
// See Mu2e-doc-4741-v6
//
// DUSAF coordinates in US feet (see footnote on page 4 for US foot vs International foot).
// Mu2e coordinates in mm
//
// Original author Rob Kutschke
//

#include "Offline/Mu2eInterfaces/inc/Detector.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

namespace mu2e {

  class DUSAFMu2eConverter: virtual public Detector{

  public:

    // Unit conversion factors
    constexpr static double us_foot_mm   = 304.80061;  // Number of mm in 1 US foot (AKA surveyors foot).
    constexpr static double intl_foot_mm = 304.8;      // Number of mm in 1 international foot.
    constexpr static double us_foot_m    = 0.30480061; // Number of m in 1 US foot (AKA surveyors foot).
    constexpr static double intl_foot_m  = 0.3048;     // Number of m in 1 international foot.

    DUSAFMu2eConverter();

    // Convert Mu2e in mm to DUSAF in us feet.
    CLHEP::Hep3Vector Mu2e_to_DUSAF( CLHEP::Hep3Vector const & ) const;

    // Convert DUSAF in us feet to mm.
    CLHEP::Hep3Vector DUSAF_to_Mu2e( CLHEP::Hep3Vector const & ) const;

    CLHEP::HepRotation const& rotationToMu2e()  const {return _toMu2e;  }
    CLHEP::HepRotation const& rotationToDUSAF() const {return _toDUSAF; }

  private:

    // Tie in points inthe two coordinate systems
    const CLHEP::Hep3Vector _tieIn_inDUSAF_usft; // in US feet
    const CLHEP::Hep3Vector _tieIn_inMu2e_mm;    // in mm

    // Rotation component of the transformation.
    const CLHEP::HepRotation _toDUSAF;
    const CLHEP::HepRotation _toMu2e;

  };

}

#endif
