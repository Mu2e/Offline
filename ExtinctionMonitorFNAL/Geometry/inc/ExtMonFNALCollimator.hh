// A class to keep information to construct an ExtMon collimator.
//
// Andrei Gaponenko, 2012
// Sam Fine and Andrei Gaponenko, 2024

#ifndef EXTMONFNALCOLLIMATOR_HH
#define EXTMONFNALCOLLIMATOR_HH

#include <string>
#include <vector>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

namespace mu2e {

  class ProtonBeamDump;
  class ExtMonFNALFilterMaker;

  //----------------------------------------------------------------
  class ExtMonFNALCollimator {
    friend class ExtMonFNALBuildingMaker; // FIXME: tmp
    friend class ExtMonFNALFilterMaker;
    // Private ctr: the class should be only obtained via the maker
    ExtMonFNALCollimator();
    // but the Filter can start with an uninitialized instance
    friend class ExtMonFNALFilter;


    std::string _name;
    std::vector<double> _channelRadius;
    std::vector<double> _alignmentPlugRadius;
    std::vector<double> _alignmentPlugInnerShellThickness;
    std::vector<double> _alignmentPlugOuterShellThickness;
    std::vector<double> _shotLinerInnerRadius;
    std::vector<double> _shotLinerInnerThickness;
    double _shotLinerOuterRadius;
    double _shotLinerOuterThickness;
    double _length;
    double _radiusTransitiondZ;

    CLHEP::Hep3Vector _centerInMu2e;
    CLHEP::HepRotation _rotationInMu2e;
    double _angleH_inBeamDump;
    double _angleV;

    void setFromDumpAngles(double angleH_inBeamDump, double angleV, const ProtonBeamDump& dump);

  public:
    std::string name() const { return _name; }

    //----------------------------------------------------------------
    const CLHEP::Hep3Vector& centerInMu2e() const { return _centerInMu2e; }
    const CLHEP::HepRotation& rotationInMu2e() const { return _rotationInMu2e; }

    // Point on the straight line along the collimator axis.
    // The position parameter is zero at the collimator center,
    // and increases towards the collimator entrance
    CLHEP::Hep3Vector trajectoryPointInMu2e(double pos) const;

    CLHEP::Hep3Vector entranceInMu2e() const { return trajectoryPointInMu2e(+0.5*_length); }
    CLHEP::Hep3Vector exitInMu2e() const { return trajectoryPointInMu2e(-0.5*_length); }

    double angleH_inBeamDump() const { return _angleH_inBeamDump; }
    double angleV() const { return _angleV; }

    CLHEP::Hep2Vector dxdzdydz() const;

    //----------------------------------------------------------------
    const std::vector<double> &channelRadius() const { return _channelRadius; }

    //alignment plug radius contains inner and outer shell thicknesses
    //alignmentPlugInnerShellThickness corresponds to the collimator channel steel shell
    //alignmentPlugOuterShellThickness corresopnds to the steel shell surrounding the concrete
    const std::vector<double> &alignmentPlugRadius() const { return _alignmentPlugRadius; }
    const std::vector<double> &alignmentPlugInnerShellThickness() const { return _alignmentPlugInnerShellThickness; }
    const std::vector<double> &alignmentPlugOuterShellThickness() const { return _alignmentPlugOuterShellThickness; }

    const std::vector<double> &shotLinerInnerRadius() const { return _shotLinerInnerRadius; }
    const std::vector<double> &shotLinerInnerThickness() const { return _shotLinerInnerThickness; }
    double shotLinerOuterRadius() const { return _shotLinerOuterRadius; }
    double shotLinerOuterThickness() const { return _shotLinerOuterThickness; }
    double length() const { return _length; }

    // 0 means sharp jumb between the radii, >0 is linear change from r1 at -dz to r2 at +dz
    double radiusTransitiondZ() const { return _radiusTransitiondZ; }

  };


} // namespace mu2e

#endif/*EXTMONFNALCOLLIMATOR_HH*/
