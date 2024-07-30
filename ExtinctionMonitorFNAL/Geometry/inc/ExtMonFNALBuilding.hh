// Geometrical info about conventional construction pieces for the extinction monitor.
//
// Andrei Gaponenko, 2012

#ifndef EXTMONFNALBUILDING_HH
#define EXTMONFNALBUILDING_HH

#include <vector>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "canvas/Persistency/Common/Wrapper.h"

#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"

namespace mu2e {

  class ExtMonFNALBuildingMaker;

  class ExtMonFNALBuilding : virtual public Detector {
  public:

    //----------------------------------------------------------------
    class CollimatorExtMonFNAL {
      friend class ExtMonFNALBuildingMaker;
      friend class ExtMonFNALBuilding;
      // Private ctr: the class should be only obtained via the maker
      CollimatorExtMonFNAL();

      std::string _name;
      double _horizontalLength;
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
      double _angleH;
      double _angleV;

    public:
      std::string name() const { return _name; }
      // NB: "horizontal" is historic name, this is actually the length of projection on dump Z axis.
      double horizontalLength() const { return _horizontalLength; }

      // Half length along axis
      double halfLength() const;

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

      // these two are not just "positioning" parameters but also affect the shape
      // thus they belong to this class
      double angleH() const { return _angleH; }
      double angleV() const { return _angleV; }
    };

    // The axes of the entrance and exit collimator lie in a vertical plane.
    // This is the angle between that plane and the "z-axis" of the beam dump.
    double filterAngleH() const { return _filterAngleH; }

    // The angle between the collimator axis and the horizontal plane
    double filterEntranceAngleV() const { return _filterEntranceAngleV; }

    const ExtMonFNALMagnet& filterMagnet() const { return _filterMagnet; }
    const CollimatorExtMonFNAL& collimator1() const { return _collimator1; }
    const CollimatorExtMonFNAL& collimator2() const { return _collimator2; }

    const CLHEP::Hep3Vector& collimator1CenterInMu2e() const { return _collimator1CenterInMu2e; }
    const CLHEP::HepRotation& collimator1RotationInMu2e() const { return _collimator1RotationInMu2e; }

    const CLHEP::Hep3Vector& collimator2CenterInMu2e() const { return _collimator2CenterInMu2e; }
    const CLHEP::HepRotation& collimator2RotationInMu2e() const { return _collimator2RotationInMu2e; }

    // convenience accessor for event generators
    // returns the point of intersection of collimator1 axis with the dump shielding face
    CLHEP::Hep3Vector filterEntranceInMu2e() const;
    // collimator2 exit
    CLHEP::Hep3Vector filterExitInMu2e() const;

    //----------------------------------------------------------------
    // Civil construction details

    double roomInsideFullHeight() const { return roomInsideFullHeight_; }
    double magnetRoomLength() const { return magnetRoomLength_; }

    double roomInsideYmin() const { return roomInsideYmin_; }
    double roomInsideYmax() const { return roomInsideYmax_; }

    double HVACductRadius() const { return _HVACductRadius; }
    double HVACductHalfLength() const { return _HVACductHalfLength; }
    CLHEP::Hep3Vector HVACductCenterInMu2e() const { return HVACductCenterInMu2e_; }

    //----------------------------------------------------------------
    // Shielding that is not part of the civil construction

    CLHEP::HepRotation shieldingRotationInMu2e() const { return shieldingRotationInMu2e_; }
    CLHEP::Hep3Vector shieldingNCenterInMu2e() const { return shieldingNCenterInMu2e_; }
    CLHEP::Hep3Vector shieldingSCenterInMu2e() const { return shieldingSCenterInMu2e_; }
    CLHEP::Hep3Vector shieldingBCenterInMu2e() const { return shieldingBCenterInMu2e_; }
    const std::vector<double>& shieldingNHalfSize() const { return shieldingNHalfSize_; }
    const std::vector<double>& shieldingSHalfSize() const { return shieldingSHalfSize_; }
    const std::vector<double>& shieldingBHalfSize() const { return shieldingBHalfSize_; }

    CLHEP::Hep3Vector coll2ShieldingCenterInMu2e() const { return coll2ShieldingCenterInMu2e_; }
    CLHEP::HepRotation coll2ShieldingRotationInMu2e() const { return coll2ShieldingRotationInMu2e_; }
    std::vector<CLHEP::Hep2Vector> coll2ShieldingOutline() const { return coll2ShieldingOutline_; }

    //----------------------------------------------------------------
    // The "detector room" here does not correspond to a physical object.
    // This is an intermediate box volume that fits between the non-rectangular
    // concrete outer wall and collimator2 shielding and contains the
    // detectors and their supports.

    std::vector<double> detectorRoomHalfSize() const { return detectorRoomHalfSize_; }
    CLHEP::Hep3Vector detectorRoomCenterInMu2e() const { return detectorRoomCenterInMu2e_; }
    CLHEP::HepRotation detectorRoomRotationInMu2e() const { return detectorRoomRotationInMu2e_; }


  private:
    friend class ExtMonFNALBuildingMaker;
    // Private ctr: the class should be only obtained via the maker
    ExtMonFNALBuilding();
    // Or read back from persistent storage
    template<class T> friend class art::Wrapper;

    double _filterAngleH;
    double _filterEntranceAngleV;

    ExtMonFNALMagnet _filterMagnet;
    CollimatorExtMonFNAL _collimator1;
    CollimatorExtMonFNAL _collimator2;

    CLHEP::Hep3Vector _collimator1CenterInMu2e;
    CLHEP::HepRotation _collimator1RotationInMu2e;

    CLHEP::Hep3Vector _collimator2CenterInMu2e;
    CLHEP::HepRotation _collimator2RotationInMu2e;

    double roomInsideFullHeight_;
    double magnetRoomLength_;

    double roomInsideYmin_;
    double roomInsideYmax_;

    double _HVACductRadius;
    double _HVACductHalfLength;
    CLHEP::Hep3Vector HVACductCenterInMu2e_;

    CLHEP::HepRotation shieldingRotationInMu2e_;
    CLHEP::Hep3Vector shieldingNCenterInMu2e_;
    CLHEP::Hep3Vector shieldingSCenterInMu2e_;
    CLHEP::Hep3Vector shieldingBCenterInMu2e_;
    std::vector<double> shieldingNHalfSize_;
    std::vector<double> shieldingSHalfSize_;
    std::vector<double> shieldingBHalfSize_;

    CLHEP::Hep3Vector coll2ShieldingCenterInMu2e_;
    CLHEP::HepRotation coll2ShieldingRotationInMu2e_;
    std::vector<CLHEP::Hep2Vector> coll2ShieldingOutline_;

    std::vector<double> detectorRoomHalfSize_;
    CLHEP::Hep3Vector detectorRoomCenterInMu2e_;
    CLHEP::HepRotation detectorRoomRotationInMu2e_;
  };


} // namespace mu2e

#endif/*EXTMONFNALBUILDING_HH*/
