// Geometrical info about conventional construction pieces for the extinction monitor.
//
// Andrei Gaponenko, 2012

#ifndef EXTMONFNALBUILDING_HH
#define EXTMONFNALBUILDING_HH

#include <vector>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "art/Persistency/Common/Wrapper.h"

#include "Mu2eInterfaces/inc/Detector.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"

namespace mu2e {

  class ExtMonFNALBuildingMaker;

  class ExtMonFNALBuilding : virtual public Detector {
  public:

    //----------------------------------------------------------------
    class CollimatorExtMonFNAL {
      friend class ExtMonFNALBuildingMaker;

      std::string _name;
      double _horizontalLength;
      std::vector<double> _channelWidth;
      std::vector<double> _channelHeight;
      std::vector<double> _alignmentPlugRadius;
      std::vector<double> _alignmentHoleRClearance;
      double _radiusTransitiondZ;
      double _angleH;
      double _angleV;

    public:
      std::string name() const { return _name; }
      // NB: "horizontal" is historic name, this is actually the length of projection on dump Z axis.
      double horizontalLength() const { return _horizontalLength; }

      // Half length along axis
      double halfLength() const;

      const std::vector<double> &channelWidth() const { return _channelWidth; }
      const std::vector<double> &channelHeight() const { return _channelHeight; }
      const std::vector<double> &alignmentPlugRadius() const { return _alignmentPlugRadius; }
      const std::vector<double> &alignmentHoleRClearance() const { return _alignmentHoleRClearance; }
      double radiusTransitiondZ() const { return _radiusTransitiondZ; }

      // these two are not just "positioning" parameters but also affect the shape
      // thus they belong to this class
      double angleH() const { return _angleH; }
      double angleV() const { return _angleV; }
    };

    const CLHEP::Hep3Vector& roomRefPointInMu2e() const { return roomRefPointInMu2e_; }
    const CLHEP::HepRotation& roomRotationInMu2e() const { return roomRotationInMu2e_; }

    // The offsets are w.r.t. the dump core center, in the plane of the
    // dump shielding face.
    double filterEntranceOffsetX() const { return _filterEntranceOffsetX; }
    double filterEntranceOffsetY() const { return _filterEntranceOffsetY; }

    double filterAngleH() const { return _filterAngleH; }
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

    // (x,z) points in Mu2e coordinates (not a direct copy of inputs!)
    const std::vector<CLHEP::Hep2Vector>& roomInsideOutline() const { return roomInsideOutline_; }
    const std::vector<CLHEP::Hep2Vector>& wallOutsideOutline() const { return wallOutsideOutline_; }
    // the ceiling is the same as the outside wall
    const std::vector<CLHEP::Hep2Vector>& ceilingOutline() const { return wallOutsideOutline_; }
    // the floor is different because of interference with core back shielding
    const std::vector<CLHEP::Hep2Vector>& floorOutline() const { return floorOutsideOutline_; }

    double roomInsideFullHeight() const { return roomInsideFullHeight_; }

    double roomWallThickness() const { return roomWallThickness_; }
    double roomFloorThickness() const { return roomFloorThickness_; }
    double roomCeilingThickness() const { return roomCeilingThickness_; }

    double dirtOverheadThickness() const { return dirtOverheadThickness_; }
    // min. extension outide of the ceiling outline
    double dirtOverheadHorizontalMargin() const { return dirtOverheadHorizontalMargin_; }

    double magnetRoomLength() const { return magnetRoomLength_; }
    double coll2ShieldingDumpXmin() const { return coll2ShieldingDumpXmin_; }
    double coll2ShieldingDumpXmax() const { return coll2ShieldingDumpXmax_; }

    // computed:
    double roomInsideYmin() const { return roomInsideYmin_; }
    double roomInsideYmax() const { return roomInsideYmax_; }

    CLHEP::Hep3Vector ceilingRefPointInMu2e() const;
    CLHEP::Hep3Vector floorRefPointInMu2e() const;

    CLHEP::Hep3Vector coll2ShieldingCenterInMu2e() const { return coll2ShieldingCenterInMu2e_; }
    CLHEP::HepRotation coll2ShieldingRotationInMu2e() const { return coll2ShieldingRotationInMu2e_; }
    const std::vector<double>& coll2ShieldingHalfSize() const { return coll2ShieldingHalfSize_; }

  private:
    friend class ExtMonFNALBuildingMaker;
      // Private ctr: the class should be only obtained via the maker
    ExtMonFNALBuilding();
    // Or read back from persistent storage
    template<class T> friend class art::Wrapper;

    CLHEP::Hep3Vector roomRefPointInMu2e_;
    CLHEP::HepRotation roomRotationInMu2e_;

    double _filterEntranceOffsetX;
    double _filterEntranceOffsetY;
    double _filterAngleH;
    double _filterEntranceAngleV;

    ExtMonFNALMagnet _filterMagnet;
    CollimatorExtMonFNAL _collimator1;
    CollimatorExtMonFNAL _collimator2;

    CLHEP::Hep3Vector _collimator1CenterInMu2e;
    CLHEP::HepRotation _collimator1RotationInMu2e;

    CLHEP::Hep3Vector _collimator2CenterInMu2e;
    CLHEP::HepRotation _collimator2RotationInMu2e;

    std::vector<CLHEP::Hep2Vector> roomInsideOutline_;
    std::vector<CLHEP::Hep2Vector> wallOutsideOutline_;
    std::vector<CLHEP::Hep2Vector> floorOutsideOutline_;
    double roomInsideFullHeight_;
    double roomWallThickness_;
    double roomFloorThickness_;
    double roomCeilingThickness_;
    double dirtOverheadThickness_;
    double dirtOverheadHorizontalMargin_;

    double magnetRoomLength_;
    double coll2ShieldingDumpXmin_;
    double coll2ShieldingDumpXmax_;

    double roomInsideYmin_;
    double roomInsideYmax_;

    CLHEP::Hep3Vector coll2ShieldingCenterInMu2e_;
    CLHEP::HepRotation coll2ShieldingRotationInMu2e_;
    std::vector<double> coll2ShieldingHalfSize_;
  };


} // namespace mu2e

#endif/*EXTMONFNALBUILDING_HH*/
