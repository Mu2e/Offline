// Geometrical info about conventional construction pieces for the extinction monitor.
//
// Andrei Gaponenko, 2012

#ifndef EXTMONFNALBUILDING_HH
#define EXTMONFNALBUILDING_HH

#include <vector>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALFilter.hh"

namespace mu2e {

  class ExtMonFNALBuildingMaker;

  class ExtMonFNALBuilding : virtual public Detector {
  public:

    const ExtMonFNALFilter& filter() const { return filter_; }

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
    CLHEP::Hep3Vector shieldingNSteelCenterInMu2e() const { return shieldingNSteelCenterInMu2e_; }
    const std::vector<double>& shieldingNSteelHalfSize() const { return shieldingNSteelHalfSize_; }
    CLHEP::Hep3Vector shieldingSSteelCenterInMu2e() const { return shieldingSSteelCenterInMu2e_; }
    const std::vector<double>& shieldingSSteelHalfSize() const { return shieldingSSteelHalfSize_; }
    CLHEP::Hep3Vector shieldingNConcreteCenterInMu2e() const { return shieldingNConcreteCenterInMu2e_; }
    const std::vector<double>& shieldingNConcreteHalfSize() const { return shieldingNConcreteHalfSize_; }
    CLHEP::Hep3Vector shieldingSConcreteCenterInMu2e() const { return shieldingSConcreteCenterInMu2e_; }
    const std::vector<double>& shieldingSConcreteHalfSize() const { return shieldingSConcreteHalfSize_; }
    CLHEP::Hep3Vector shieldingFloorConcreteCenterInMu2e() const { return shieldingFloorConcreteCenterInMu2e_; }
    const std::vector<double>& shieldingFloorConcreteHalfSize() const { return shieldingFloorConcreteHalfSize_; }

    const CLHEP::Hep3Vector& coll2ShieldingCenterInMu2e() const { return coll2ShieldingCenterInMu2e_; }
    const CLHEP::HepRotation& coll2ShieldingRotationInMu2e() const { return coll2ShieldingRotationInMu2e_; }
    const std::vector<CLHEP::Hep2Vector>& coll2ShieldingOutline() const { return coll2ShieldingOutline_; }

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

    ExtMonFNALFilter filter_;

    double roomInsideFullHeight_;
    double magnetRoomLength_;

    double roomInsideYmin_;
    double roomInsideYmax_;

    double _HVACductRadius;
    double _HVACductHalfLength;
    CLHEP::Hep3Vector HVACductCenterInMu2e_;

    CLHEP::HepRotation shieldingRotationInMu2e_;
    CLHEP::Hep3Vector shieldingNSteelCenterInMu2e_;
    std::vector<double> shieldingNSteelHalfSize_;
    CLHEP::Hep3Vector shieldingSSteelCenterInMu2e_;
    std::vector<double> shieldingSSteelHalfSize_;
    CLHEP::Hep3Vector shieldingNConcreteCenterInMu2e_;
    std::vector<double> shieldingNConcreteHalfSize_;
    CLHEP::Hep3Vector shieldingSConcreteCenterInMu2e_;
    std::vector<double> shieldingSConcreteHalfSize_;
    CLHEP::Hep3Vector shieldingFloorConcreteCenterInMu2e_;
    std::vector<double> shieldingFloorConcreteHalfSize_;

    CLHEP::Hep3Vector coll2ShieldingCenterInMu2e_;
    CLHEP::HepRotation coll2ShieldingRotationInMu2e_;
    std::vector<CLHEP::Hep2Vector> coll2ShieldingOutline_;

    std::vector<double> detectorRoomHalfSize_;
    CLHEP::Hep3Vector detectorRoomCenterInMu2e_;
    CLHEP::HepRotation detectorRoomRotationInMu2e_;
  };


} // namespace mu2e

#endif/*EXTMONFNALBUILDING_HH*/
