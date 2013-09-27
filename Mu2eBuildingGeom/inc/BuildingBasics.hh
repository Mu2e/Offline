// Some basic civil construction parameters should be outside of
// Mu2eBuilding class.  The latter describes the complicated shape of
// the hall and depends on e.g. ProtonBeamDump position and rotation.
// Since ProtonBeamDump itself needs to know e.g. hall floor level,
// this building parameter must be in an independent class.
//
// Andrei Gaponenko, 2011

#ifndef BUILDINGBASICS_HH
#define BUILDINGBASICS_HH

#include "Mu2eInterfaces/inc/Detector.hh"
#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {
  class BuildingBasicsMaker;
  class Mu2eBuilding;

  class BuildingBasics : virtual public Detector {
  public:
    // mu2e.origin.heightAboveHallFloor == -detectorHallFloorTopY
    double detectorHallFloorTopY() const { return detectorHallFloorTopY_; }

    double detectorHallInsideFullHeight() const { return detectorHallInsideFullHeight_; }
    double detectorHallCeilingThickness() const { return detectorHallCeilingThickness_; }
    double detectorHallInnerTSCeilingThickness() const { return detectorHallInnerTSCeilingThickness_; }
    double detectorHallFloorThickness()   const { return detectorHallFloorThickness_; }

    double detectorHallFloorTopDepthBelowGrade() const { return detectorHallFloorTopDepthBelowGrade_; }
    // the level of dirt outside the formal box volume
    double yFlatEarth() const { return detectorHallFloorTopDepthBelowGrade_ + detectorHallFloorTopY_; }

  private:
    friend class Mu2eBuilding;
    friend class BuildingBasicsMaker;

    // Private ctr: the class should be only obtained via the maker
    BuildingBasics()
      : detectorHallFloorTopY_(0.)
      , detectorHallInsideFullHeight_(0.)
      , detectorHallCeilingThickness_(0.)
      , detectorHallInnerTSCeilingThickness_(0.)
      , detectorHallFloorThickness_(0.)
      , detectorHallFloorTopDepthBelowGrade_(0.)
    {}

    // Or read back from persistent storage
    template<class T> friend class art::Wrapper;

    double detectorHallFloorTopY_;
    double detectorHallInsideFullHeight_;
    double detectorHallCeilingThickness_;
    double detectorHallInnerTSCeilingThickness_;
    double detectorHallFloorThickness_;
    double detectorHallFloorTopDepthBelowGrade_;
  };
}

#endif/*BUILDINGBASICS_HH*/
