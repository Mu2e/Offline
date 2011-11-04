// G4 specific geometry info.  Not available in non-geant jobs.
// 
// Andrei Gaponennko, 2011

#ifndef WORLDG4_HH
#define WORLDG4_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "GeometryService/inc/Detector.hh"

namespace mu2e {

  class WorldG4Maker;

  class WorldG4 : public Detector {
  public:

    const std::vector<double>& halfLengths() const { return _halfLengths; }

    // All the coordinates are in the G4 world reference frame
    const CLHEP::Hep3Vector& mu2eOriginInWorld() const { return _mu2eOriginInWorld; }
    const CLHEP::Hep3Vector& cosmicReferencePoint() const { return _cosmicReferencePoint; }
    double dirtG4Ymin() const { return _dirtG4Ymin; }
    double dirtG4Ymax() const { return _dirtG4Ymax; }

    // implement Detector's method
    virtual std::string name() const { return "WorldG4"; }

    //----------------------------------------------------------------
  private: 
    friend class WorldG4Maker;

    // Private ctr: the class should be only obtained via the maker
    WorldG4() {}

    std::vector<double> _halfLengths;
    CLHEP::Hep3Vector _mu2eOriginInWorld;
    CLHEP::Hep3Vector _cosmicReferencePoint;
    double _dirtG4Ymin;
    double _dirtG4Ymax;

  };

}

#endif/*WORLDG4_HH*/
