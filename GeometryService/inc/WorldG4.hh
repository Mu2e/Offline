// G4 specific geometry info.  Not available in non-geant jobs.
//
// Geometry parameters that have meaning in the physical world, such
// as object sizes or positions, should be defined in other places.
// This class collects artifacts resulting from the use of a finite
// "world" volume by G4, and from the use by G4 of a coordinate system
// that is different from the Mu2e system.
//
// Andrei Gaponenko, 2011

#ifndef WORLDG4_HH
#define WORLDG4_HH

#include <vector>
#include <ostream>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "art/Persistency/Common/Wrapper.h"

#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class WorldG4Maker;

  class WorldG4 : virtual public Detector {
  public:

    const std::vector<double>& halfLengths() const { return _halfLengths; }

    // The formal hall box
    const std::vector<double>& hallFormalHalfSize() const { return _hallFormalHalfSize; }
    const CLHEP::Hep3Vector&   hallFormalCenterInWorld() const { return _hallFormalCenterInWorld; }

    // All the coordinates are in the G4 world reference frame
    const CLHEP::Hep3Vector& mu2eOriginInWorld() const { return _mu2eOriginInWorld; }

    double dirtG4Ymax() const { return _dirtG4Ymax; }
    double dirtG4Ymin() const { return _dirtG4Ymin; }  // FIXME: not set

    // The argument is a point in Mu2e coordinates.  Return true if this point
    // is inside the G4 world, false otherwise.
    bool inWorld( CLHEP::Hep3Vector const& x0 ) const;

    // The argument is a point in Mu2e coordinates.  Throw an exception if this point
    // is outside the G4 world.
    void inWorldOrThrow( CLHEP::Hep3Vector const& x0 ) const;

    //----------------------------------------------------------------
  private:
    friend class WorldG4Maker;
    template<class T> friend class art::Wrapper; // Needed for persistency

    // Private ctr: the class should be only obtained via the maker
    WorldG4();

    std::vector<double> _halfLengths;
    std::vector<double> _hallFormalHalfSize;
    CLHEP::Hep3Vector _hallFormalCenterInWorld;
    CLHEP::Hep3Vector _mu2eOriginInWorld;
    CLHEP::Hep3Vector _cosmicReferencePoint;
    double _dirtG4Ymin;
    double _dirtG4Ymax;
  };

  std::ostream& operator<<(std::ostream& os, const WorldG4& w);
}

#endif/*WORLDG4_HH*/
