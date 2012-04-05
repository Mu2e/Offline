// G4 specific geometry info.  Not available in non-geant jobs.
//
// Andrei Gaponenko, 2011

#ifndef WORLDG4_HH
#define WORLDG4_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "art/Persistency/Common/Wrapper.h"

#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class WorldG4Maker;

  class WorldG4 : virtual public Detector {
  public:

    const std::vector<double>& halfLengths() const { return _halfLengths; }

    // The formal boundary of the "HallAir" volume, in Mu2e coordinates.
    // Set to accommodate ProtonBeamDump and ExtMonFNAL volumes inside.
    double hallFormalZminInMu2e() const { return _hallFormalZminInMu2e; }

    // All the coordinates are in the G4 world reference frame
    const CLHEP::Hep3Vector& mu2eOriginInWorld() const { return _mu2eOriginInWorld; }
    const CLHEP::Hep3Vector& cosmicReferencePoint() const { return _cosmicReferencePoint; }
    double dirtG4Ymin() const { return _dirtG4Ymin; }
    double dirtG4Ymax() const { return _dirtG4Ymax; }

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
    WorldG4() {}

    std::vector<double> _halfLengths;
    CLHEP::Hep3Vector _mu2eOriginInWorld;
    CLHEP::Hep3Vector _cosmicReferencePoint;
    double _dirtG4Ymin;
    double _dirtG4Ymax;
    double _hallFormalZminInMu2e;

  };

}

#endif/*WORLDG4_HH*/
