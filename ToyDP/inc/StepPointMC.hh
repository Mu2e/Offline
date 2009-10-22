#ifndef StepPointMC_h
#define StepPointMC_h 1
//
// A persistable class representing a point that is on a track and
// is also inside, or on the boundary of, some G4 volume.  This can be
// used for saving points on the trajectory of the tracking and 
// cosmic ray veto systems and for non-senstive material that we wish 
// to record for purposes of debugging fitters.  We may need a different 
// class to hold the corresponding information for calorimeters.
//
// $Id: StepPointMC.hh,v 1.4 2009/10/22 19:56:10 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/22 19:56:10 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <ostream>

// Mu2e includes
#include "LTrackerGeom/inc/StrawIndex.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e { 

  class StepPointMC{

  public:

    // This might change some day.
    typedef unsigned long VolumeId_type;

    StepPointMC():
      _trackId(-1),
      _volumeId(0),
      _edep(0.),
      _time(0.),
      _position(),
      _momentum(){
    }
    
    StepPointMC( int                      trackId,
		 VolumeId_type            volumeId,
		 double                   edep,
		 double                   time,
		 CLHEP::Hep3Vector const& position,
		 CLHEP::Hep3Vector const& momentum
		 ):
      _trackId(trackId),
      _volumeId(volumeId),
      _edep(edep),
      _time(time),
      _position(position),
      _momentum(momentum){
    }
    
    // Accept compiler generated versions of:
    //  d'tor
    //  copy c'tor 
    //  assignment operator
    
    void print( std::ostream& ost, bool doEndl = true ) const;
    void print() const { print(std::cout); }

    int                      trackId()  const { return _trackId; }
    double                   eDep()     const { return _edep;    } 
    double                   time()     const { return _time;    }
    VolumeId_type            volumeId() const { return _volumeId; }
    CLHEP::Hep3Vector const& position() const { return _position; }
    CLHEP::Hep3Vector const& momentum() const { return _momentum; }

    // If the volume is a straw, then return the volumeId as a StrawIndex.
    // It's the user's job to know if this is a reasonable thing to do.
    StrawIndex strawIndex() const { return StrawIndex::fromInt(_volumeId); }

  private:
  
    int               _trackId;
    VolumeId_type     _volumeId;
    double            _edep;
    CLHEP::Hep3Vector _position;
    CLHEP::Hep3Vector _momentum;
    double            _time;
    
  };

  inline std::ostream& operator<<( std::ostream& ost,
				   StepPointMC const& h){
    h.print(ost, false);
    return ost;
  }
  

} // namespace mu2e

#endif
