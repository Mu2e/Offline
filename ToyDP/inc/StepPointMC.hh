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
// $Id: StepPointMC.hh,v 1.10 2010/04/04 20:35:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/04 20:35:13 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <ostream>

// Mu2e includes
#include "TrackerGeom/inc/StrawIndex.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e { 

  class StepPointMC{

  public:

    // This might change some day.
    typedef unsigned long VolumeId_type;

    StepPointMC():
      _trackId(0),
      _volumeId(0),
      _totalEDep(0.),
      _position(),
      _momentum(),
      _time(0.),
      _stepLength(0.){
    }
    
    StepPointMC( uint32_t                 trackId,
                 VolumeId_type            volumeId,
                 double                   totalEDep,
                 double                   time,
                 CLHEP::Hep3Vector const& position,
                 CLHEP::Hep3Vector const& momentum,
                 double                   stepLength
                 ):
      _trackId(trackId),
      _volumeId(volumeId),
      _totalEDep(totalEDep),
      _position(position),
      _momentum(momentum),
      _time(time),
      _stepLength(stepLength){
    }
    
    // Accept compiler generated versions of:
    //  d'tor
    //  copy c'tor 
    //  assignment operator
    
    void print( std::ostream& ost, bool doEndl = true ) const;
    void print() const { print(std::cout); }

    uint32_t                 trackId()    const { return _trackId;   }
    VolumeId_type            volumeId()   const { return _volumeId;  }
    double                   totalEDep()  const { return _totalEDep; } 
    CLHEP::Hep3Vector const& position()   const { return _position;  }
    CLHEP::Hep3Vector const& momentum()   const { return _momentum;  }
    double                   time()       const { return _time;      }
    double                   stepLength() const { return _stepLength;}

    // Kept for backwards compatibility.
    double                   eDep()     const { return _totalEDep;    } 

    // Return the volumeId as a StrawIndex.
    // It's the user's job to know if this is a reasonable thing to do.
    StrawIndex strawIndex() const { return StrawIndex(_volumeId); }

  private:
  
    uint32_t          _trackId;
    VolumeId_type     _volumeId;
    double            _totalEDep;
    CLHEP::Hep3Vector _position;
    CLHEP::Hep3Vector _momentum;
    double            _time;
    double            _stepLength;
    
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StepPointMC const& h){
    h.print(ost, false);
    return ost;
  }
  

} // namespace mu2e

#endif
