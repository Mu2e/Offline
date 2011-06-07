#ifndef MCDataProducts_StepPointMC_hh
#define MCDataProducts_StepPointMC_hh
//
// A persistable class representing a point that is on a track and
// is also inside, or on the boundary of, some G4 volume.  This can be
// used for saving points on the trajectory of the tracking and
// cosmic ray veto systems and for non-senstive material that we wish
// to record for purposes of debugging fitters.  We may need a different
// class to hold the corresponding information for calorimeters.
//
// $Id: StepPointMC.hh,v 1.3 2011/06/07 21:32:21 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/07 21:32:21 $
//
// Original author Rob Kutschke
//

#include "art/Persistency/Common/Ptr.h"
#include "cetlib/map_vector.h"

#include "MCDataProducts/inc/ProcessCode.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "TrackerGeom/inc/StrawIndex.hh"

#include "art/Persistency/Common/OrphanHandle.h"

#include "CLHEP/Vector/ThreeVector.h"

#include <ostream>

namespace mu2e {

  class StepPointMC{

  public:

    // This might change some day.
    typedef unsigned long VolumeId_type;

    StepPointMC():
      _trackId(cet::map_vector_key()),
      _track(),
      _volumeId(0),
      _totalEnergyDeposit(0.),
      _nonIonizingEnergyDeposit(0.),
      _position(),
      _momentum(),
      _time(0.),
      _proper(0.),
      _stepLength(0.),
      _endProcessCode(){
    }

    StepPointMC( unsigned                 trackId,
                 VolumeId_type            volumeId,
                 double                   totalEDep,
                 double                   nonIonizingEDep,
                 double                   time,
                 double                   proper,
                 CLHEP::Hep3Vector const& position,
                 CLHEP::Hep3Vector const& momentum,
                 double                   stepLength,
                 ProcessCode              endProcessCode
                 ):
      _trackId(cet::map_vector_key(trackId)),
      _track(),
      _volumeId(volumeId),
      _totalEnergyDeposit(totalEDep),
      _nonIonizingEnergyDeposit(nonIonizingEDep),
      _position(position),
      _momentum(momentum),
      _time(time),
      _proper(proper),
      _stepLength(stepLength),
      _endProcessCode(endProcessCode){
    }

    // Accept compiler generated versions of:
    //  d'tor
    //  copy c'tor
    //  assignment operator

    void print( std::ostream& ost, bool doEndl = true ) const;
    void print() const { print(std::cout); }

    // Accesors.
    cet::map_vector_key          trackId()          const { return _trackId;   }
    art::Ptr<SimParticle> const& simParticle()      const { return _track;     }
    VolumeId_type                volumeId()         const { return _volumeId;  }
    double                       totalEDep()        const { return _totalEnergyDeposit; }
    double                       nonIonizingEDep()  const { return _nonIonizingEnergyDeposit; }
    double                       ionizingEdep()     const { return _totalEnergyDeposit-_nonIonizingEnergyDeposit; }
    CLHEP::Hep3Vector const&     position()         const { return _position;  }
    CLHEP::Hep3Vector const&     momentum()         const { return _momentum;  }
    double                       time()             const { return _time;      }
    double                       properTime()       const { return _proper;      }
    double                       stepLength()       const { return _stepLength;}
    ProcessCode                  endProcessCode()   const { return _endProcessCode;}

    // Modifier.
    void setSimParticlePtr( art::OrphanHandle<SimParticleCollection>& handle ){
      _track = art::Ptr<SimParticle>(handle,_trackId.asInt());
    }

    // Kept for backwards compatibility.
    double eDep()     const { return _totalEnergyDeposit;    }

    // Return the volumeId as a StrawIndex.
    // It's the user's job to know if this is a reasonable thing to do.
    StrawIndex strawIndex() const { return StrawIndex(_volumeId); }

  private:

    cet::map_vector_key   _trackId;
    art::Ptr<SimParticle> _track;
    VolumeId_type         _volumeId;
    double                _totalEnergyDeposit;
    double                _nonIonizingEnergyDeposit;
    CLHEP::Hep3Vector     _position;
    CLHEP::Hep3Vector     _momentum;
    double                _time;
    double                _proper;
    double                _stepLength;
    ProcessCode           _endProcessCode;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StepPointMC const& h){
    h.print(ost, false);
    return ost;
  }


} // namespace mu2e

#endif /* MCDataProducts_StepPointMC_hh */
