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
// $Id: StepPointMC.hh,v 1.9 2013/09/06 16:14:00 gandr Exp $
// $Author: gandr $
// $Date: 2013/09/06 16:14:00 $
//
// Original author Rob Kutschke
//

#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib/map_vector.h"

#include "MCDataProducts/inc/ProcessCode.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "DataProducts/inc/StrawIndex.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <ostream>

namespace mu2e {

  class StepPointMC{

  public:

    // This might change some day.
    typedef uint16_t VolumeId_type;

    StepPointMC():
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

    StepPointMC( art::Ptr<SimParticle> const& atrack,
                 VolumeId_type                volumeId,
                 double                       totalEDep,
                 double                       nonIonizingEDep,
                 double                       time,
                 double                       proper,
                 CLHEP::Hep3Vector const&     position,
                 CLHEP::Hep3Vector const&     momentum,
                 double                       stepLength,
                 ProcessCode                  endProcessCode
                 ):
      _track(atrack),
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

    // Accept compiler generated versions of: d'tor, copy c'tor, assignment operator

    // Printing PDG ID of the particle requires that the _track Ptr is resolvable.
    // This is not the case until the SimParticle collection is committed, i.e.
    // printing PDG ID should not be enabled from withing G4 module.
    void print( std::ostream& ost, bool doEndl = true, bool printPDGId = true ) const;
    void print() const { print(std::cout); }

    // Accesor and modifier; the modifier is needed for event mixing.
    art::Ptr<SimParticle> const& simParticle() const { return _track; }
    art::Ptr<SimParticle>&       simParticle()       { return _track; }

    // Accesors.
    cet::map_vector_key trackId() const {
      return ( _track.isNonnull() ) ? cet::map_vector_key(_track.key()): cet::map_vector_key(0);
    }

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

    // Kept for backwards compatibility.
    double eDep() const { return _totalEnergyDeposit;    }

    // Return the volumeId as a StrawIndex.
    // This only makes sense for StepPointMCs from the tracker collection.
    // It's the user's job to know if this is a reasonable thing to do.
    StrawIndex strawIndex() const { return StrawIndex(_volumeId); }

    // Return the volumeId as a VirtualDetectorId.
    // This only makes sense for StepPointMCs from the virtual detector collection.
    // It's the user's job to know if this is a reasonable thing to do.
    VirtualDetectorId virtualDetectorId() const { return VirtualDetectorId(_volumeId); }

    // Return the volumeId as a CRSScintillatorBarIndex
    // This only makes sense for StepPointMCs from the virtual detector collection.
    // It's the user's job to know if this is a reasonable thing to do.
    CRSScintillatorBarIndex barIndex() const { return CRSScintillatorBarIndex(_volumeId); }

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
