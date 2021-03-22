#ifndef MCDataProducts_StepPointMC_hh
#define MCDataProducts_StepPointMC_hh
//
// A persistable class representing a point that is on a track and
// is also inside, or on the boundary of, some G4 volume.  This is used
// for saving information for G4Steps inside sensitive detctors; it is
// also used for the time virtual detector.
//
// More specifically, the point, 3-momentum, time and proper time
// stored in the StepPointMC are the attributes of the G4Track
// at the *beginning* of the just completed G4Step.
//
// The track that made the step is identified by the data member
// _track, which is an art::Ptr to a SimParticle inside a
// SimParticleCollection.
//
// The data members _totalEnergyDeposit, _nonIonizingEnergyDeposit,
// and _stepLength are attributes of the just completed G4Step.
//
// The data member  _endProcessCode is an enum type that describes
// why G4 stopped the current step when it did.  The legal values
// of this data member are found in:
//     MCDataProducts/inc/ProcessCode.hh
//
// The meaning of the data member _volumeId is a little subtle and
// is described below in Note 1).
//
// The Mu2eG4 module can produce many different data products of type
// StepPointMCCollection.  These are distinguished by assigning different
// instance names to the data products.  The Mu2eG4 module can be configured
// to produce all such data products, none or an arbitrary subset.
//
// Each such data product contains StepPointMCs from within a single
// subsystem, such as the straw gas, the calorimeter crystals, the
// CRV planks, the virtual detectors and so on.  The list of all
// legal instance names is found in the enum in:
//     MCDataProducts/inc/StepInstanceName.hh
//
// There is also one special StepPointMCCollection, the timeVD
// or time virtual detector.  When this data product is enabled,
// the user must supply a set of time values.  If a G4Step spans
// one of these times, then a StepPointMC is recorded for that
// G4Step; as always the recorded track attributes are at the
// start of the track.
//
// Original author Rob Kutschke
//
// Notes:
// 1) About the data member _volumeId
//    The meaning of is data member depends on the StepPointMCCollection
//    from which the StepPointMC is taken.  There is no way to
//    know the meaning from the StepPointMC itself.  The table below
//    describes the meaning of _volumId for some of the StepPointMCCollection
//    instances:
//
//    Instance name:     Meaning:
//                       The variable should be used in a call to:
//    tracker          - const Straw& Tracker::getStraw( const StrawId& ) const;
//    calorimeter      - const Crystal& Calorimeter::crystal(int) const
//    calorimeterRO    - Identifies a SiPM attached to a crystal.
//    virtualdetector  - CLHEP::Hep3Vector  const& VirtualDetector::getLocal(int) const;
//                       CLHEP::Hep3Vector  const& VirtualDetector::getGlobal(int) const;
//                       const CLHEP::HepRotation* VirtualDetector::getRotation(int) const;
//    stoppingtarget   - TargetFoil const& StoppingTarget::foil( unsigned int ) const;
//    CRV              - const CRSScintillatorBar& CosmicRayShield::getBar ( CRSScintillatorBarIndex) const;
//
//   In a few cases the argument of the call to the geometry function is a simple integer
//   or unsigned integer.  In other cases, the argument is an opaque integer type
//   whose purpose is to cause compiler errors if you use an object of that type to index into the
//   wrong thing.  Examples of this are StrawId, CRSScintillatorBarIndex and VirtualDetectorId.
//
//   The class StepPointMC also provides convenience functions to cast the _volumeId into one
//   of these opaque types.
//
//   The orignal design was to define an opaque integer type for every geometry subsystem
//   and to use that to design a StepPointMC system with compile time saftey against taking
//   a StepPointMC from the calorimeter and using its volumeId to index into the tracker straws.
//   The effort was stillborn.
//

#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib/map_vector.h"

#include "MCDataProducts/inc/ProcessCode.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <ostream>
#include <vector>

namespace mu2e {

  class StepPointMC{

  public:

    // This might change some day.
    typedef unsigned long VolumeId_type;

    StepPointMC():
      _track(),
      _volumeId(0),
      _totalEnergyDeposit(0.),
      _nonIonizingEnergyDeposit(0.),
      _visibleEnergyDeposit(0.),
      _position(),
      _postPosition(),
      _momentum(),
      _postMomentum(),
      _time(0.),
      _proper(0.),
      _stepLength(0.),
      _endProcessCode(){
    }

    StepPointMC( art::Ptr<SimParticle> const& atrack,
                 VolumeId_type                volumeId,
                 double                       totalEDep,
                 double                       nonIonizingEDep,
                 double                       visEDep,
                 double                       time,
                 double                       proper,
                 CLHEP::Hep3Vector const&     position,
                 CLHEP::Hep3Vector const&     postPosition,
                 CLHEP::Hep3Vector const&     momentum,
                 CLHEP::Hep3Vector const&     postMomentum,
                 double                       stepLength,
                 ProcessCode                  endProcessCode
                 ):
      _track(atrack),
      _volumeId(volumeId),
      _totalEnergyDeposit(totalEDep),
      _nonIonizingEnergyDeposit(nonIonizingEDep),
      _visibleEnergyDeposit(visEDep),
      _position(position),
      _postPosition(postPosition),
      _momentum(momentum),
      _postMomentum(postMomentum),
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
    double                       visibleEDep()      const { return _visibleEnergyDeposit; }
    CLHEP::Hep3Vector const&     position()         const { return _position;  }
    CLHEP::Hep3Vector const&     postPosition()     const { return _postPosition;  }
    CLHEP::Hep3Vector const&     momentum()         const { return _momentum;  }
    CLHEP::Hep3Vector const&     postMomentum()     const { return _postMomentum;  }
    double                       time()             const { return _time;      }
    double                       properTime()       const { return _proper;      }
    double                       stepLength()       const { return _stepLength;}
    ProcessCode                  endProcessCode()   const { return _endProcessCode;}

    // Kept for backwards compatibility.
    double eDep() const { return _totalEnergyDeposit;    }

    // Return the volumeId as a StrawId.
    // This only makes sense for StepPointMCs from the tracker collection.
    // It's the user's job to know if this is a reasonable thing to do.
    StrawId    strawId()    const { return static_cast<StrawId>(_volumeId); }

    // Return the volumeId as a VirtualDetectorId.
    // This only makes sense for StepPointMCs from the virtual detector collection.
    // It's the user's job to know if this is a reasonable thing to do.
    VirtualDetectorId virtualDetectorId() const { return VirtualDetectorId(_volumeId); }

    // Return the volumeId as a CRSScintillatorBarIndex
    // This only makes sense for StepPointMCs from the virtual detector collection.
    // It's the user's job to know if this is a reasonable thing to do.
    CRSScintillatorBarIndex barIndex() const { return CRSScintillatorBarIndex(_volumeId); }

  private:

    // note that the energy deposits and step length are internally floats

    art::Ptr<SimParticle> _track;
    VolumeId_type         _volumeId;
    float                 _totalEnergyDeposit;
    float                 _nonIonizingEnergyDeposit;
    float                 _visibleEnergyDeposit; // used in scintillators
    CLHEP::Hep3Vector     _position;
    CLHEP::Hep3Vector     _postPosition;
    CLHEP::Hep3Vector     _momentum;
    CLHEP::Hep3Vector     _postMomentum;
    double                _time;
    double                _proper;
    float                 _stepLength;
    ProcessCode           _endProcessCode;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   StepPointMC const& h){
    h.print(ost, false);
    return ost;
  }
   typedef std::vector<mu2e::StepPointMC> StepPointMCCollection;

} // namespace mu2e

#endif /* MCDataProducts_StepPointMC_hh */
