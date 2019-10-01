#ifndef MCDataProducts_SimParticle_hh
#define MCDataProducts_SimParticle_hh

//
// Information about particles created by Geant4.
//
// $Id: SimParticle.hh,v 1.12 2013/09/27 16:03:41 gandr Exp $
// $Author: gandr $
// $Date: 2013/09/27 16:03:41 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) Internally G4 numbers tracks 1...N.  An earlier version of TrackingAction
//    renumbered them 0...(N-1); this was an artifact of the SimParticleCollection
//    class being a std::vector, which starts at 0. But now SimParticleCollection
//    is a cet::map_vector, so it is no longer necessary to do the renumbering.
//    Therefore the correct test to see if a particle has a parent is
//
// 2) The trackId, parentIds and daughterIds are all of the correct type to be
//    used to find the corresponding information in SimParticleCollection
//
// 3) The trackId, parentIds and daughterIds are all redundant now that the Ptr
//    versions are available.  We will get rid of them as soon as we check
//    backwards compatibility.

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "DataProducts/inc/PDGCode.hh"

#include "canvas/Persistency/Common/Ptr.h"

#include "cetlib/map_vector.h"

#include <iostream>
#include <vector>

namespace mu2e {

  struct SimParticle {

    typedef cet::map_vector_key key_type;

    // A default c'tor is required for ROOT.
    SimParticle():
      _id(),
      _stageOffset(),
      _parentSim(),
      _pdgId(),
      _genParticle(),
      _startPosition(),
      _startMomentum(),
      _startGlobalTime(0.),
      _startProperTime(0.),
      _startVolumeIndex(0),
      _startG4Status(),
      _creationCode(),
      _endPosition(),
      _endMomentum(),
      _endGlobalTime(0.),
      _endProperTime(0.),
      _endVolumeIndex(0),
      _endG4Status(),
      _stoppingCode(),
      _preLastStepKE(-1.),
      _endKE(-1.),
      _nSteps(0),
      _trackLength(-1.),
      _daughterSims(),
      _endDefined(false){
    }

    SimParticle( key_type                       aid,
                 unsigned                       stageOffset,
                 art::Ptr<SimParticle> const&   aparentSim,
                 PDGCode::type                  apdgId,
                 art::Ptr<GenParticle> const&   agenParticle,
                 const CLHEP::Hep3Vector&       aposition,
                 const CLHEP::HepLorentzVector& amomentum,
                 double                         astartGlobalTime,
                 double                         astartProperTime,
                 unsigned                       astartVolumeIndex,
                 unsigned                       astartG4Status,
                 ProcessCode                    acreationCode):
      _id(aid),
      _stageOffset(stageOffset),
      _parentSim(aparentSim),
      _pdgId(apdgId),
      _genParticle(agenParticle),
      _startPosition(aposition),
      _startMomentum(amomentum),
      _startGlobalTime(astartGlobalTime),
      _startProperTime(astartProperTime),
      _startVolumeIndex(astartVolumeIndex),
      _startG4Status(astartG4Status),
      _creationCode(acreationCode),
      _endPosition(),
      _endMomentum(),
      _endGlobalTime(),
      _endProperTime(),
      _endVolumeIndex(),
      _endG4Status(),
      _stoppingCode(),
      _preLastStepKE(-1),
      _endKE(-1),
      _nSteps(0),
      _trackLength(-1.),
      _daughterSims(),
      _endDefined(false)
    {}

    // Accept compiler generated d'tor, copy c'tor and assignment operator.

    void addEndInfo( CLHEP::Hep3Vector       aendPosition,
                     CLHEP::HepLorentzVector aendMomentum,
                     double                  aendGlobalTime,
                     double                  aendProperTime,
                     unsigned                aendVolumeIndex,
                     unsigned                aendG4Status,
                     ProcessCode             astoppingCode,
                     float                   endKE,
                     int                     nSteps,
		     double                  trackLength){
      _endDefined      = true;
      _endPosition     = aendPosition;
      _endMomentum     = aendMomentum;
      _endGlobalTime   = aendGlobalTime;
      _endProperTime   = aendProperTime;
      _endVolumeIndex  = aendVolumeIndex;
      _endG4Status     = aendG4Status;
      _stoppingCode    = astoppingCode;
      _preLastStepKE   = -1.0;      
      _endKE           = endKE;
      _nSteps          = nSteps;
      _trackLength     = trackLength;
    }

    void addDaughter( art::Ptr<SimParticle> const& p ){
      _daughterSims.push_back(p);
    }

    // Some Accessors/Modifier pairs.
    // The modifiers are needed by the event mixing code.

    // Key of this track within the SimParticleCollection.  See notes 1 and 2.
    key_type  id() const {return _id;}
    key_type& id()       { return _id;}

    unsigned stageOffset() const { return _stageOffset; }

    // The parent of this track; may be null.
    art::Ptr<SimParticle> const& parent() const { return _parentSim; }
    art::Ptr<SimParticle>&       parent()       { return _parentSim; }

    // work back through the genealogy to find the earliest representation of this
    // physical particle. This bridges the gaps created in staged MC production, where particles
    // are stopped at a detector volume.
    SimParticle const& originParticle() const {
      return selfParent() ? parent()->originParticle() : *this;
    }

    // The genparticle corresponding to this track; may be null.
    art::Ptr<GenParticle> const& genParticle() const { return _genParticle;}
    art::Ptr<GenParticle>&       genParticle()       { return _genParticle;}

    // Most members have only accessors

    // PDG particle ID code.
    PDGCode::type pdgId() const {return _pdgId;}

    // Where was this particle created: in the event generator or in G4?
    bool isSecondary()   const { return _parentSim.isNonnull(); }
    bool isPrimary()     const { return _genParticle.isNonnull(); }
    bool selfParent()    const { return _parentSim.isNonnull() && _creationCode == ProcessCode::mu2ePrimary; }

    // Some synonyms for the previous two accessors.
    bool hasParent()     const { return _parentSim.isNonnull(); }
    bool fromGenerator() const { return _genParticle.isNonnull(); }
    bool madeInG4()      const { return _genParticle.isNull();    }

    // Information at the start of the track.
    CLHEP::Hep3Vector const& startPosition()       const { return _startPosition;}
    CLHEP::HepLorentzVector const& startMomentum() const { return _startMomentum;}
    double      startGlobalTime()  const { return _startGlobalTime;}
    double      startProperTime()  const { return _startProperTime;}
    unsigned    startVolumeIndex() const { return _startVolumeIndex;}
    unsigned    startG4Status()    const { return _startG4Status;}
    ProcessCode creationCode()     const { return _creationCode;   }

    // Information at the end of the track.
    CLHEP::Hep3Vector const& endPosition() const { return _endPosition;}
    CLHEP::HepLorentzVector const& endMomentum() const { return _endMomentum;}
    double       endGlobalTime()  const { return _endGlobalTime; }
    double       endProperTime()  const { return _endProperTime; }
    unsigned     endVolumeIndex() const { return _endVolumeIndex;}
    unsigned     endG4Status()    const { return _endG4Status;   }
    ProcessCode  stoppingCode()   const { return _stoppingCode;  }
    double       preLastStepKineticEnergy() const { return _preLastStepKE; }
    float        endKineticEnergy() const { return _endKE; }
    int          nSteps()         const { return _nSteps;        }
    double       trackLength()    const { return _trackLength;   }

    // SimParticle daughters of this track.
    std::vector<art::Ptr<SimParticle> > const& daughters()   const { return _daughterSims; }
    std::vector<art::Ptr<SimParticle> >&       daughters()         { return _daughterSims; }

    // SimParticle indices of daughters of this track.
    // DO NOT USE - this is an expensive (at run time) crutch for legacy code.
    std::vector<key_type>                      daughterIds() const;

    // Is the second half defined?
    bool endDefined() const { return _endDefined;}

    // Modifiers;
    void setDaughterPtrs  ( std::vector<art::Ptr<SimParticle> > const& ptr){
      _daughterSims.clear();
      _daughterSims.reserve(ptr.size());
      _daughterSims.insert( _daughterSims.begin(), ptr.begin(), ptr.end() );
    }

    // Two older accessors that will soon be removed from the interface.

    // Index of the parent of this track; may be null.
    key_type                     parentId() const {
      return ( _parentSim.isNonnull()) ? key_type(_parentSim.key()) : key_type(0);
    }

    // Index into the container of generated tracks;
    // -1 if there is no corresponding generated track.
    int                          generatorIndex() const {
      return ( _genParticle.isNonnull()) ? _genParticle.key() : -1;
    }

  private:
    // G4 ID number of this track and of its parent.
    // See notes 1 and 2.
    key_type _id;

    // The offset between G4 track number and _id above, which
    // uniquely identifies the Mu2eG4 simulation stage that produced
    // this particle
    unsigned _stageOffset;

    art::Ptr<SimParticle> _parentSim;

    // PDG particle ID code.  See note 1.
    PDGCode::type _pdgId;

    art::Ptr<GenParticle>  _genParticle;

    // Information at the start of the track.
    CLHEP::Hep3Vector       _startPosition;
    CLHEP::HepLorentzVector _startMomentum;
    double                  _startGlobalTime;
    double                  _startProperTime;
    unsigned                _startVolumeIndex;
    unsigned                _startG4Status;
    ProcessCode             _creationCode;

    // Information at the end of the track.
    CLHEP::Hep3Vector       _endPosition;
    CLHEP::HepLorentzVector _endMomentum;
    double                  _endGlobalTime;
    double                  _endProperTime;
    unsigned                _endVolumeIndex;
    unsigned                _endG4Status;
    ProcessCode             _stoppingCode;
    double                  _preLastStepKE;
    float                   _endKE;
    int                     _nSteps;
    double                  _trackLength;

    // SimParticle IDs of daughters of this track.
    std::vector<art::Ptr<SimParticle> > _daughterSims;

    // Is the second half defined?
    bool _endDefined;

  };

}

#endif /* MCDataProducts_SimParticle_hh */
