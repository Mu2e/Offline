#ifndef ToyDP_SimParticle_hh
#define ToyDP_SimParticle_hh

//
// Information about particles created by Geant4.
//
// $Id: SimParticle.hh,v 1.7 2010/12/29 18:15:41 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/12/29 18:15:41 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) Internally G4 numbers tracks 1...N.  An earlier version of TrackingAction
//    renumbered them 0...(N-1); this was an artifact of the SimParticleCollection 
//    class being a std::vector, which starts at 0. But now SimParticleCollection 
//    is a MapVector, so it is no longer necessary to do the renumbering.  
//    Therefore the correct test to see if a particle has a parent is
//
// 2) The trackId, parentIds and daughterIds are all of the correct type to be
//    used to find the corresponding information in SimParticleCollection
// 
// 3) I would like to make the PDG id code of type PDGCode::type.  
//    However I am worried that this might screw up the persistency
//    mechanism if it resolves to different types on different machines.
//

#include <vector>

// Mu2e includes
#include "ToyDP/inc/StoppingCode.hh"

// Includes from external packages.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

namespace mu2e {

  struct SimParticle {

    typedef MapVectorKey key_type;

    // This c'tor is required for ROOT.
    SimParticle(){};

    SimParticle( key_type                       aid,
                 key_type                       aparentId,
                 int32_t                        apdgId,
                 int32_t                        agenIndex,
                 const CLHEP::Hep3Vector&       aposition,
                 const CLHEP::HepLorentzVector& amomentum,
                 double                         astartGlobalTime,
                 double                         astartProperTime,
                 uint32_t                       astartVolumeIndex,
                 uint32_t                       astartG4Status,
                 double                         aweight=1.):
      _id(aid),
      _parentId(aparentId),
      _pdgId(apdgId),
      _genIndex(agenIndex),
      _startPosition(aposition),
      _startMomentum(amomentum),
      _startGlobalTime(astartGlobalTime),
      _startProperTime(astartProperTime),
      _startVolumeIndex(astartVolumeIndex),
      _startG4Status(astartG4Status),
      _weight(aweight),
      _endDefined(false)
    {}

    // Accept compiler generated d'tor, copy c'tor and assignment operator.

    void addEndInfo( CLHEP::Hep3Vector       aendPosition,
                     CLHEP::HepLorentzVector aendMomentum,
                     double                  aendGlobalTime,
                     double                  aendProperTime,
                     uint32_t                aendVolumeIndex,
                     uint32_t                aendG4Status,
                     StoppingCode            astoppingCode){
      _endDefined      = true;
      _endPosition     = aendPosition;
      _endMomentum     = aendMomentum;
      _endGlobalTime   = aendGlobalTime;
      _endProperTime   = aendProperTime;
      _endVolumeIndex  = aendVolumeIndex;
      _endG4Status     = aendG4Status;
      _stoppingCode    = astoppingCode;
    }

    void addDaughter( key_type id ){
      _daughterIds.push_back(id);
    }

    // Accessors

    // Index of this track.  See notes 1 and 2.
    key_type id() const {return _id;}

    // Index of the parent of this track.
    key_type  parentId()  const { return _parentId;}
    bool      hasParent() const { return (_parentId != key_type(0)); }

    // PDG particle ID code.  See note 3.
    int32_t pdgId() const {return _pdgId;}

    // Index into the container of generated tracks; 
    // -1 if there is no corresponding generated track.
    int32_t generatorIndex() const { return _genIndex;}
    bool fromGenerator() const { return (_genIndex != -1); }
    bool madeInG4()      const { return (_genIndex == -1); }

    // Information at the start of the track.
    CLHEP::Hep3Vector const& startPosition()       const { return _startPosition;}
    CLHEP::HepLorentzVector const& startMomentum() const { return _startMomentum;}
    double   startGlobalTime()  const { return _startGlobalTime;}
    double   startProperTime()  const { return _startProperTime;}
    uint32_t startVolumeIndex() const { return _startVolumeIndex;}
    uint32_t startG4Status()    const { return _startG4Status;}

    // Information at the end of the track.
    CLHEP::Hep3Vector const& endPosition() const { return _endPosition;}
    CLHEP::HepLorentzVector const& endMomentum() const { return _endMomentum;}
    double       endGlobalTime()  const { return _endGlobalTime; }
    double       endProperTime()  const { return _endProperTime; }
    uint32_t     endVolumeIndex() const { return _endVolumeIndex;}
    uint32_t     endG4Status()    const { return _endG4Status;   }
    StoppingCode stoppingCode()   const { return _stoppingCode;  }

    // SimParticle indices of daughters of this track.
    std::vector<key_type> const& daughterIds() const { return _daughterIds;}

    // Weight
    double weight() const { return  _weight;}

    // Is the second half defined?
    bool endDefined() const { return _endDefined;}

  private:
    // G4 ID number of this track and of its parent.  
    // See notes 1 and 2.
    key_type _id;
    key_type _parentId;

    // PDG particle ID code.  See note 1.
    int32_t _pdgId;

    // Index into the container of generated tracks; 
    // -1 if there is no corresponding generated track.
    int32_t _genIndex;

    // Information at the start of the track.
    CLHEP::Hep3Vector       _startPosition;
    CLHEP::HepLorentzVector _startMomentum;
    double                  _startGlobalTime;
    double                  _startProperTime;
    uint32_t                _startVolumeIndex;
    uint32_t                _startG4Status;

    // Information at the end fo the track.
    CLHEP::Hep3Vector       _endPosition;
    CLHEP::HepLorentzVector _endMomentum;
    double                  _endGlobalTime;
    double                  _endProperTime;
    uint32_t                _endVolumeIndex;
    uint32_t                _endG4Status;

    // The reason that the particle stopped.
    StoppingCode            _stoppingCode;

    // SimParticle IDs of daughters of this track.
    std::vector<key_type>  _daughterIds;

    // Weight
    double _weight;

    // Is the second half defined?
    bool _endDefined;
 
  };

}

#endif
