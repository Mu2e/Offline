#ifndef ToyDP_SimParticle_hh
#define ToyDP_SimParticle_hh

//
// A temporary class to hold particles created by Geant4.
//
// $Id: SimParticle.hh,v 1.3 2010/03/31 19:54:06 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/03/31 19:54:06 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) I would like to make the PDG id code of type PDGCode::type.  
//    However I am worried that this might screw up the persistency
//    mechanism if it resolves to different types on different machines.
//

#include <vector>

// Mu2e includes
#include "Mu2eUtilities/inc/PDGCode.hh"

// Includes from external packages.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

namespace mu2e {

  struct SimParticle {

    // This c'tor is required for ROOT.
    SimParticle(){};

    SimParticle( uint32_t                       aid,
                 int32_t                        aparentId,
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
                     uint32_t                aendG4Status){
      _endDefined      = true;
      _endPosition     = aendPosition;
      _endMomentum     = aendMomentum;
      _endGlobalTime   = aendGlobalTime;
      _endProperTime   = aendProperTime;
      _endVolumeIndex  = aendVolumeIndex;
      _endG4Status     = aendG4Status;
    }

    void addDaughter( uint32_t id ){
      _daughterIds.push_back(id);
    }

    // Accessors

    // Index of this track.
    uint32_t id() const {return _id;}

    // Index of the parent of this track.
    int32_t  parentId()  const { return _parentId;}
    bool     hasParent() const { return (_parentId != -1); }

    // PDG particle ID code.  See note 1.
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

    // Information at the start of the track.
    CLHEP::Hep3Vector const& endPosition() const { return _endPosition;}
    CLHEP::HepLorentzVector const& endMomentum() const { return _endMomentum;}
    double   endGlobalTime()  const { return _endGlobalTime;}
    double   endProperTime()  const { return _endProperTime;}
    uint32_t endVolumeIndex() const { return _endVolumeIndex;}
    uint32_t endG4Status()    const { return _endG4Status;}

    // SimParticle indices of daughters of this track.
    std::vector<uint32_t> const& daughterIds() const { return _daughterIds;}

    // Weight
    double weight() const { return  _weight;}

    // Is the second half defined?
    bool endDefined() const { return _endDefined;}

  private:
    // Id (serial number) of this track and of its parent.  
    // Tracks do not reach TrackingAction in order.
    uint32_t _id;
    int32_t  _parentId;

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

    CLHEP::Hep3Vector       _endPosition;
    CLHEP::HepLorentzVector _endMomentum;
    double                  _endGlobalTime;
    double                  _endProperTime;
    uint32_t                _endVolumeIndex;
    uint32_t                _endG4Status;

    // SimParticle indices of daughters of this track.
    std::vector<uint32_t>   _daughterIds;

    // Weight
    double _weight;

    // Is the second half defined?
    bool _endDefined;
 
  };

}

#endif
