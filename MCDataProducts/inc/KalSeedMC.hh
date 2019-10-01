#ifndef MCDataProducts_KalSeedMC_hh
#define MCDataProducts_KalSeedMC_hh
//
//  MC truth match to a KalSeed (Kalman fit track)
//  Original Author: Dave Brown (LBNL) 5 Feb. 2019
//
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "DataProducts/inc/PDGCode.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "art/Framework/Principal/Handle.h"
#include "cetlib/map_vector.h"
#include <Rtypes.h>
#include <utility>
#include <vector>

namespace mu2e {
// define some structs and types
  // small stub of SimParticle for quick reference to basic information, plus summarize genealogy
  struct SimPartStub {
    typedef art::Ptr<SimParticle> SPPtr;
    typedef art::Handle<SimParticleCollection> SPCH;
    PDGCode::type _pdg; // code of this particle
    ProcessCode _proc; // particle creation process
    GenId _gid; // generator code
    MCRelationship _rel; // relationship of this particle to its primary
    uint16_t _nhits; // number of associated StrawHits
    uint16_t _nactive; // number of associated active hits
    XYZVec _mom; // initial momentum 
    cet::map_vector_key _spkey; // key to the SimParticle
    // construct a Ptr from Handle and key
    SPPtr simParticle(SPCH spcH) const { return SPPtr(spcH,_spkey.asUint()); }
    SimPartStub() : _pdg(PDGCode::null), _nhits(0), _nactive(0) {}
    // partial constructor from a SimParticle;
    SimPartStub(SPPtr const& spp)  : _pdg(spp->pdgId()),
    _proc(spp->creationCode()), _gid(GenId::unknown), _rel(MCRelationship::none),
    _nhits(0), _nactive(0), _mom(Geom::toXYZVec(spp->startMomentum())), _spkey(spp.key()){
    // dig down to the GenParticle
      auto simPtr = spp;
      while (simPtr->genParticle().isNull() && simPtr->parent().isNonnull()) {
	simPtr = simPtr->parent();
      }
      if(simPtr->genParticle().isNonnull())_gid = simPtr->genParticle()->generatorId();
    }
  };
  // sampled pair of momentum and position (tracker system) of the primary matched particle
  // These come from the virtual detectors
  struct VDStep {
    XYZVec _pos;  // postion in DETECTOR COORDINATES
    XYZVec _mom;
    double _time;
    VirtualDetectorId _vdid;
    VDStep() : _time(0.0) {}
    VDStep(CLHEP::Hep3Vector const& pos,CLHEP::Hep3Vector const& mom, double time, VirtualDetectorId const& vdid) :
      
      _pos(Geom::toXYZVec(pos)),
    _mom(Geom::toXYZVec(mom)),
    _time(time),
    _vdid(vdid) {}
  };
//
// MC information for TrackStrawHits on this fit
  struct TrkStrawHitMC {
    StrawHitIndex strawDigiMCIndex() const { return _sdmcindex; }
    StrawHitIndex simPartStubIndex() const { return _spindex; }
    StrawId const& strawid() const { return _strawId; }
    float energySum() const { return _energySum; }
    float stepTime() const { return _time; }
    XYZVec const& clusterPosition() const { return _cpos; }
    XYZVec const& particleMomentum() const { return _mom; }
    StrawHitIndex _sdmcindex; // index into the original StrawDigiMC collection
    StrawHitIndex _spindex; // index into the associated SimPartStub of this DigiMC
    StrawId _strawId; // the ID of the straw that was hit
    float _energySum; // sum of all MC true energy deposited by trigger particles
    float _time; // time of trigger StepPoint with time maps applied, wrapped to the beam
    XYZVec _cpos; // trigger cluster position in detector coordinates
    XYZVec _mom; // momentum of particle at point where digi created
  };

  struct KalSeedMC { 
    SimPartStub const& simParticle(size_t index=0) const { return _simps.at(index); }
    std::vector<SimPartStub> const& simParticles() const { return _simps; }
    std::vector<TrkStrawHitMC> const & trkStrawHitMCs() const { return _tshmcs; }
    TrkStrawHitMC const& trkStrawHitMC(size_t index) const { return _tshmcs.at(index); }
    SimPartStub const& simParticle(TrkStrawHitMC const& tshmc) const { return simParticle(tshmc.simPartStubIndex()); }
    // data products
    std::vector<SimPartStub> _simps; // associated sim particles, and their relationship
    std::vector<TrkStrawHitMC> _tshmcs;  // MC info for each TrkStrawHitSeed
    std::vector<VDStep> _vdsteps; // sampling of true momentum and position from VDS
  };
  typedef std::vector<KalSeedMC> KalSeedMCCollection;
  typedef art::Assns<KalSeed,KalSeedMC> KalSeedMCAssns;
}
#endif
