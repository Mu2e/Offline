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
#include <Rtypes.h>
#include <utility>
#include <vector>

namespace mu2e {
// define some structs and types
  // small stub of SimParticle for quick reference to basic information, plus summarize genealogy
  struct SimPartStub {
    typedef art::Ptr<SimParticle> SPPtr;
    PDGCode::type _pdg; // code of this particle
    ProcessCode _proc; // particle creation process
    GenId _gid; // generator code
    MCRelationship _rel; // relationship of this particle to its primary
    uint16_t _nhits; // number of associated StrawHits
    uint16_t _nactive; // number of associated active hits
    SimPartStub() : _pdg(PDGCode::null), _nhits(0), _nactive(0) {}
    // partial constructor from a SimParticle;
    SimPartStub(SPPtr const& spp)  : _pdg(spp->pdgId()),
    _proc(spp->creationCode()), _rel(MCRelationship::none),
    _nhits(0), _nactive(0){
    // dig down to the GenParticle
      if(spp->genParticle().isNonnull()) _gid = spp->genParticle()->generatorId().id();
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
    int16_t simPartStubIndex() const { return _spindex; }
    int16_t _spindex; // index into the associated SimPartStub of this DigiMC
    double _energySum; // sum of all MC true energy deposited
    XYZVec _pos; // in WORLD coordinates
    XYZVec _mom; // momentum of particle at point where digi created
    double _time; // with time maps applied
    double _wireEndTime; // the time that the signal fires TDC
    StrawId _strawId; // the ID of the straw that was hit
    int _gen; // generator ID
    bool _xtalk; // flag if this was a cross-talk hit
  };

  struct KalSeedMC { 
    SimPartStub const& simParticle(size_t index=0) const { return _simps.at(index); }
    std::vector<SimPartStub> const& simParticles() const { return _simps; }
    std::vector<TrkStrawHitMC> const & trkStrawHitMC() const { return _tshmcs; }
    TrkStrawHitMC const& trkStrawHitMC(size_t index) const { return _tshmcs.at(index); }
    SimPartStub const& simParticle(TrkStrawHitMC const& tshmc) const { return simParticle(tshmc.simPartStubIndex()); }
    std::vector<TrkStrawHitMC> const & trkStrawHitMCs() const { return _tshmcs; }
    // data products
    std::vector<SimPartStub> _simps; // associated sim particles, and their relationship
    std::vector<TrkStrawHitMC> _tshmcs;  // MC info for each TrkStrawHitSeed
    std::vector<VDStep> _vdsteps; // sampling of true momentum and position from VDS
  };
  typedef std::vector<KalSeedMC> KalSeedMCCollection;
  typedef art::Assns<KalSeed,KalSeedMC> KalSeedMCAssns;
}
#endif
