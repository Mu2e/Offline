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
    MCRelationship _rel; // relationship of this particle to its primary
    uint16_t _nhits; // number of associated StrawHits
    uint16_t _nactive; // number of associated active hits
    float _ecalo; // Associated calorimeter cluster energy
    SimPartStub() : _pdg(PDGCode::null), _nhits(0), _nactive(0), _ecalo(0.0) {}
    // partial constructor from a SimParticle;
    SimPartStub(SPPtr const& spp)  : _pdg(spp->pdgId()),
    _proc(spp->creationCode()), _rel(MCRelationship::none),
    _nhits(0), _nactive(0), _ecalo(0.0) {}
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
