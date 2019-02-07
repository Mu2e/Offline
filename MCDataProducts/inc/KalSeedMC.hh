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
    XYZVec _pos;
    XYZVec _mom;
    double _time;
    VirtualDetectorId _vdid;
    VDStep() : _time(0.0) {}
    VDStep(StepPointMC const& vdstep) : _pos(Geom::toXYZVec(vdstep.position())),
    _mom(Geom::toXYZVec(vdstep.momentum())),
    _time(vdstep.time()),
    _vdid(vdstep.virtualDetectorId()) {}
  };

  struct KalSeedMC { 
    SimPartStub const& simParticle(size_t index=0) const { return _simps.at(index); }
    std::vector<SimPartStub> const& simParticles() const { return _simps; }
    // data products
    std::vector<SimPartStub> _simps; // associated sim particles, and their relationship
    std::vector<int16_t> _digisimps; // reference into sim particles for each TrkStrawHitSeed
    std::vector<VDStep> _vdsteps; // sampling of true momentum and position from VDS
  };
  typedef std::vector<KalSeedMC> KalSeedMCCollection;
  typedef art::Assns<KalSeed,KalSeedMC> KalSeedMCAssns;
}
#endif
