#ifndef MCDataProducts_CaloClusterMC_hh
#define MCDataProducts_CaloClusterMC_hh
//
// Summarize the truth content of a CaloCluster
//
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include <vector>
#include <utility>
namespace mu2e {
  struct CaloMCEDep {
    art::Ptr<SimParticle> const& sim() const { return _simp; }
    float energyDeposit() const { return _edep; }
    float time() const { return _time; }
    art::Ptr<SimParticle> _simp; // sim particle
    MCRelationship _rel; // relationship of this particle to its primary
    float _edep; // energy deposited by this particle in this cluster
    float _time; // average time of the energy deposition by this particle; includes all offsets!
    // should include an energy-weighted position FIXME!
  };
    
  struct CaloClusterMC {
    std::vector<CaloMCEDep> const& energyDeposits() const { return _edeps; }
    float totalEnergyDeposit() const { return _edep; }
    float averageTime() const { return _time; }
    
    std::vector<CaloMCEDep> _edeps; // energy contributions from individual SimParticles, sorted by contribution
    float _edep; // total edep
    float _time; // energy-weighted time average, including all offsets!
    // should include an energy-weighted position FIXME!
  };
  typedef std::vector<mu2e::CaloClusterMC> CaloClusterMCCollection;
  typedef art::Assns<CaloCluster,CaloClusterMC> CaloClusterMCAssns;
}
#endif
