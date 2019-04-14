//
// Struct describing MC truth for calorimeter cluster 
// original author: Dave Brown (LBNL), Jan 2019
//
#ifndef CaloClusterInfoMC_HH
#define CaloClusterInfoMC_HH
#include "Rtypes.h"
#include "DataProducts/inc/XYZVec.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
namespace mu2e 
{
  struct CaloClusterInfoMC {
    Int_t _nsim; // # of sim particles associated with this cluster
    Float_t _etot; // total true energy from all particles in this cluster
    Float_t _tavg; // average time over all particles
    Float_t _eprimary; // primary particle true energy in this cluster
    Float_t _tprimary; // primary particle average time
    MCRelationship _prel; // relationship of the cluster primary particle to the event primary
    CaloClusterInfoMC() : _nsim(0), _etot(0.0), _tavg(0.0), _eprimary(0.0), _tprimary(0.0){}
    void reset() { *this = CaloClusterInfoMC(); }
    static std::string const& leafnames() { 
      static const std::string leaves = 
	std::string("nsim/I:etot/F:tavg/F:eprimary/F:tprimary/F:prel/B:prem/B");
      return leaves;
    }
  };
}
#endif
