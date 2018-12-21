//
// MC relationships of some objects.  Extracted from KalDiag
// Dave Brown, LBNL, 20 Jun 2016
//
#ifndef MCRelationship_HH
#define MCRelationship_HH
// art
#include "canvas/Persistency/Common/Ptr.h"
// MC data
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
namespace mu2e 
{  
 // some convenient typedefs    
  typedef art::Ptr<SimParticle> SPPtr;
  struct MCRelationship {
    // values of relationship of 2 MC objects
    enum relation {none=-1,same,daughter,mother,sibling,udaughter,umother,usibling};
    relation _rel;
    // convenience operators
    bool operator ==(relation rval) const { return _rel == rval; }
    bool operator !=(relation rval) const { return _rel != rval; }
    // trivial constructor
    MCRelationship(relation rval=none) : _rel(rval) {}
    // construct from SimParticles
    MCRelationship(SPPtr const& sppi,SPPtr const& sppj);
    // construct from StrawDigiMCs
    MCRelationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2);
    // from the mixture
    MCRelationship(StrawDigiMC const& mcd, SPPtr const& spp);
  };
}
#endif
