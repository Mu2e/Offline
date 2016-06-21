//
// MC relationships of some objects.  Extracted from KalDiag
// Dave Brown, LBNL, 20 Jun 2016
//
#ifndef MCRelationship_HH
#define MCRelationship_HH
// art
#include "art/Persistency/Common/Ptr.h"
// MC data
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
namespace mu2e 
{  
 // some convenient typedefs    
  typedef art::Ptr<SimParticle> SPPtr;
  namespace MCRelationship {
    // values of relationship of 2 MC objects
    enum relation {none=-1,same,daughter,mother,sibling,udaughter,umother,usibling};

// relationship information
    relation relationship(SPPtr const& sppi,SPPtr const& sppj);
//    static relation relationship(SimParticle const& spi,SimParticle const& spj);
    relation relationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2);
  }
}
#endif
