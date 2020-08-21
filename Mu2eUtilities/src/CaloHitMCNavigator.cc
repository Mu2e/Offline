//
//
//
// Original author Rob Kutschke
//

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"


namespace mu2e {

  CaloHitMCNavigator::CaloHitMCNavigator( CaloHitCollection const&               hits,
                                          CaloHitMCTruthCollection const&        truth,
                                          CaloHitSimPartMCCollection const&      sims ):
    _hits(&hits),
    _truth(&truth),
    _sims(&sims)
   {}

}
