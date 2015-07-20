//
//
// $Id: CaloHitMCNavigator.cc,v 1.1 2013/03/08 01:22:32 echenard Exp $
// $Author: echenard $
// $Date: 2013/03/08 01:22:32 $
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
