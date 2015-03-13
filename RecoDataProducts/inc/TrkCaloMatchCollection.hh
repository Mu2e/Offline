//container of the trajectories which have an intersection with the Calorimeter
//
//
// $Id: TrkToCaloExtrapolCollection.hh,v 1.1 2012/07/10 00:02:19 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:19 $
//
// Original author B. Echenard
//


#ifndef RecoDataProducts_TrkCaloMatchCollection_hh
#define RecoDataProducts_TrkCaloMatchCollection_hh


#include <vector>
#include "RecoDataProducts/inc/TrkCaloMatch.hh"


namespace mu2e {
  typedef std::vector<mu2e::TrkCaloMatch> TrkCaloMatchCollection;
}

#endif

