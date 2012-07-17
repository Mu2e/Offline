#ifndef TrackCaloMatching_TrackClusterLink_hh
#define TrackCaloMatching_TrackClusterLink_hh

//
// $Id: TrackClusterLink.hh,v 1.1 2012/07/17 20:00:26 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/17 20:00:26 $
//
// Original author Gianantonio Pezzullo

#include "art/Persistency/Common/Assns.h"
#include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

namespace mu2e {

  typedef art::Assns<TrkToCaloExtrapol, CaloCluster>  TrackClusterLink;

}

#endif/*TrackCaloMatching_TrackClusterLink_hh*/
