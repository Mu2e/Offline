// $Id: TrackClusterLink.hh,v 1.1 2012/07/10 00:02:20 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:20 $
//
// Original author Gianantonio Pezzullo


#ifndef TrackClusterLink_hh
#define TrackClusterLink_hh

#include "art/Persistency/Common/Assns.h"
#include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"

namespace mu2e {
  class TrkToCaloExtrapol;
  class CaloCluster;

  typedef art::Assns<TrkToCaloExtrapol, CaloCluster>  TrackClusterLink;
}

#endif/*TrackClusterLink_hh*/
