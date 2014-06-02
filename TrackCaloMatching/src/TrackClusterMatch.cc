//
// $Id: TrackClusterMatch.cc,v 1.1 2014/06/02 22:26:15 murat Exp $
// $Author: murat $
// $Date: 2014/06/02 22:26:15 $
//
// Original author Ivan Logashenko
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "TrackCaloMatching/inc/TrackClusterMatch.hh"

using namespace std;

namespace mu2e {

  TrackClusterMatch::TrackClusterMatch() {
  }

  TrackClusterMatch::TrackClusterMatch(TrkToCaloExtrapolPtr& Tex, CaloClusterPtr& Cluster, double Chi2) 
//     _textrapol(Tex),
//     _cluster  (Cluster)
  {
    _textrapol = Tex;
    _cluster   = Cluster;
    _chi2      = Chi2;
  }

  TrackClusterMatch::~TrackClusterMatch() {
  }


//-----------------------------------------------------------------------------
  void TrackClusterMatch::print(const char* Option) const {
    printf(" >>> WARNING: TrackClusterMatch::print not implemented yet\n");
  }

} // namespace mu2e
