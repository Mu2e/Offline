//
// $Id: TrackClusterMatch.cc,v 1.3 2014/06/11 16:21:03 murat Exp $
// $Author: murat $
// $Date: 2014/06/11 16:21:03 $
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

  TrackClusterMatch::TrackClusterMatch(TrkToCaloExtrapolPtr& Tex, CaloClusterPtr& Cluster, Data_t* Data) 
//     _textrapol(Tex),
//     _cluster  (Cluster)
  {
    _textrapol = Tex;
    _cluster   = Cluster;
    _dx        = Data->dx;
    _dy        = Data->dy;
    _dz        = Data->dz;
    _du        = Data->du;
    _dv        = Data->dv;
    _dt        = Data->dt;
    _ep        = Data->ep;
    _chi2      = Data->chi2;
    _chi2_time = Data->chi2_time;
  }

  TrackClusterMatch::~TrackClusterMatch() {
  }


//-----------------------------------------------------------------------------
  void TrackClusterMatch::print(const char* Option) const {
    printf(" >>> WARNING: TrackClusterMatch::print not implemented yet\n");
  }

} // namespace mu2e
