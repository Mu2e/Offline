//
// $Id: TrackClusterMatch.cc,v 1.4 2014/06/24 19:08:39 murat Exp $
// $Author: murat $
// $Date: 2014/06/24 19:08:39 $
//
// Original author Ivan Logashenko
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib_except/exception.h"

// Mu2e includes
#include "RecoDataProducts/inc/TrackClusterMatch.hh"

using namespace std;

namespace mu2e {

  TrackClusterMatch::TrackClusterMatch() {
  }

  TrackClusterMatch::TrackClusterMatch(TrkCaloIntersectPtr& Tex, CaloClusterPtr& Cluster, Data_t* Data) 
  {
    _icl       = Data->icl;
    _iex       = Data->iex;

    _textrapol = Tex;
    _cluster   = Cluster;

    _xtrk      = Data->xtrk;
    _ytrk      = Data->ytrk;
    _ztrk      = Data->ztrk;
    _ttrk      = Data->ttrk;

    _nx        = Data->nx;
    _ny        = Data->ny;
    _nz        = Data->nz;

    _dx        = Data->dx;
    _dy        = Data->dy;
    _dz        = Data->dz;

    _du        = Data->du;
    _dv        = Data->dv;

    _dt        = Data->dt;
    _ep        = Data->ep;

    _chi2      = Data->chi2;
    _chi2_time = Data->chi2_time;

    _int_depth = Data->int_depth;
    _ds        = Data->ds;
    _dr        = Data->dr;
    _sint      = Data->sint;
  }

  TrackClusterMatch::~TrackClusterMatch() {
  }


//-----------------------------------------------------------------------------
  void TrackClusterMatch::print(const char* Option) const {
    printf(" >>> WARNING: TrackClusterMatch::print not implemented yet\n");
  }

} // namespace mu2e
