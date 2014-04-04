///////////////////////////////////////////////////////////////////////////////
// utilities for the Module to perform BaBar Kalman fit
// 2014-04-04 P.Murat: cloned from TrkPatRec/inc/TrkPatRec.hh
//
// $Id: CalTimePeak.hh,v 1.1 2014/04/04 21:23:34 murat Exp $
// $Author: murat $
// $Date: 2014/04/04 21:23:34 $
///////////////////////////////////////////////////////////////////////////////
#ifndef __CalPatRec_CalTimePeak_hh__
#define __CalPatRec_CalTimePeak_hh__

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "KalmanTests/inc/TrkDef.hh"
// C++
#include <vector>

using namespace std; 

namespace mu2e {

// struct to keep track of hits in a time peak

  struct CalTimePeak {
    const CaloCluster*    _cluster;	// not owned, just a pointer
    double                _z;
    double                _tpeak;
    std::vector<hitIndex> _trkptrs;

    CalTimePeak(const CaloCluster* Cl,double Z) {
      _cluster = Cl;
      _z       = Z;
    }

    bool operator < (CalTimePeak const& other ) const { return _trkptrs.size() < other._trkptrs.size(); }
    bool operator > (CalTimePeak const& other ) const { return _trkptrs.size() > other._trkptrs.size(); }
  };
}

#endif
