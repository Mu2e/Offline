///////////////////////////////////////////////////////////////////////////////
// utilities for the Module to perform BaBar Kalman fit
// 2014-04-04 P.Murat: cloned from TrkPatRec/inc/TrkPatRec.hh
//
// $Id: CalTimePeak.hh,v 1.3 2014/06/06 21:35:08 murat Exp $
// $Author: murat $
// $Date: 2014/06/06 21:35:08 $
///////////////////////////////////////////////////////////////////////////////
#ifndef __CalPatRec_CalTimePeak_hh__
#define __CalPatRec_CalTimePeak_hh__

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"

// C++
#include <vector>

using namespace std; 

namespace mu2e {

// struct to keep track of hits in a time peak

//  typedef std::vector<hitIndex> hitIndexCollection;
  
  struct CalTimePeak {
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    const CaloCluster*         _cluster;     // cached, not owned
    int                        _cprIndex;    // CalPatRec track index or -1, if no track
    double                     _x;           // cluster X coordinate in the detector system
    double                     _y;		// cluster Y coordinate in the detector system
    double                     _z;		// cluster Z coordinate in the detector system
    double                     _tpeak;       // cluster time ?
    std::vector<StrawHitIndex> _index;       // selects subset of _shcol hits
    double                     _tmin;        // lower bound of the timing window
    double                     _tmax;        // upper bound of the timing window

    const StrawHitCollection*      _shcol;  // cached pointer to the StrawHit collection
    const StrawHitFlagCollection*  _shfcol; // cached pointer to the StrawHitFlag collection
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    CalTimePeak();
    CalTimePeak(const CaloCluster* Cl, double X, double Y, double Z);
    ~CalTimePeak();

    const CaloCluster* Cluster  ()      const { return _cluster; }
    int                CprIndex ()      const { return _cprIndex;}
    double             ClusterX ()      const { return _x;       }
    double             ClusterY ()      const { return _y;       }
    double             ClusterZ ()      const { return _z;       }
    double             ClusterT0()      const { return _cluster->time(); }
    double             TMin     ()      const { return _tmin; }
    double             TMax     ()      const { return _tmax; }
    int                NHits    ()      const { return _index.size(); }
    int                HitIndex (int I) const { return _index.at(I); }

    void               SetCprIndex(int Index) { _cprIndex = Index; }

    void               clear();

    bool operator < (CalTimePeak const& other ) const { return _index.size() < other._index.size(); }
    bool operator > (CalTimePeak const& other ) const { return _index.size() > other._index.size(); }
  };

  typedef std::vector<CalTimePeak> CalTimePeakCollection;
}

#endif
