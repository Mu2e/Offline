// utilities for the Module to perform BaBar Kalman fit
//
// $Id: TrkPatRec.hh,v 1.3 2013/04/04 01:09:21 brownd Exp $
// $Author: brownd $
// $Date: 2013/04/04 01:09:21 $
//

#include "KalmanTests/inc/TrkDef.hh"
// C++
#include <vector>

using namespace std; 

namespace mu2e 
{
// struct to keep track of hits in a time peak
  struct TrkTimePeak {
    std::vector<hitIndex> _trkptrs;
    double _tpeak;
    double _peakmax;
    TrkTimePeak(double tpeak,double ymax) : _tpeak(tpeak),_peakmax(ymax) {}
    bool operator < (TrkTimePeak const& other ) const { return _trkptrs.size() < other._trkptrs.size(); }
    bool operator > (TrkTimePeak const& other ) const { return _trkptrs.size() > other._trkptrs.size(); }
  };
}
