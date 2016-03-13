// utilities for the Module to perform BaBar Kalman fit
//
// $Id: TrkPatRec.hh,v 1.4 2014/05/31 14:28:10 brownd Exp $
// $Author: brownd $
// $Date: 2014/05/31 14:28:10 $
//

#include "TrkReco/inc/TrkDef.hh"
// C++
#include <vector>
#include "Rtypes.h"

using namespace std; 

namespace mu2e 
{
// struct to keep track of hits in a time peak
  struct TrkTimePeak {
    std::vector<hitIndex> _trkptrs;
    double _tpeak; // average peak time
    double _peakmax; // value of the peak maximum
    double _phi; // average resolved phi
    TrkTimePeak(double tpeak,double ymax) : _tpeak(tpeak),_peakmax(ymax), _phi(0.0) {}
    bool operator < (TrkTimePeak const& other ) const { return _trkptrs.size() < other._trkptrs.size(); }
    bool operator > (TrkTimePeak const& other ) const { return _trkptrs.size() > other._trkptrs.size(); }
  };

// struct for time peak hit diagnostics
  struct TimePeakHitInfo {
    Float_t _dt; // time relative to the peak
    Float_t _dphi; // resolved phi relative to the peak
    Float_t _rho; // transverse radius of the peak
    Float_t _mva; // peak MVA output
    Int_t _mcpdg, _mcgen, _mcproc; // MC truth info
  };

  struct TimePeakMVA {
    std::vector<Double_t> _pars;
    Double_t& _dt;
    Double_t& _dphi;
    Double_t& _rho;
    TimePeakMVA() : _pars(3,0.0), _dt(_pars[0]), _dphi(_pars[1]), _rho(_pars[2]) {}
  };
}

