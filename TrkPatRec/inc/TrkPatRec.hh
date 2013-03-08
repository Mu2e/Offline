// utilities for the Module to perform BaBar Kalman fit
//
// $Id: TrkPatRec.hh,v 1.2 2013/03/08 04:33:26 brownd Exp $
// $Author: brownd $
// $Date: 2013/03/08 04:33:26 $
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

// struct for flagging hits
  struct TrkHitFlag {
#define TIGHTBIT 0x1
#define LOOSEBIT 0x2
#define VERYLOOSEBIT 0x4
#define DELTABIT 0x8
#define STEREOBIT 0x10
    TrkHitFlag() : _iflag(0) {}
    void setVeryLoose() { _iflag |= VERYLOOSEBIT; }
    void setLoose() { _iflag |= LOOSEBIT; }
    void setTight() { _iflag |= TIGHTBIT; }
    void setDelta() { _iflag |= DELTABIT; }
    void setStereo() { _iflag |= STEREOBIT; }
    bool veryLoose() const { return (_iflag & VERYLOOSEBIT) != 0; }
    bool loose() const { return (_iflag & LOOSEBIT) != 0; }
    bool tight() const { return (_iflag & TIGHTBIT) != 0; }
    bool delta() const { return (_iflag & DELTABIT) != 0; }
    bool stereo() const { return(_iflag & STEREOBIT) != 0; }
    unsigned _iflag;
  };
}
