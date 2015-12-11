#ifndef RecoDataProducts_HelixVal_hh
#define RecoDataProducts_HelixVal_hh
//
// helix parameters container
//
// $Id: HelixVal.hh,v 1.1 2012/05/23 07:58:06 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/05/23 07:58:06 $
//
// Original author G. Tassielli
//

// C++ includes
#include <vector>

// Mu2e includes

namespace mu2e {

// copied form TrkDef.hh to avoid compilation problem
struct HitIndex {
  size_t _index; // index into the straw hit container
  int _ambig; // hit ambiguity.  0 means compute from track
  HitIndex() : _index(0),_ambig(0) {}
  HitIndex(size_t index,int ambig=0) : _index(index),_ambig(ambig) {}
  HitIndex& operator = (size_t index) { _index = index; return *this; }
};

struct HelixVal {
        HelixVal() :
                _d0(0.0),
                _phi0(0.0),
                _omega(0.0),
                _z0(0.0),
                _tanDip(0.0)
        {
                _covMtrx[0][0]=1.;        _covMtrx[0][1]=0.;        _covMtrx[0][2]=0.;        _covMtrx[0][3]=0.;        _covMtrx[0][4]=0.;
                _covMtrx[1][0]=0.;        _covMtrx[1][1]=0.1*0.1;   _covMtrx[1][2]=0.;        _covMtrx[1][3]=0.;        _covMtrx[1][4]=0.;
                _covMtrx[2][0]=0.;        _covMtrx[2][1]=0.;        _covMtrx[2][2]=1e-2*1e-2; _covMtrx[2][3]=0.;        _covMtrx[2][4]=0.;
                _covMtrx[3][0]=0.;        _covMtrx[3][1]=0.;        _covMtrx[3][2]=0.;        _covMtrx[3][3]=1.;        _covMtrx[3][4]=0.;
                _covMtrx[4][0]=0.;        _covMtrx[4][1]=0.;        _covMtrx[4][2]=0.;        _covMtrx[4][3]=0.;        _covMtrx[4][4]=0.1*0.1;
        }
        std::vector<HitIndex> _selectedTrackerHitsIdx;
        double _d0;
        double _phi0;
        double _omega;
        double _z0;
        double _tanDip;
        double _covMtrx[5][5];
};

} // namespace mu2e

#endif /* RecoDataProducts_HelixVal_hh */
