#ifndef RecoDataProducts_DFlag_hh
#define RecoDataProducts_DFlag_hh
//
//  Simple struct to flag a datum within a collection for later selection.
// 
// $Id: DFlag.hh,v 1.1 2013/01/26 18:18:44 brownd Exp $
// $Author: brownd $
// $Date: 2013/01/26 18:18:44 $
//
// Original author David Brown
//
// Mu2e includes
// general includes
#include <algorithm>
#include <math.h>

namespace mu2e {
  struct DFlag {
  // set default values on construction
    DFlag() : _qflag(-1),_pflag(0){}
    // allow merging: properties are 
    void merge(DFlag const& other) {
      _qflag = std::max(_qflag,other._qflag);
      _pflag |= other._pflag;
    }
// quality
    int _qflag;
// property (unspecified)
    unsigned _pflag;
    bool hasProperties(unsigned pmask) const { return (_pflag & pmask) == pmask; }
    bool hasQuality(int qval) const { return _qflag >= qval; }
    bool select(unsigned pmask,int qval) const { return hasProperties(pmask) && hasQuality(qval); }
  };
};
#endif

