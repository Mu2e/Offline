//
// struct to place variables we want to test in TrkQual
// 
#ifndef TrkQualTestInfo_HH
#define TrkQualTestInfo_HH
#include "DataProducts/inc/threevec.hh"
#include "TrkDiag/inc/helixpar.hh"
#include "Rtypes.h"
namespace mu2e
{
// general information about a track
  struct TrkQualTestInfo {
    Float_t _longCon;      // consistency between the longitudinal position from the fit with the previous time division measurement
    TrkQualTestInfo() { reset(); }
    void reset() { 
      _longCon=-1000.0;
    }
    static std::string const& leafnames() { static const std::string leaves =
	std::string("longCon/F");
     return leaves;
    }
  };
}
#endif

