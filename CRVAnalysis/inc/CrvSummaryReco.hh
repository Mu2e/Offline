#ifndef CrvSummaryReco_hh
#define CrvSummaryReco_hh

#include "Rtypes.h"

namespace mu2e
{
  struct CrvSummaryReco
  {
    Int_t               _totalPEs;
    Int_t               _nHitCounters; 
    CrvSummaryReco(int totalPEs, int nHitCounters) :
              _totalPEs(totalPEs),
              _nHitCounters(nHitCounters)
              {}
    CrvSummaryReco() :
              _totalPEs(0),
              _nHitCounters(0)
              {}
  };

}
#endif


