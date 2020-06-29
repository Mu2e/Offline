#ifndef CrvSummaryMC_hh
#define CrvSummaryMC_hh

#include "Rtypes.h"

namespace mu2e
{
  struct CrvSummaryMC
  {
    Double_t            _totalEnergyDeposited;
    Int_t               _nHitCounters; 
    CrvSummaryMC(int totalEnergyDeposited, int nHitCounters) :
              _totalEnergyDeposited(totalEnergyDeposited),
              _nHitCounters(nHitCounters)
              {}
    CrvSummaryMC() :
              _totalEnergyDeposited(0),
              _nHitCounters(0)
              {}
  };

}
#endif


