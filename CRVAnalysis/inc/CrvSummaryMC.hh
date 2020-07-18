#ifndef CrvSummaryMC_hh
#define CrvSummaryMC_hh

#include "Rtypes.h"

namespace mu2e
{
  struct CrvSummaryMC
  {
    Double_t            _totalEnergyDeposited;
    Double_t            _minPathLayer;
    Double_t            _maxPathLayer;
    Int_t               _nHitCounters;
    CrvSummaryMC(double totalEnergyDeposited, double minPathLayer,  double maxPathLayer, int nHitCounters) :
              _totalEnergyDeposited(totalEnergyDeposited),
              _minPathLayer(minPathLayer),
              _maxPathLayer(maxPathLayer),
              _nHitCounters(nHitCounters)
              {}
    CrvSummaryMC() :
              _totalEnergyDeposited(0),
	      _minPathLayer(0),
	      _maxPathLayer(0),
              _nHitCounters(0)
              {}
  };

}
#endif
