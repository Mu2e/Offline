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
    Float_t             _x, _y, _z;       //position of the first MC particle in CRV
    Int_t               _crvSectorNumber; //CRV sector number of the first MC particle in CR
    Int_t               _crvSectorType;   //CRV sector type of the first MC particle in CRV
    Int_t               _pdgId;           //PDG ID of the first MC particle in CRV

    CrvSummaryMC(double totalEnergyDeposited, double minPathLayer,  double maxPathLayer, int nHitCounters, CLHEP::Hep3Vector pos, int crvSectorNumber, int crvSectorType, int pdgId) :
              _totalEnergyDeposited(totalEnergyDeposited),
              _minPathLayer(minPathLayer),
              _maxPathLayer(maxPathLayer),
              _nHitCounters(nHitCounters),
              _x(pos.x()), _y(pos.y()), _z(pos.z()),
              _crvSectorNumber(crvSectorNumber),
              _crvSectorType(crvSectorType),
              _pdgId(pdgId)
              {}
    CrvSummaryMC() :
              _totalEnergyDeposited(0),
              _minPathLayer(0),
              _maxPathLayer(0),
              _nHitCounters(0),
              _x(-99999), _y(-99999), _z(-99999),
              _crvSectorNumber(-1),
              _crvSectorType(-1),
              _pdgId(-99999)
              {}
  };

}
#endif
