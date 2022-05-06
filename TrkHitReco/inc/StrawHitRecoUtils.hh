#ifndef TrkHitReco_StrawHitRecoUtils_hh
#define TrkHitReco_StrawHitRecoUtils_hh

#include <vector>
#include "TH1F.h"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"

namespace mu2e {
  class StrawHitRecoUtils {
    public:
      StrawHitRecoUtils(double pbtOffset, mu2e::TrkHitReco::FitType fittype, unsigned npre, float invnpre, float invgainAvg, float* invgain, int diagLevel, TH1F* maxiter,
          mu2e::StrawIdMask mask, size_t nplanes, size_t npanels, bool writesh, float minT, float maxT, float minE, float maxE, bool filter, bool flagXT,
          float ctE, float ctMinT, float ctMaxT, bool usecc, float clusterDt, size_t numDigis) :
        _pbtOffset(pbtOffset), _fittype(fittype), _npre(npre), _invnpre(invnpre), _invgainAvg(invgainAvg), _invgain(invgain), _diagLevel(diagLevel),
        _maxiter(maxiter), _mask(mask), _npanels(npanels), _writesh(writesh), _minT(minT), _maxT(maxT), _minE(minE), _maxE(maxE),
        _filter(filter), _flagXT(flagXT), _ctE(ctE), _ctMinT(ctMinT), _ctMaxT(ctMaxT), _usecc(usecc), _clusterDt(clusterDt)
    {
      if (!_filter && _flagXT){
        hits_by_panel = std::vector<std::vector<size_t> >(nplanes*npanels,std::vector<size_t>());
        largeHits.clear();
        largeHitPanels.clear();
        largeHits.reserve(numDigis);
        largeHitPanels.reserve(numDigis);
      }
    };

      void flagCrossTalk(std::unique_ptr<mu2e::StrawHitCollection> const& shCol,
          std::unique_ptr<mu2e::ComboHitCollection> const& chCol);
      bool createComboHit(size_t isd, std::unique_ptr<mu2e::ComboHitCollection> const& chCol,
          std::unique_ptr<mu2e::StrawHitCollection> const& shCol,
          const mu2e::CaloClusterCollection *caloClusters,
          mu2e::StrawId const& sid, mu2e::TrkTypes::TDCValues const& tdc, mu2e::TrkTypes::TOTValues const& tot,
          mu2e::TrkTypes::ADCValue const& pmp, mu2e::TrkTypes::ADCWaveform const& waveform,
          mu2e::TrackerStatus const& trackerStatus,  mu2e::StrawResponse const& srep, mu2e::Tracker const& tt);

      float peakMinusPedAvg(mu2e::TrkTypes::ADCWaveform const& adcData) const;
      float peakMinusPed(mu2e::StrawId id, mu2e::TrkTypes::ADCWaveform const& adcData) const;
      float peakMinusPedFirmware(mu2e::StrawId id, mu2e::TrkTypes::ADCValue const& pmp) const;

    private:
      float _pbtOffset;
      mu2e::TrkHitReco::FitType _fittype;
      unsigned _npre;
      float _invnpre;
      float _invgainAvg;
      float *_invgain;
      int _diagLevel;
      TH1F *_maxiter;
      mu2e::StrawIdMask _mask;
      size_t _npanels;
      bool _writesh;
      float _minT, _maxT, _minE, _maxE;
      bool _filter, _flagXT;
      float _ctE, _ctMinT, _ctMaxT;
      bool _usecc;
      float _clusterDt;

      std::vector<std::vector<size_t> > hits_by_panel;
      std::vector<size_t> largeHits;
      std::vector<size_t> largeHitPanels;
  };

}

#endif
