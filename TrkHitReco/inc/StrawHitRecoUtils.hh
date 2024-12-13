#ifndef TrkHitReco_StrawHitRecoUtils_hh
#define TrkHitReco_StrawHitRecoUtils_hh

#include <vector>
#include <memory>
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"
#include "Offline/TrkHitReco/inc/PeakFit.hh"

namespace mu2e {
  class EventWindowMarker;
  class StrawHitRecoUtils {
    public:
      using ADCWFIter = TrkTypes::ADCWaveform::const_iterator;
      StrawHitRecoUtils(TrkHitReco::FitType fittype, int diagLevel,
          StrawIdMask mask, bool writesh,
          float minT, float maxT,
          bool overrideminTOff, float minTOff, bool overridemaxTOff, float maxTOff,
          float minE, float maxE, float minR, float maxR,
          bool filter,
          float ctE, float ctMinT, float ctMaxT, bool usecc, float clusterDt);

      void flagCrossTalk(std::unique_ptr<StrawHitCollection> const& shCol,
          std::unique_ptr<ComboHitCollection> const& chCol) const;

      bool createComboHit(EventWindowMarker const& ewm, size_t isd, std::unique_ptr<ComboHitCollection> const& chCol,
          std::unique_ptr<StrawHitCollection> const& shCol,
          const CaloClusterCollection *caloClusters,
          double pbtOffset,
          StrawId const& sid, TrkTypes::TDCValues const& tdc, TrkTypes::TOTValues const& tot,
          double pmp,
          TrackerStatus const& trackerStatus, StrawResponse const& srep, Tracker const& tt) const;

      double peakMinusPedWF(TrkTypes::ADCWaveform const& adcData, StrawResponse const& srep, ADCWFIter& maxiter) const;

    private:
      // algorith parameters
      TrkHitReco::FitType _fittype;
      int _diagLevel;
      StrawIdMask _mask;
      bool _writesh;
      // flag/filter parameters
      float _minT, _maxT;
      bool _overrideminTOff;
      float _minTOff;
      bool _overridemaxTOff;
      float _maxTOff;
      float _minE, _maxE, _minR, _maxR;
      bool _filter;
      // cross-talk parameters
      float _ctE, _ctMinT, _ctMaxT;
      // Calo-cluster parameters
      bool _usecc;
      float _clusterDt;

      std::vector<std::vector<size_t> > hits_by_panel;
      std::vector<size_t> largeHits;
      std::vector<size_t> largeHitPanels;
  };
}
#endif
