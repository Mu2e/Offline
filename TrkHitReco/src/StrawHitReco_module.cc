//
// This module transforms StrawDigi objects into StrawHit objects
//
// Original author David Brown, LBNL
// Merged with flag and position creation B. Echenard, CalTech
//
// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Run.h"

// conditions
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/ConditionsBase/inc/TrackerCalibrationStructs.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"

#include "Offline/TrkHitReco/inc/PeakFit.hh"
#include "Offline/TrkHitReco/inc/StrawHitRecoUtils.hh"

#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

#include "TH1F.h"

#include <memory>
#include <numeric>


namespace mu2e {
  using namespace TrkTypes;

  class StrawHitReco : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      using ADCWFIter = TrkTypes::ADCWaveform::const_iterator;
      struct Config {
        fhicl::Atom<int> diag{ Name("diagLevel"), Comment("Diag level"), 0};
        fhicl::Atom<int> print{ Name("printLevel"), Comment("Print level"), 0};
        fhicl::Atom<int> fittype { Name( "FitType"), Comment("Waveform Fit Type")};
        fhicl::Atom<bool> usecc{ Name("UseCalorimeter"), Comment("Use Calo cluster times to filter" )};
        fhicl::Atom<float>clusterDt{ Name("clusterDt"), Comment("Calo cluster time 1/2 window")};
        fhicl::Atom<float>minE{ Name("MinimumEnergy"), Comment("Minimum straw energy deposit (MeV)")};
        fhicl::Atom<float>maxE{ Name("MaximumEnergy"), Comment("Maximum straw energy deposit (MeV)")};
        fhicl::Atom<float>minR{ Name("MinimumRadius"), Comment("Minimum transverse radius (mm)")};
        fhicl::Atom<float>maxR{ Name("MaximumRadius"), Comment("Maximum transverse radius (mm)")};
        fhicl::Atom<float>minT{ Name("MinimumTime"), Comment("Earliest StrawDigi time to process (nsec)")};
        fhicl::Atom<float>maxT{ Name("MaximumTime"), Comment("Latest StrawDigi time to process (nsec)")};
        fhicl::OptionalAtom<float>minTOff{ Name("MinimumTimeOffSpill"), Comment("Earliest StrawDigi time to process (nsec)")};
        fhicl::OptionalAtom<float>maxTOff{ Name("MaximumTimeOffSpill"), Comment("Latest StrawDigi time to process (nsec)")};
        fhicl::Atom<float>ctE{ Name("crossTalkEnergy"), Comment("Energy to filter cross-talk in adajcent straws (MeV)")};
        fhicl::Atom<float>ctMinT{ Name("crossTalkMinimumTime"), Comment("Earliest time for cross-talk filter (nsec)")};
        fhicl::Atom<float>ctMaxT{ Name("crossTalkMaximumTime"), Comment("Latest time for cross-talk filter (nsec)")};
        fhicl::Atom<bool>filter{ Name("FilterHits"), Comment("Filter hits (alternative is to just flag)") };
        fhicl::Atom<bool>writesh{ Name("WriteStrawHitCollection"), Comment("Save StrawHitCollection")};
        fhicl::Atom<bool>flagXT{ Name("FlagCrossTalk"), Comment("Search for cross-talk")};
        fhicl::Atom<art::InputTag> sdcTag{ Name("StrawDigiCollectionTag"), Comment("StrawDigiCollection producer")};
        fhicl::Atom<art::InputTag> sdadcTag{ Name("StrawDigiADCWaveformCollectionTag"), Comment("StrawDigiADCWaveformCollection producer")};
        fhicl::Atom<art::InputTag> cccTag{ Name("CaloClusterCollectionTag"), Comment("CaloClusterCollection producer")};
        fhicl::Atom<art::InputTag> pbttoken{ Name("ProtonBunchTimeTag"), Comment("ProtonBunchTime producer")};
        fhicl::Atom<art::InputTag> EWM { Name("EventWindowMarker"), Comment("EventWindowMarker")};
     };

      using Parameters = art::EDProducer::Table<Config>;
      explicit StrawHitReco(Parameters const& config);
      void produce( art::Event& e) override;
      void beginJob() override;

    private:
      float _minTOff, _maxTOff;
      bool _overrideminTOff;
      bool _overridemaxTOff;
      StrawHitRecoUtils _shrUtils;
      bool  _writesh;                // write straw hits or not
      bool  _flagXT; // flag cross-talk
      bool _usecc;
      bool  _useADCWF;
      int   _printLevel;
      int   _diagLevel;

      art::ProductToken<StrawDigiCollection> const _sdctoken;
      art::ProductToken<StrawDigiADCWaveformCollection> const _sdadctoken;
      art::ProductToken<CaloClusterCollection> const _ccctoken;
      art::ProductToken<ProtonBunchTime> const _pbttoken;
      art::ProductToken<EventWindowMarker> const _ewmtoken;
      std::unique_ptr<TrkHitReco::PeakFit> _pfit; // peak fitting algorithm
      // diagnostic
      TH1F* _maxiter;
      // handles
      ProditionsHandle<StrawResponse> _strawResponse_h;
      ProditionsHandle<TrackerStatus> _trackerStatus_h;
      ProditionsHandle<Tracker> _alignedTracker_h;
  };

  StrawHitReco::StrawHitReco(Parameters const& config) :
    art::EDProducer{config},
    _minTOff(config().minT()),
    _maxTOff(config().maxT()),
    _overrideminTOff(config().minTOff(_minTOff)),
    _overridemaxTOff(config().maxTOff(_maxTOff)),
    _shrUtils ((TrkHitReco::FitType) config().fittype(),
        config().diag(),
        StrawIdMask::uniquestraw, // this module produces individual straw ComboHits
        config().writesh(),
        config().minT(),
        config().maxT(),
        _overrideminTOff,
        _minTOff,
        _overridemaxTOff,
        _maxTOff,
        config().minE(),
        config().maxE(),
        config().minR(),
        config().maxR(),
        config().filter(),
        config().ctE(),
        config().ctMinT(),
        config().ctMaxT(),
        config().usecc(),
        config().clusterDt()),
    _writesh(config().writesh()),
    _flagXT(config().flagXT()),
    _usecc(config().usecc()),
    _useADCWF(config().fittype() != mu2e::TrkHitReco::FitType::firmwarepmp ),
    _printLevel(config().print()),
    _diagLevel(config().diag()),
    _sdctoken{consumes<StrawDigiCollection>(config().sdcTag())},
    _sdadctoken{mayConsume<StrawDigiADCWaveformCollection>(config().sdadcTag())},
    _ccctoken{mayConsume<CaloClusterCollection>(config().cccTag())},
    _pbttoken{consumes<ProtonBunchTime>(config().pbttoken())},
    _ewmtoken{consumes<EventWindowMarker>(config().EWM())}

  {
    produces<ComboHitCollection>();
    produces<IntensityInfoTrackerHits>();
    if (_writesh) produces<StrawHitCollection>();
    if (_printLevel > 0) std::cout << "In StrawHitReco constructor " << std::endl;
  }

  //------------------------------------------------------------------------------------------
  void StrawHitReco::beginJob()
  {
    if(_diagLevel > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _maxiter   = tfs->make<TH1F>( "maxiter",  "ADC max",16,-0.5,15.5 );
    }
  }

  //------------------------------------------------------------------------------------------
  void StrawHitReco::produce(art::Event& event)
  {
    if (_printLevel > 0) std::cout << "In StrawHitReco produce " << std::endl;

    const Tracker& tt = _alignedTracker_h.get(event.id());
    auto const& srep = _strawResponse_h.get(event.id());
    auto sdH = event.getValidHandle(_sdctoken);
    const StrawDigiCollection& sdcol(*sdH);

    const StrawDigiADCWaveformCollection* sdadcc(0);
    if (_useADCWF){
      auto sdawH = event.getValidHandle(_sdadctoken);
      sdadcc = sdawH.product();
    }

    const CaloClusterCollection* caloClusters(0);
    if(_usecc){
      auto ccH = event.getValidHandle(_ccctoken);
      caloClusters = ccH.product();
    }

    auto pbtH = event.getValidHandle(_pbttoken);
    const ProtonBunchTime& pbt(*pbtH);
    double pbtOffset = pbt.pbtime_;
    auto ewmH = event.getValidHandle(_ewmtoken);
    const EventWindowMarker& ewm(*ewmH);

    std::unique_ptr<StrawHitCollection> shCol;
    if(_writesh){
      shCol = std::unique_ptr<StrawHitCollection>(new StrawHitCollection);
      shCol->reserve(sdcol.size());
    }
    std::unique_ptr<ComboHitCollection> chCol(new ComboHitCollection());
    std::unique_ptr<IntensityInfoTrackerHits>  intInfo(new IntensityInfoTrackerHits());
    chCol->reserve(sdcol.size());

    TrackerStatus const& trackerStatus = _trackerStatus_h.get(event.id());

    double pmp(0.0);
    for (size_t isd=0;isd<sdcol.size();++isd) {
      const StrawDigi& digi = sdcol[isd];
      // compute peak-pedestal
      if(!_useADCWF){
        pmp = digi.PMP();
      } else {
        auto const& adcwf = sdadcc->at(isd).samples();
        ADCWFIter maxiter;
        pmp = _shrUtils.peakMinusPedWF(adcwf,srep,maxiter);
        if(_diagLevel > 0)_maxiter->Fill( std::distance(adcwf.begin(),maxiter));
      }
      _shrUtils.createComboHit(ewm, isd, chCol, shCol, caloClusters, pbtOffset,
          digi.strawId(), digi.TDC(), digi.TOT(), pmp,
          trackerStatus,  srep, tt);
      //flag straw and electronic cross-talk
      if(_flagXT){
        _shrUtils.flagCrossTalk(shCol, chCol);
      }
    }
    if(_writesh)event.put(std::move(shCol));
    intInfo->setNTrackerHits(chCol->size());
    event.put(std::move(intInfo));
    event.put(std::move(chCol));
  }
}

using mu2e::StrawHitReco;
DEFINE_ART_MODULE(StrawHitReco)
