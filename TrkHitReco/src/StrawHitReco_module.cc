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
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Run.h"

// conditions
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ConditionsService/inc/AcceleratorParams.hh"
#include "Offline/ConditionsBase/inc/TrackerCalibrationStructs.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"

#include "Offline/TrkHitReco/inc/PeakFit.hh"
#include "Offline/TrkHitReco/inc/PeakFitRoot.hh"
#include "Offline/TrkHitReco/inc/PeakFitFunction.hh"
#include "Offline/TrkHitReco/inc/ComboPeakFitRoot.hh"
#include "Offline/TrkHitReco/inc/StrawHitRecoUtils.hh"

#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"


#include "TH1F.h"

#include <memory>
#include <numeric>


namespace mu2e {
  using namespace TrkTypes;

  class StrawHitReco : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int> diag{ Name("diagLevel"), Comment("Diag level"), 0};
        fhicl::Atom<int> print{ Name("printLevel"), Comment("Print level"), 0};
        fhicl::Atom<int> fittype { Name( "FitType"), Comment("Waveform Fit Type")};//, TrkHitReco::FitType::firmwarepmp};
        fhicl::Atom<bool> usecc{ Name("UseCalorimeter"), Comment("Use Calo cluster times to filter" ),false};
        fhicl::Atom<float>clusterDt{ Name("clusterDt"), Comment("Calo cluster time 1/2 window"),100};
        fhicl::Atom<float>minE{ Name("minimumEnergy"), Comment("Minimum straw energy deposit (MeV)"),0.0};
        fhicl::Atom<float>maxE{ Name("maximumEnergy"), Comment("Maximum straw energy deposit (MeV)"),0.0035};
        fhicl::Atom<float>ctE{ Name("crossTalkEnergy"), Comment("Energy to filter cross-talk in adajcent straws (MeV)"),0.007};
        fhicl::Atom<float>ctMinT{ Name("crossTalkMinimumTime"), Comment("Earliest time for cross-talk filter (nsec)"),-1};
        fhicl::Atom<float>ctMaxT{ Name("crossTalkMaximumTime"), Comment("Latest time for cross-talk filter (nsec)"),100};
        fhicl::Atom<float>minT{ Name("minimumTime"), Comment("Earliest StrawDigi time to process (nsec)")};
        fhicl::Atom<float>maxT{ Name("maximumTime"), Comment("Latest StrawDigi time to process (nsec)")};
        fhicl::Atom<bool>filter{ Name("FilterHits"), Comment("Filter hits (alternative is to just flag)") };
        fhicl::Atom<bool>writesh{ Name("WriteStrawHitCollection"), Comment("Save StrawHitCollection")};
        fhicl::Atom<bool>flagXT{ Name("FlagCrossTalk"), Comment("Search for cross-talk"),false};
        fhicl::Atom<art::InputTag> sdcTag{ Name("StrawDigiCollectionTag"), Comment("StrawDigiCollection producer")};
        fhicl::Atom<art::InputTag> sdadcTag{ Name("StrawDigiADCWaveformCollectionTag"), Comment("StrawDigiADCWaveformCollection producer")};
        fhicl::Atom<art::InputTag> cccTag{ Name("CaloClusterCollectionTag"), Comment("CaloClusterCollection producer")};
        fhicl::Atom<art::InputTag> pbttoken{ Name("ProtonBunchTimeTag"), Comment("ProtonBunchTime producer")};
  };

  using Parameters = art::EDProducer::Table<Config>;
  explicit StrawHitReco(Parameters const& config);
  void produce( art::Event& e) override;
  void beginRun( art::Run& run ) override;
  void beginJob() override;


  private:
  TrkHitReco::FitType _fittype; // peak Fitter
  bool  _usecc;                   // use calorimeter cluster filtering
  float _clusterDt;               // maximum hit-calo lcuster time difference
  float _minE;             // energy range (MeV)
  float _maxE;             // energy range (MeV)
  float _ctE;                     // minimum charge to flag neighbors as cross talk
  float _ctMinT;                  // time relative to proton hit to flag cross talk (ns)
  float _ctMaxT;                  // time relative to proton hit to flag cross talk (ns)
  float _minT, _maxT;             // time range
  bool  _filter;                // filter the output, or just flag
  bool  _writesh;                // write straw hits or not
  bool  _flagXT; // flag cross-talk
  int   _printLevel;
  int   _diagLevel;
  StrawIdMask _mask;
  StrawEnd _end[2]; // helper
  float _invnpre; // cache
  float _invgainAvg; // cache
  float _invgain[96]; // cache
  unsigned _npre; //cache

  art::ProductToken<StrawDigiCollection> const _sdctoken;
  art::ProductToken<StrawDigiADCWaveformCollection> const _sdadctoken;
  art::ProductToken<CaloClusterCollection> const _ccctoken;
  art::ProductToken<ProtonBunchTime> const _pbttoken; // name of the module that makes eventwindowmarkers
  std::unique_ptr<TrkHitReco::PeakFit> _pfit; // peak fitting algorithm
  // diagnostic
  TH1F* _maxiter;
  // helper function
  ProditionsHandle<StrawResponse> _strawResponse_h;
  ProditionsHandle<TrackerStatus> _trackerStatus_h;
  ProditionsHandle<Tracker> _alignedTracker_h;
};

StrawHitReco::StrawHitReco(Parameters const& config) :
  art::EDProducer{config},
  _fittype((TrkHitReco::FitType) config().fittype()),
  _usecc(config().usecc()),
  _clusterDt(config().clusterDt()),
  _minE(config().minE()),
  _maxE(config().maxE()),
  _ctE(config().ctE()),
  _ctMinT(config().ctMinT()),
  _ctMaxT(config().ctMaxT()),
  _minT(config().minT()),
  _maxT(config().maxT()),
  _filter(config().filter()),
  _writesh(config().writesh()),
  _flagXT(config().flagXT()),
  _printLevel(config().print()),
  _diagLevel(config().diag()),
  _mask(StrawIdMask::uniquestraw), // this module produces individual straw ComboHits
  _end{StrawEnd::cal,StrawEnd::hv}, // this should be in a general place, FIXME!
  _sdctoken{consumes<StrawDigiCollection>(config().sdcTag())},
  _sdadctoken{mayConsume<StrawDigiADCWaveformCollection>(config().sdadcTag())},
  _ccctoken{mayConsume<CaloClusterCollection>(config().cccTag())},
  _pbttoken{consumes<ProtonBunchTime>(config().pbttoken())}
{
  produces<ComboHitCollection>();
  if(_writesh)produces<StrawHitCollection>();
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

void StrawHitReco::beginRun(art::Run& run)
{
  auto const& srep = _strawResponse_h.get(run.id());
  // set cache for peak-ped calculation (default)
  _npre = srep.nADCPreSamples();
  _invnpre = 1.0/(float)_npre;
  _invgainAvg = srep.adcLSB()*srep.peakMinusPedestalEnergyScale()/srep.strawGain();
  for (int i=0;i<96;i++){
    StrawId dummyId(0,0,i);
    _invgain[i] = srep.adcLSB()*srep.peakMinusPedestalEnergyScale(dummyId)/srep.strawGain();
  }

  // Detailed histogram-based waveform fits are no longer supported TODO!
  if (_fittype != TrkHitReco::FitType::peakminusped && _fittype != TrkHitReco::FitType::peakminuspedavg && _fittype != TrkHitReco::FitType::firmwarepmp)
    throw cet::exception("RECO")<<"TrkHitReco: Peak fit " << _fittype << " not implemented " <<  std::endl;
}

//------------------------------------------------------------------------------------------
void StrawHitReco::produce(art::Event& event)
{
  if (_printLevel > 0) std::cout << "In StrawHitReco produce " << std::endl;

  const Tracker& tt = _alignedTracker_h.get(event.id());

  size_t nplanes = tt.nPlanes();
  size_t npanels = tt.getPlane(0).nPanels();
  auto const& srep = _strawResponse_h.get(event.id());
  auto sdH = event.getValidHandle(_sdctoken);
  const StrawDigiCollection& sdcol(*sdH);

  const StrawDigiADCWaveformCollection *sdadccol(0);
  if (_fittype != TrkHitReco::FitType::firmwarepmp) {
    auto sdawH = event.getValidHandle(_sdadctoken);
    sdadccol = sdawH.product();
  }
  StrawDigiADCWaveform adcwaveform;

  const CaloClusterCollection* caloClusters(0);
  if(_usecc){
    auto ccH = event.getValidHandle(_ccctoken);
    caloClusters = ccH.product();
  }

  double pbtOffset = 0;
  art::Handle<ProtonBunchTime> pbtHandle;
  auto pbtH = event.getValidHandle(_pbttoken);
  const ProtonBunchTime& pbt(*pbtH);
  pbtOffset = pbt.pbtime_;

  std::unique_ptr<StrawHitCollection> shCol;
  if(_writesh){
    shCol = std::unique_ptr<StrawHitCollection>(new StrawHitCollection);
    shCol->reserve(sdcol.size());
  }
  std::unique_ptr<ComboHitCollection> chCol(new ComboHitCollection());
  chCol->reserve(sdcol.size());



  TrackerStatus const& trackerStatus = _trackerStatus_h.get(event.id());

  mu2e::StrawHitRecoUtils shrUtils(pbtOffset, _fittype, _npre, _invnpre, _invgainAvg, _invgain,
      _diagLevel, _maxiter, _mask, nplanes, npanels, _writesh, _minT, _maxT, _minE, _maxE, _filter, _flagXT,
      _ctE, _ctMinT, _ctMaxT, _usecc, _clusterDt, sdcol.size());

  for (size_t isd=0;isd<sdcol.size();++isd) {
    const StrawDigi& digi = sdcol[isd];

    if (_fittype != TrkHitReco::FitType::firmwarepmp)
      adcwaveform = sdadccol->at(isd);

    shrUtils.createComboHit(isd, chCol, shCol, caloClusters, digi.strawId(), digi.TDC(), digi.TOT(), digi.PMP(), adcwaveform.samples(),
        trackerStatus,  srep, tt);
  }
  //flag straw and electronic cross-talk
  if(!_filter && _flagXT){
    shrUtils.flagCrossTalk(shCol, chCol);
  }

  if(_writesh)event.put(std::move(shCol));
  event.put(std::move(chCol));
}

}


using mu2e::StrawHitReco;
DEFINE_ART_MODULE(StrawHitReco);
