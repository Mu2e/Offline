//
// This module transforms StrawDigi objects into StrawHit objects
//
// Original author David Brown, LBNL
// Merged with flag and position creation B. Echenard, CalTech
//
// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

// conditions
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsBase/inc/TrackerCalibrationStructs.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/TrackerStatus.hh"

#include "TrkHitReco/inc/PeakFit.hh"
#include "TrkHitReco/inc/PeakFitRoot.hh"
#include "TrkHitReco/inc/PeakFitFunction.hh"
#include "TrkHitReco/inc/ComboPeakFitRoot.hh"

#include "RecoDataProducts/inc/ProtonBunchTime.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

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
	fhicl::Atom<float>minT{ Name("minimumTime"), Comment("Earliest StrawDigi time to process (nsec)"),500};
	fhicl::Atom<float>maxT{ Name("maximumTime"), Comment("Latest StrawDigi time to process (nsec)"),2000};
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
      float peakMinusPedAvg(TrkTypes::ADCWaveform const& adcData) const;
      float peakMinusPed(StrawId id, TrkTypes::ADCWaveform const& adcData) const;
      float peakMinusPedFirmware(StrawId id, TrkTypes::ADCValue const& pmp) const;
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

      std::vector<std::vector<size_t> > hits_by_panel(nplanes*npanels,std::vector<size_t>());
      std::vector<size_t> largeHits, largeHitPanels;
      largeHits.reserve(sdcol.size());
      largeHitPanels.reserve(sdcol.size());

      TrackerStatus const& trackerStatus = _trackerStatus_h.get(event.id());

      for (size_t isd=0;isd<sdcol.size();++isd) {
	const StrawDigi& digi = sdcol[isd];
	// flag digis that shouldn't be here or we don't want
	StrawHitFlag flag;
	if (trackerStatus.noSignal(digi.strawId()) || trackerStatus.suppress(digi.strawId())) {
	  flag.merge(StrawHitFlag::dead);
	}

	// start by reconstructing the times
	TDCTimes times;
	srep.calibrateTimes(digi.TDC(),times,digi.strawId());
	// find the end with the earliest time
	StrawEnd eend(StrawEnd::cal);
	if(times[StrawEnd::hv] < times[StrawEnd::cal])
	  eend = StrawEnd(StrawEnd::hv);
	// take the earliest of the 2 end times
	float time = times[eend.end()] - pbtOffset;
	if (time < _minT || time > _maxT ){
	  if(_filter)continue;
	} else
	  flag.merge(StrawHitFlag::timesel);

	//calorimeter filtering
	if (_usecc && caloClusters) {
	  bool outsideCaloTime(true);
	  for (const auto& cluster : *caloClusters)
	    if (std::abs(time-cluster.time())<_clusterDt) {outsideCaloTime=false; break;}
	  if (outsideCaloTime){
	    if(_filter)continue;
	  } else
	    flag.merge(StrawHitFlag::calosel);
	}
	// filter based on waveform shape (xtalk, undershoot, etc).  FIXME!
	//extract energy from waveform
	float energy(0.0);
	if (_fittype == TrkHitReco::FitType::peakminuspedavg){
          auto adcwaveform = sdadccol->at(isd);
	  float charge = peakMinusPedAvg(adcwaveform.samples());
	  energy = srep.ionizationEnergy(charge);
	} else if (_fittype == TrkHitReco::FitType::peakminusped){
          auto adcwaveform = sdadccol->at(isd);
	  float charge = peakMinusPed(digi.strawId(),adcwaveform.samples());
	  energy = srep.ionizationEnergy(charge);
	} else if (_fittype == TrkHitReco::FitType::firmwarepmp){
          float charge = peakMinusPedFirmware(digi.strawId(), digi.PMP());
          energy = srep.ionizationEnergy(charge);
        } else {
          auto adcwaveform = sdadccol->at(isd);
	  TrkHitReco::PeakFitParams params;
	  _pfit->process(adcwaveform.samples(),params);
	  energy = srep.ionizationEnergy(params._charge/srep.strawGain());
	  if (_printLevel > 1) std::cout << "Fit status = " << params._status << " NDF = " << params._ndf << " chisquared " << params._chi2
	    << " Fit charge = " << params._charge << " Fit time = " << params._time << std::endl;
	}
	// energy selection
	if( energy > _maxE || energy < _minE ) {
	  if(_filter) continue;
	} else
	  flag.merge(StrawHitFlag::energysel);
	// time-over-threshold
	TOTTimes tots{0.0,0.0};
	for(size_t iend=0;iend<2;++iend){
	  tots[iend] = digi.TOT(_end[iend])*srep.totLSB();
	}
	// choose earliest end TOT: maybe average later?
	float tot = tots[eend.end()];
	// filter on specific ionization FIXME!
	// filter based on composite e/P separation FIXME!
	const Straw& straw  = tt.getStraw( digi.strawId() );
	double dw, dwerr;
	double dt = times[StrawEnd::hv] - times[StrawEnd::cal];
        double halfpv;
	// get distance along wire from the straw center and it's estimated error
	bool td = srep.wireDistance(straw,energy,dt, dw,dwerr,halfpv);
        float propd = straw.halfLength()+dw;
        if (eend == StrawEnd(StrawEnd::cal))
          propd = straw.halfLength()-dw;
	XYZVec pos = Geom::toXYZVec(straw.getMidPoint()+dw*straw.getDirection());
	// create combo hit
	static const XYZVec _zdir(0.0,0.0,1.0);
	ComboHit ch;
	ch._nsh = 1; // 'combo' of 1 hit
	ch._pos = pos;
	ch._wdir = straw.getDirection();
	ch._sdir = _zdir.Cross(ch._wdir);
	ch._wdist = dw;
	ch._wres = dwerr;
	ch._time = time;
	ch._edep = energy;
	ch._sid = straw.id();
	ch._dtime = srep.driftTime(straw,tot,energy);
        ch._ptime = propd/(2*halfpv);
	ch._pathlength = srep.pathLength(straw,tot);
	ch.addIndex(isd); // reference the digi; this allows MC truth matching to work
	// crude initial estimate of the transverse error
	static const float invsqrt12 = 1.0/sqrt(12.0);
	ch._tres = tt.strawOuterRadius()*invsqrt12;
	// set flags
	ch._mask = _mask;
	ch._flag = flag;
	if (td) ch._flag.merge(StrawHitFlag::tdiv);
	ch._tend = eend;
	if(!_filter && _flagXT){
	  //buffer large hit for cross-talk analysis
	  size_t iplane       = straw.id().getPlane();
	  size_t ipnl         = straw.id().getPanel();
	  size_t global_panel = ipnl + iplane*npanels;
	  hits_by_panel[global_panel].push_back(shCol->size());
	  if (energy >= _ctE) {largeHits.push_back(shCol->size()); largeHitPanels.push_back(global_panel);}
	}
	chCol->push_back(std::move(ch));
	// optionally create legacy straw hit (for diagnostics and calibration)
	if(_writesh){
	  StrawHit hit(digi.strawId(),times,tots,energy);
	  shCol->push_back(std::move(hit));
	}
      }
      //flag straw and electronic cross-talk
      if(!_filter && _flagXT){
	for (size_t ilarge=0; ilarge < largeHits.size();++ilarge)
	{
	  const StrawHit& sh = (*shCol)[largeHits[ilarge]];
	  for (size_t jsh : hits_by_panel[largeHitPanels[ilarge]])
	  {
            if (jsh==largeHits[ilarge]) continue;
            const StrawHit& sh2 = (*shCol)[jsh];
            if (sh2.time()-sh.time() > _ctMinT && sh2.time()-sh.time() < _ctMaxT)
            {
              if (sh.strawId().samePreamp(sh2.strawId())) (*chCol)[jsh]._flag.merge(StrawHitFlag::elecxtalk);
              if (sh.strawId().nearestNeighbor(sh2.strawId())) (*chCol)[jsh]._flag.merge(StrawHitFlag::strawxtalk);
            }
          }
        }
      }

      if(_writesh)event.put(std::move(shCol));
      event.put(std::move(chCol));
  }

  float StrawHitReco::peakMinusPedAvg(TrkTypes::ADCWaveform const& adcData) const {
    auto wfstart = adcData.begin() + _npre;
    float pedestal = std::accumulate(adcData.begin(), wfstart, 0)*_invnpre;
//    auto maxIter = std::max_element(wfstart,adcData.end());
    auto maxIter = wfstart;
    auto nextIter = maxIter; nextIter++;
    while(nextIter != adcData.end() && *nextIter > *maxIter){
      ++maxIter;
      ++nextIter;
    }
    float peak = *maxIter;
    if(_diagLevel > 0)_maxiter->Fill(std::distance(wfstart,maxIter));
    return (peak-pedestal)*_invgainAvg;
  }

  float StrawHitReco::peakMinusPed(StrawId id, TrkTypes::ADCWaveform const& adcData) const {
    auto wfstart = adcData.begin() + _npre;
    float pedestal = std::accumulate(adcData.begin(), wfstart, 0)*_invnpre;
//    auto maxIter = std::max_element(wfstart,adcData.end());
    auto maxIter = wfstart;
    while(maxIter != adcData.end() && *(maxIter+1) > *maxIter)
      ++maxIter;
    float peak = *maxIter;
    if(_diagLevel > 0)_maxiter->Fill(std::distance(wfstart,maxIter));
    return (peak-pedestal)*_invgain[id.getStraw()];
  }

  float StrawHitReco::peakMinusPedFirmware(StrawId id, TrkTypes::ADCValue const& pmp) const {
    return pmp*_invgain[id.getStraw()];
  }



}


using mu2e::StrawHitReco;
DEFINE_ART_MODULE(StrawHitReco);
