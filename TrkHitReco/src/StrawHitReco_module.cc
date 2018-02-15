//
// This module transforms StrawDigi objects into StrawHit objects
//
// $Id: StrawHitReco_module.cc,v 1.12 2014/03/25 22:14:39 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/03/25 22:14:39 $
//
// Original author David Brown, LBNL
// Merged with flag and position creation B. Echenard, CalTech
//
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "ConditionsBase/inc/TrackerCalibrationStructs.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "TrackerConditions/inc/StrawResponse.hh"

#include "TrkHitReco/inc/PeakFit.hh"
#include "TrkHitReco/inc/PeakFitRoot.hh"
#include "TrkHitReco/inc/PeakFitFunction.hh"
#include "TrkHitReco/inc/ComboPeakFitRoot.hh"

#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

#include <memory>

 

namespace mu2e {
  using namespace TrkTypes;

  class StrawHitReco : public art::EDProducer 
  {
     public:
       explicit StrawHitReco(fhicl::ParameterSet const& pset);
       virtual ~StrawHitReco(); 
       virtual void produce( art::Event& e);
       virtual void beginRun( art::Run& run );
       virtual void beginJob();


     private:
       TrkHitReco::FitType _fittype; // peak Fitter
       bool   _usecc;                   // use calorimeter cluster filtering
       float _clusterDt;               // maximum hit-calo lcuster time difference
       float _minE, _maxE;             // energy range (MeV)
       float _ctE;                     // minimum charge to flag neighbors as cross talk
       float _ctMinT;                  // time relative to proton hit to flag cross talk (ns)
       float _ctMaxT;                  // time relative to proton hit to flag cross talk (ns)
       float _minT, _maxT;             // time range
       float _maxR2;		      // maximum transverse radius (squared) 
       bool   _filter;                // filter the output, or just flag 
       bool  _writesh;		      // write straw hits or not
       bool _flagXT; // flag cross-talk
       int    _printLevel;
       int    _diagLevel;
       StrawIdMask _mask; 
       StrawEnd _end[2]; // helper

       art::InputTag _sdtag, _cctag;
       fhicl::ParameterSet _peakfit;  // peak fit (charge reconstruction) parameters
       std::unique_ptr<TrkHitReco::PeakFit> _pfit; // peak fitting algorithm
       
 };

  StrawHitReco::StrawHitReco(fhicl::ParameterSet const& pset) :
      _fittype((TrkHitReco::FitType) pset.get<unsigned>("FitType",TrkHitReco::FitType::peakminusped)),
      _usecc(pset.get<bool>(         "UseCalorimeter",false)),     
      _clusterDt(pset.get<float>(   "clusterDt",100)),
      _minE(pset.get<float>(        "minimumEnergy",0.0)), // Minimum deposited straw energy (MeV)
      _maxE(pset.get<float>(        "maximumEnergy",0.003)), // MeV
      _ctE(pset.get<float>(         "crossTalkEnergy",0.007)), // MeV
      _ctMinT(pset.get<float>(      "crossTalkMinimumTime",-1)), // nsec
      _ctMaxT(pset.get<float>(      "crossTalkMaximumTime",100)), // nsec
      _minT(pset.get<float>(        "minimumTime",500)), // nsec
      _maxT(pset.get<float>(        "maximumTime",2000)), // nsec
      _filter(pset.get<bool>(      "FilterHits")),
      _writesh(pset.get<bool>(      "WriteStrawHitCollection")),
      _flagXT(pset.get<bool>(      "FlagCrossTalk",false)),
      _printLevel(pset.get<int>(     "printLevel",0)),
      _diagLevel(pset.get<int>(      "diagLevel",0)),
      _end{TrkTypes::cal,TrkTypes::hv}, // this should be in a general place, FIXME!
      _sdtag  (pset.get<art::InputTag>("StrawDigiCollection","makeSD")),
      _cctag  (pset.get<art::InputTag>("caloClusterModuleLabel","CaloClusterFast")),
      _peakfit(pset.get<fhicl::ParameterSet>("PeakFitter",fhicl::ParameterSet()))
  {
      produces<ComboHitCollection>();
      if(_writesh)produces<StrawHitCollection>();
      // each hit is a unique straw
      std::vector<StrawIdMask::field> fields{StrawIdMask::plane,StrawIdMask::panel,StrawIdMask::straw};
      _mask = StrawIdMask(fields);
      float maxR = pset.get<float>("maximumRadius",650); // mm
      _maxR2 = maxR*maxR;

      if (_printLevel > 0) std::cout << "In StrawHitReco constructor " << std::endl;
  }

  StrawHitReco::~StrawHitReco() {}

  //------------------------------------------------------------------------------------------
  void StrawHitReco::beginJob()
  {
  }

  void StrawHitReco::beginRun(art::Run& run)
  {    
      ConditionsHandle<StrawElectronics> strawele = ConditionsHandle<StrawElectronics>("ignored");
      // this must be done here because strawele is not accessible at startup and pfit references it
      if (_fittype == TrkHitReco::FitType::peakminusped)
         _pfit = std::unique_ptr<TrkHitReco::PeakFit>(new TrkHitReco::PeakFit(*strawele,_peakfit) );
      else if (_fittype == TrkHitReco::FitType::combopeakfit)
	 _pfit = std::unique_ptr<TrkHitReco::PeakFit>(new TrkHitReco::ComboPeakFitRoot(*strawele,_peakfit) );
      else
	 _pfit = std::unique_ptr<TrkHitReco::PeakFit>(new TrkHitReco::PeakFitRoot(*strawele,_peakfit) );
      if (_printLevel > 0) std::cout << "In StrawHitReco begin Run " << std::endl;
  }

  //------------------------------------------------------------------------------------------
  void StrawHitReco::produce(art::Event& event)
  {        
      if (_printLevel > 0) std::cout << "In StrawHitReco produce " << std::endl;

      const Tracker& tracker = getTrackerOrThrow();
      const TTracker& tt(*GeomHandle<TTracker>());
      size_t nplanes = tt.nPlanes();
      size_t npanels = tt.getPlane(0).nPanels();
      
      ConditionsHandle<StrawElectronics> strawele = ConditionsHandle<StrawElectronics>("ignored");
      ConditionsHandle<StrawPhysics> strawphys = ConditionsHandle<StrawPhysics>("ignored");
      ConditionsHandle<StrawResponse> srep = ConditionsHandle<StrawResponse>("ignored");
      auto sdH = event.getValidHandle<StrawDigiCollection>(_sdtag);
      const StrawDigiCollection& sdcol(*sdH);
      
      const CaloClusterCollection* caloClusters(0);
      if(_usecc){
	auto ccH = event.getValidHandle<CaloClusterCollection>(_cctag);
	caloClusters = ccH.product();
      }

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
      
      for (size_t isd=0;isd<sdcol.size();++isd)
      {
	const StrawDigi& digi = sdcol[isd];
	StrawHitFlag flag;
	// start by reconstructing the times
	TDCTimes times;
	strawele->tdcTimes(digi.TDC(),times);
	// take the earliest of the 2 end times
	float time = std::min(times[0],times[1]);
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

	//extract energy from waveform
	TrkHitReco::PeakFitParams params;
	_pfit->process(digi.adcWaveform(),params);
	float energy = strawphys->ionizationEnergy(params._charge/strawphys->strawGain());
	if (_printLevel > 1) std::cout << "Fit status = " << params._status << " NDF = " << params._ndf << " chisquared " << params._chi2
	  << " Fit charge = " << params._charge << " Fit time = " << params._time << std::endl;

	if( energy < _minE  || energy > _maxE){
	  if(_filter) continue;
	} else
	  flag.merge(StrawHitFlag::energysel);
	// time-over-threshold
	TOTTimes tots{0.0,0.0};
	for(size_t iend=0;iend<2;++iend){
	  tots[iend] = digi.TOT(_end[iend])*strawele->totLSB();
	}
	//create straw hit; this is currently required by the wireDistance function FIXME!
	StrawHit hit(tt.getStrawIndex(digi.strawId()),times,tots,energy);
	// get distance along wire from the straw center and it's estimated error
	const Straw& straw  = tracker.getStraw( digi.strawId() );
	float dw, dwerr;
	bool td = srep->wireDistance(hit,straw.getHalfLength(),dw,dwerr);
	XYZVec pos = Geom::toXYZVec(straw.getMidPoint()+dw*straw.getDirection());
	float rad2 = pos.Perp2();
	if(rad2 > _maxR2) {
	  if(_filter)continue;
	} else
	  flag.merge(StrawHitFlag::radsel);
	if(_writesh){
	  shCol->push_back(std::move(hit));
	}
	// create combo hit
	ComboHit ch;
	ch._nsh = 1; // 'combo' of 1 hit
	ch._pos = pos;
	ch._wdir = straw.getDirection();
	ch._wdist = dw;
	ch._wres = dwerr;
	ch._time = time;
	ch._edep = energy;
	ch._sid = straw.id();
	ch.addIndex(isd); // reference the digi; this allows MC truth matching to work
	// crude initial estimate of the transverse error
	static const float invsqrt12 = 1.0/sqrt(12.0);
	ch._tres = straw.getRadius()*invsqrt12;
	// set flags
	ch._mask = _mask;
	ch._flag = flag;
	if (td) ch._flag.merge(StrawHitFlag::tdiv); 
	if(!_filter && _flagXT){
	  //buffer large hit for cross-talk analysis
	  size_t iplane       = straw.id().getPlane();
	  size_t ipnl         = straw.id().getPanel();
	  size_t global_panel = ipnl + iplane*npanels;
	  hits_by_panel[global_panel].push_back(shCol->size());          
	  if (energy >= _ctE) {largeHits.push_back(shCol->size()); largeHitPanels.push_back(global_panel);}
	}

	chCol->push_back(std::move(ch));
      }

      //flag straw and electronic cross-talk
      if(!_filter && _flagXT){
	for (size_t ilarge=0; ilarge < largeHits.size();++ilarge)
	{
	  const StrawHit& sh = (*shCol)[largeHits[ilarge]];
	  const Straw& straw = tracker.getStraw( sh.strawIndex() );
	  for (size_t jsh : hits_by_panel[largeHitPanels[ilarge]])
	  {
	    if (jsh==largeHits[ilarge]) continue;              
	    const StrawHit& sh2 = (*shCol)[jsh]; 
	    if (sh2.time()-sh.time() > _ctMinT && sh2.time()-sh.time() < _ctMaxT)
	    {
	      if (straw.isSamePreamp(sh2.strawIndex()))       (*chCol)[jsh]._flag.merge(StrawHitFlag::elecxtalk);
	      if (straw.isNearestNeighbour(sh2.strawIndex())) (*chCol)[jsh]._flag.merge(StrawHitFlag::strawxtalk);
	    }           
	  }
	}
      }

      if(_writesh)event.put(std::move(shCol));
      event.put(std::move(chCol));
  }

}


using mu2e::StrawHitReco;
DEFINE_ART_MODULE(StrawHitReco);

