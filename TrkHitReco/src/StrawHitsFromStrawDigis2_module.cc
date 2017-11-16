//
// This module transforms StrawDigi objects into StrawHit objects
//
// $Id: StrawHitsFromStrawDigis2_module.cc,v 1.12 2014/03/25 22:14:39 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/03/25 22:14:39 $
//
// Original author David Brown, LBNL
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

#include "TrkChargeReco/inc/PeakFit.hh"
#include "TrkChargeReco/inc/PeakFitRoot.hh"
#include "TrkChargeReco/inc/PeakFitFunction.hh"
#include "TrkChargeReco/inc/ComboPeakFitRoot.hh"

#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

#include <memory>

 

namespace mu2e {
  using namespace TrkTypes;

  class StrawHitsFromStrawDigis2 : public art::EDProducer 
  {
     public:
       explicit StrawHitsFromStrawDigis2(fhicl::ParameterSet const& pset);
       virtual ~StrawHitsFromStrawDigis2(); 
       virtual void produce( art::Event& e);
       virtual void beginRun( art::Run& run );
       virtual void beginJob();


     private:
       double mbbuffer_;                // buffer on that for ghost hits (wrapping)
       double maxdt_;                   // maximum time difference between end times
       bool   singledigi_;              // turn single-end digitizations into hits
       TrkChargeReco::FitType fittype_; // peak Fitter
       bool   usecc_;                   // use calorimeter cluster filtering
       double clusterDt_;               // maximum hit-calo lcuster time difference
       double minE_;                    // minimum charge (units??)
       double maxE_;                    // maximum charge (units??)
       double ctE_;                     // minimum charge to flag neighbors as cross talk
       double ctMinT_;                  // time relative to proton hit to flag cross talk (ns)
       double ctMaxT_;                  // time relative to proton hit to flag cross talk (ns)
       double minT_;                    // minimum hit time
       double maxT_;                    // maximum hit time
       bool   trigMode_;                // trigger mode cut on dt and other thing 
       int    printLevel_;
       int    diagLevel_;
       StrawEnd end_[2]; // helper
 
       std::string strawDigis_;
       std::string caloClusterModuleLabel_;
       fhicl::ParameterSet peakfit_;       
       
       std::unique_ptr<TrkChargeReco::PeakFit> pfit_;
       
 };


  StrawHitsFromStrawDigis2::StrawHitsFromStrawDigis2(fhicl::ParameterSet const& pset) :
      mbbuffer_(pset.get<double>(    "TimeBuffer",100.0)), 
      maxdt_(pset.get<double>(       "MaxTimeDifference",8.0)), 
      singledigi_(pset.get<bool>(    "UseSingleDigis",false)), 
      fittype_((TrkChargeReco::FitType) pset.get<unsigned>("FitType",TrkChargeReco::FitType::peakminusped)),
      usecc_(pset.get<bool>(         "UseCalorimeter",false)),     
      clusterDt_(pset.get<double>(   "clusterDt",100)),
      minE_(pset.get<double>(        "minimumEnergy",0.0)), // Minimum deposited straw energy (MeV)
      maxE_(pset.get<double>(        "maximumEnergy",0.01)), // MeV
      ctE_(pset.get<double>(         "crossTalkEnergy",0.007)), // MeV
      ctMinT_(pset.get<double>(      "crossTalkMinimumTime",-1)), // nsec
      ctMaxT_(pset.get<double>(      "crossTalkMaximumTime",100)), // nsec
      minT_(pset.get<double>(        "minimumTime",500)), // nsec
      maxT_(pset.get<double>(        "maximumTime",2000)), // nsec
      trigMode_(pset.get<bool>(      "trigMode",true)),
      printLevel_(pset.get<int>(     "printLevel",0)),
      diagLevel_(pset.get<int>(      "diagLevel",0)),
      end_{TrkTypes::cal,TrkTypes::hv},
      strawDigis_(pset.get<std::string>("StrawDigis","makeSD")),
      caloClusterModuleLabel_(pset.get<std::string>("caloClusterModuleLabel","CaloClusterFast")),
      peakfit_(pset.get<fhicl::ParameterSet>("PeakFitter",fhicl::ParameterSet()))
  {
      produces<StrawHitCollection>();
      produces<StrawHitFlagCollection>();
      produces<StrawHitPositionCollection>();
      
      if (printLevel_ > 0) std::cout << "In StrawHitsFromStrawDigis2 constructor " << std::endl;
  }

  StrawHitsFromStrawDigis2::~StrawHitsFromStrawDigis2() {}

  
  //------------------------------------------------------------------------------------------
  void StrawHitsFromStrawDigis2::beginJob()
  {
  }

  void StrawHitsFromStrawDigis2::beginRun(art::Run& run)
  {    
      ConditionsHandle<StrawElectronics> strawele = ConditionsHandle<StrawElectronics>("ignored");
                 
      // this must be done here because strawele is not accessible at startup and it 
      // contains a const refenence to pfit_, so this can't be instanciated earliere
      if (fittype_ == TrkChargeReco::FitType::peakminusped)
         pfit_ = std::unique_ptr<TrkChargeReco::PeakFit>(new TrkChargeReco::PeakFit(*strawele,peakfit_) );
      else if (fittype_ == TrkChargeReco::FitType::combopeakfit)
	 pfit_ = std::unique_ptr<TrkChargeReco::PeakFit>(new TrkChargeReco::ComboPeakFitRoot(*strawele,peakfit_) );
      else
	 pfit_ = std::unique_ptr<TrkChargeReco::PeakFit>(new TrkChargeReco::PeakFitRoot(*strawele,peakfit_) );
                      
      if (printLevel_ > 0) std::cout << "In StrawHitsFromStrawDigis2 begin Run " << std::endl;
  }



  //------------------------------------------------------------------------------------------
  void StrawHitsFromStrawDigis2::produce(art::Event& event)
  {        
      if (printLevel_ > 0) std::cout << "In StrawHitsFromStrawDigis2 produce " << std::endl;

      const Tracker& tracker = getTrackerOrThrow();
      const TTracker& tt(*GeomHandle<TTracker>());
      size_t nplanes = tt.nPlanes();
      size_t npanels = tt.getPlane(0).nPanels();
      
      ConditionsHandle<StrawElectronics> strawele = ConditionsHandle<StrawElectronics>("ignored");
      ConditionsHandle<StrawPhysics> strawphys = ConditionsHandle<StrawPhysics>("ignored");
      ConditionsHandle<TrackerCalibrations> tcal("ignored");

      
      art::Handle<StrawDigiCollection> strawdigisHandle;
      event.getByLabel(strawDigis_,strawdigisHandle);
      const StrawDigiCollection& strawdigis(*strawdigisHandle);
      
      const CaloClusterCollection* caloClusters(0);
      art::Handle<CaloClusterCollection> caloClusterHandle;
      if (event.getByLabel(caloClusterModuleLabel_, caloClusterHandle)) caloClusters = caloClusterHandle.product();

      std::unique_ptr<StrawHitCollection> strawHits(new StrawHitCollection);
      strawHits->reserve(strawdigis.size());
      std::unique_ptr<StrawHitFlagCollection> strawHitFlags(new StrawHitFlagCollection);
      strawHitFlags->reserve(strawdigis.size());
      std::unique_ptr<StrawHitPositionCollection> strawHitPositions(new StrawHitPositionCollection);
      strawHitPositions->reserve(strawdigis.size());      


      std::vector<std::vector<size_t> > hits_by_panel(nplanes*npanels,std::vector<size_t>());    
      std::vector<size_t> largeHits, largeHitPanels;
      largeHits.reserve(strawdigis.size());
      largeHitPanels.reserve(strawdigis.size());
      
      SHInfo shinfo;



      for (size_t isd=0;isd<strawdigis.size();++isd)
      {
          const StrawDigi& digi = strawdigis[isd];          
          TDCTimes times;
          strawele->tdcTimes(digi.TDC(),times);
	  TOTTimes tots{0.0,0.0};
	  for(size_t iend=0;iend<2;++iend){
	    tots[iend] = digi.TOT(end_[iend])*strawele->totLSB();
	  }
	  // take the earliest of the 2 end times
	  float time = std::min(times[0],times[1]);

          //calorimeter filtering
          if (usecc_ && caloClusters)
          {
             bool outsideCaloTime(true);
             for (const auto& cluster : *caloClusters) 
               if (std::abs(time-cluster.time())<clusterDt_) {outsideCaloTime=false; break;}
             if (outsideCaloTime) continue;
          }
	  
          
          //prefiltering on time if needed
          if (trigMode_ && (time < minT_ ||time > maxT_)) continue;


          //extract energy from waveform
	  // note: pedestal is being subtracting inside strawele, in the real experiment we will need
	  // per-channel version of this FIXME!!!
	  TrkChargeReco::PeakFitParams params;
	  pfit_->process(digi.adcWaveform(),params);
 	  double energy = strawphys->ionizationEnergy(params._charge/strawphys->strawGain());
	  if (printLevel_ > 1) std::cout << "Fit status = " << params._status << " NDF = " << params._ndf << " chisquared " << params._chi2
	                                 << " Fit charge = " << params._charge << " Fit time = " << params._time << std::endl;
                                                  
          //create straw hit, mc info if requested
          const Straw& straw  = tracker.getStraw( digi.strawIndex() );
          StrawHit hit(digi.strawIndex(),times,tots,energy);
          
          StrawHitFlag flag;
          if (energy > minE_ && energy < maxE_) flag.merge(StrawHitFlag::energysel);
          if (time > minT_ && time < maxT_)     flag.merge(StrawHitFlag::timesel);
          if (usecc_)                           flag.merge(StrawHitFlag::calosel);
          
          StrawHitPosition shp;
          tcal->StrawHitInfo(straw,hit,shinfo);
          shp._pos   = shinfo._pos;
          shp._phi   = shinfo._pos.phi();
          shp._wdir  = straw.getDirection();
          shp._wdist = shinfo._tddist;
          shp._wres  = shinfo._tdres;
          shp._tres  = straw.getRadius()*0.288675135; //0.28867 = 1/sqrt(12)          
          if (shinfo._tdiv) shp._flag.merge(StrawHitFlag::tdiv); 
                         

          //buffer large hit for cross-talk analysis
          size_t iplane       = straw.id().getPlane();
          size_t ipnl         = straw.id().getPanel();
          size_t global_panel = ipnl + iplane*npanels;

          hits_by_panel[global_panel].push_back(strawHits->size());          
          if (energy >= ctE_) {largeHits.push_back(strawHits->size()); largeHitPanels.push_back(global_panel);}        

          strawHits->push_back(std::move(hit));          
          strawHitFlags->push_back(std::move(flag));          
          strawHitPositions->push_back(std::move(shp));
      }
      
      
      //flag straw and electronic cross-talk
      for (size_t ilarge=0; ilarge < largeHits.size();++ilarge)
      {
          const StrawHit& sh = strawHits->at(largeHits[ilarge]);
          const Straw& straw = tracker.getStraw( sh.strawIndex() );
          for (size_t jsh : hits_by_panel[largeHitPanels[ilarge]])
          {
              if (jsh==largeHits[ilarge]) continue;              
              const StrawHit& sh2 = strawHits->at(jsh);              
              if (sh2.time()-sh.time() > ctMinT_ && sh2.time()-sh.time() < ctMaxT_)
              {
                 if (straw.isSamePreamp(sh2.strawIndex()))       strawHitFlags->at(jsh).merge(StrawHitFlag::elecxtalk);
                 if (straw.isNearestNeighbour(sh2.strawIndex())) strawHitFlags->at(jsh).merge(StrawHitFlag::strawxtalk);
              }           
          }
      }
             

      event.put(std::move(strawHits));
      event.put(std::move(strawHitFlags));
      event.put(std::move(strawHitPositions));


   }

}


using mu2e::StrawHitsFromStrawDigis2;
DEFINE_ART_MODULE(StrawHitsFromStrawDigis2);

