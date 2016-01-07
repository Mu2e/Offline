//
// An EDProducer Module that reads CaloHit objects and turns them into
// RecoCaloCrystalHit objects, collection
//
// $Id: CaloHitsFromCaloDigis_module.cc,v 1.17 2014/08/01 20:57:45 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:45 $
//
// Original author KLG
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// ROOT includes
#include "TH1F.h"
#include "TF1.h"
#include "TSpline.h"
#include "TFile.h"
#include "CalPatRec/inc/THackData.hh"
#include "TROOT.h"
#include "TFolder.h"
#include "TTree.h"
#include "TH2F.h"

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/sort_functors.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"


#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/RecoCaloDigi.hh"
#include "RecoDataProducts/inc/RecoCaloDigiCollection.hh"
#include "MCDataProducts/inc/CaloDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloDigiMC.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/PDGCode.hh"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGaussQ.h"


using namespace std;
using art::Event;


namespace mu2e {

  struct Hist_t {
    TH1F*  _hEdep;
    TH1F*  _hAmplitude;
    TH1F*  _hEdepMerged;
    TH1F*  _hTime;

    TH2F*  _hTimeMostE;
    TH2F*  _hTimeVsE;
    TH2F*  _hTimeVsAmp;
    TH2F*  _hAmpVsE;
    TH2F*  _hPileUpVsAmp;
    
    TH1F*  _hTimeMerged;
    TH1F*  _hChi2PileUp;
    TH1F*  _hChi2Time;
    TH1F*  _hNDofTime;
    TH1F*  _hFitEta;
    TH1F*  _hFitPeak;
    TH1F*  _hFitSigma;
    TH1F*  _hFitNorm;
    TH2F*  _hEdepAmp;

 //some histograms for the occupancy and
    //estimation of the total data flux
    TH1F*  _hNSamples[5];
    TH1F*  _hNHits   [5];   
    TH1F*  _hWFLength[5];

    TH2F*  _hWFLengthVsAmp;
    
    TH2F*  _DAQRateMap;
  };

  class CaloCrystalHitsFromCaloHits : public art::EDProducer {

  public:
    typedef art::Ptr<RecoCaloDigi>           RecoCaloDigiPtr;

    explicit CaloCrystalHitsFromCaloHits(fhicl::ParameterSet const& pset) :

      // Parameters

      _diagLevel           (pset.get<int>                ("diagLevel")),		      
      _debugMode           (pset.get<int>                ("debugMode")),		      
      _time4Merge          (pset.get<double>             ("time4Merge")), // MeV	      
      _caloDigisModuleLabel(pset.get<std::string>        ("caloDigisModuleLabel"))
    {
  
      // Tell the framework what we make.
      produces<CaloCrystalHitCollection>();
    }
    virtual ~CaloCrystalHitsFromCaloHits() { }


    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(art::Run &);
    void produce( art::Event& e);

  private:

    int    _diagLevel;
    int    _debugMode;

    double _time4Merge; 
    
    string _caloDigisModuleLabel; // Name of the module that made the calo hits.

    int               _nHits[5], _nSamples[5];
    int               _wfWithPileUp;
    std::vector<int> * _hitMap[2000]; 
    
    //some diagnostic histograms
    Hist_t              _hist;
    TH1F*               _hDt;
    
    TTree*              _tree;
    
    Int_t               _evt, _run, _nROHits, _hitCounter;
    Float_t             _time[100000], _Chi2Time[100000], _nDof[100000], _fitEta[100000], _fitNorm[100000], _fitSigma[100000], _fitPeak[100000]
      , _amp[100000], _charge[100000];

    

    //----------------------------------------------------------------------//
    void              makeCaloHits   (CaloCrystalHitCollection& caloHits, 
				      const art::Handle<RecoCaloDigiCollection>   RecoCaloDigisHandle);
    //    const RecoCaloDigiCollection  recoCaloHits) ;
    int               pairRO         (int ROid);

    double            _eHitMax, _timeReco;

    const Calorimeter*   _calorimeter; // cached pointer to the calorimeter geometry

  };

  void CaloCrystalHitsFromCaloHits::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
    _hist._hEdep         = tfs->make<TH1F>("hEdep","Hit energy deposition",200,0.,500);
    _hist._hAmplitude    = tfs->make<TH1F>("hAmplitude","Waveform amplitude", 2500,0.,2500);
    _hist._hEdepMerged   = tfs->make<TH1F>("hEdepMerged","Crystal energy deposition",400, 0., 200);
    _hist._hTime         = tfs->make<TH1F>("hTime","Hit time ", 12000, 0., 2000);
    _hist._hTimeMerged   = tfs->make<TH1F>("hTimeMerged","Crystal time ", 12000, 0., 2000);
    _hist._hTimeMostE    = tfs->make<TH2F>("hTimeMostE",
					   "Hit time of the most energetic hit vs reco charge",
					   500 ,   0., 500,
					   1000, -20.,  20);
    _hist._hTimeVsE      = tfs->make<TH2F>("hTimeVsE",
					   "Hit time resolution vs reco charge",
					   500 ,   0., 500,
					   1000, -20.,  20);
    _hist._hTimeVsAmp    = tfs->make<TH2F>("hTimeVsAmplitude",
					      "Hit amplitude resolution vs reco charge",
					   1250 ,   0., 2500,
					   1000, -20.,  20);

    _hist._hAmpVsE       = tfs->make<TH2F>("hAmpVsE","Amplitude versus charge; Amplitue [mV]; Energy [MeV]",
					   2500,0.,2500,
					   500 ,   0., 500);
    _hist._hPileUpVsAmp  = tfs->make<TH2F>("hPileUpVsAmp","Pile-up disriminator versus amplitude",
					   500 ,   0.,  500,
					   2500,   0., 2500);
    

    _hist._hChi2PileUp   = tfs->make<TH1F>("hChi2PileUp","#Chi^2 from pileup finder  ", 1000, 0., 10);
    _hist._hChi2Time     = tfs->make<TH1F>("hChi2Time","#Chi^2 from time reco alg:#Chi2^2/ndof"   , 4000, 0., 200);
    _hist._hNDofTime     = tfs->make<TH1F>("hNDofTime","nDof from time reconstruction fit; nDof [#]", 200, 0., 10);
    
    //histograms for the time reconstruction fit results
    _hist._hFitEta       = tfs->make<TH1F>("hFitEta", "Fit - #eta distribution; #eta", 100, -5, 5);
    _hist._hFitPeak      = tfs->make<TH1F>("hFitPeak","Fit - peak distribution; peak [ns]",4000, 0.,  2000);
    _hist._hFitSigma     = tfs->make<TH1F>("hFitSigma","Fit - #sigma distribution; #sigma [ns]", 100, 0, 100);
    _hist._hFitNorm      = tfs->make<TH1F>("hFitNorm","Fit - Norm distribution; Norm [a.u.]", 5000, 0, 10000);



    _hDt                 = tfs->make<TH1F>("hDt","#Deltat distribution; #Delta t = t_{reco} t_{MC} [ns]", 1600, -200., 200);
 

    for (int i=0; i<5; ++i){
      _hist._hNSamples[i]   = tfs->make<TH1F>(Form("hNSamples%i",i),
					      "Numerb of samples / event distribution; nSamples/event [#] ", 
					      10000, 1e2, 1e5); 
      _hist._hNHits   [i]   = tfs->make<TH1F>(Form("hNHits%i",i),
					      "Numerb of hits / event distribution; nHits/event [#] ",
					      10000, 0., 5000);
      _hist._hWFLength[i]   = tfs->make<TH1F>(Form("hWFLength%i",i),
					      "wavefrom length distribution; waveform length [# samples]",
					      1000, 0, 1000);

    }
    

    _hist._hWFLengthVsAmp = tfs->make<TH2F>("hWFLengthVsAmp",
					    "wavefrom length vs amplitude distribution; waveform length [# samples]; Amp [ADC counts]",
					    1000, 0, 1000, 
					    2000, 0, 2000);


    _hist._DAQRateMap = tfs->make<TH2F>("DAQRateMap",
					"rate per channel; channel is [#]; Expected rate [# samples/event]",
					3500, 0, 3500, 
					200, 0, 200);

    //------------------------------------------------------------------------------------------//
    
    //clean the vector used for unpackaging the waveforms
    for (int i=0; i<2000; ++i){
      _hitMap[i] = new std::vector<int>();
    }
    
    if (_diagLevel > 1){
      //create the  TTree      
       _tree  = tfs->make<TTree>("Calo", "Calo");

       _tree->Branch("evt",          &_evt ,        "evt/I");
       _tree->Branch("run",          &_run ,        "run/I");
       _tree->Branch("nROHits",      &_nROHits,     "nROHits/I");

       _tree->Branch("time",         &_time ,        "time[nROHits]/F");
       _tree->Branch("Chi2Time",     &_Chi2Time ,    "Chi2Time[nROHits]/F");
       _tree->Branch("nDof",         &_nDof ,        "nDof[nROHits]/F");
       _tree->Branch("fitEta",       &_fitEta ,      "fitEta[nROHits]/F");
       _tree->Branch("fitNorm",      &_fitNorm ,     "fitNorm[nROHits]/F");
       _tree->Branch("fitSigma",     &_fitSigma ,    "fitSigma[nROHits]/F");
       _tree->Branch("fitPeak",      &_fitPeak ,     "fitPeak[nROHits]/F");

       _tree->Branch("amp",          &_amp ,         "amp[nROHits]/F");
       _tree->Branch("charge",       &_charge ,      "charge[nROHits]/F");
    }
    
    _eHitMax = -1.;
  }

  void  CaloCrystalHitsFromCaloHits::endJob(){
    
    for (int i=0; i<2000; ++i){
      delete _hitMap[i];
    }
    
  }


//-----------------------------------------------------------------------------
  void CaloCrystalHitsFromCaloHits::beginRun(art::Run &){
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
  }

//-----------------------------------------------------------------------------

  void CaloCrystalHitsFromCaloHits::produce(art::Event& event) {
    if ( _debugMode > 0 ) {
      printf("[CaloCrystalHitsFromCaloHits::produce] begins \n");
    }

    //clear the hitMap
    for (int i=0; i<2000; ++i){
      _hitMap[i]->clear();
    }

    //Get handles to calorimeter RO (aka APD) collection
    art::Handle<RecoCaloDigiCollection>       recoCaloDigisHandle;
    event.getByLabel(_caloDigisModuleLabel, recoCaloDigisHandle);
    if ( !recoCaloDigisHandle.isValid()) return;
    const RecoCaloDigiCollection* recoCaloDigis = recoCaloDigisHandle.product();   
     
    _evt        = event.id().event();
    _run        = event.run();
    _nROHits    = 0;
    _hitCounter = 0;

    //Create a new RecoCaloCrystalHit collection and fill it
    unique_ptr<CaloCrystalHitCollection> caloHits(new CaloCrystalHitCollection);
  
    int    nRecoDigis = recoCaloDigis->size();
    int    roId;

    // fill the map that associate for each crystal all the  corresponding RecoCaloDigi indexes
    for (int i=0; i<nRecoDigis; ++i){
      roId = recoCaloDigis->at(i).ROid();
      _hitMap[_calorimeter->crystalByRO(roId)]->push_back(i);
    }

    //reate the CaloCrystalHitCollection
    //    makeCaloHits(*caloHits, *recoCaloDigis);
    makeCaloHits(*caloHits, recoCaloDigisHandle);

    if (_diagLevel > 0){
      _nROHits =  _hitCounter;

      if (_diagLevel > 1){
	_tree->Fill();
      }
    }

    if ( _debugMode > 0 ) {
      printf("[CaloCrystalHitsFromCaloHits::produce] produced RecoCrystalHits ");
      printf(": caloHits.size()  = %i \n", int(caloHits->size()));
    }

    event.put(std::move(caloHits));

    return;
  }

  void CaloCrystalHitsFromCaloHits::makeCaloHits(CaloCrystalHitCollection   & CaloHits, 
						 //					   const RecoCaloDigiCollection      RecoCaloDigis) {
						 const art::Handle<RecoCaloDigiCollection>   RecoCaloDigisHandle) {
         
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());
      
    const RecoCaloDigiCollection* RecoCaloDigis = RecoCaloDigisHandle.product();   

    //define few helping variables
    
    for(unsigned int icry = 0; icry<cal.nCrystal();icry++){
	
      int nhits(0);
      if (_hitMap[icry] != 0) {
	nhits = _hitMap[icry]->size();
      }

      if (nhits == 0)             continue;// goto CRYSTAL_LOOP_ENS;// continue; 

      int pairingVec[nhits] = {0};    
	
                                              //loop over the RO
      for(int ipv1 = 0; ipv1<nhits; ipv1++){	    
	if(pairingVec[ipv1] == -1)              continue;
      	    
	int paired = 0;
	
	int indexWF1 = _hitMap[icry]->at(ipv1);
	int indexWF2(-1);

	for(int ipv2 = ipv1+1; ipv2<nhits; ipv2++){
	  
	  if (_diagLevel > 0){
	    _hist._hTime->Fill(RecoCaloDigis->at(ipv2).time());
	    _hist._hEdep->Fill(RecoCaloDigis->at(ipv2).edep());
	  }

	  if(pairingVec[ipv2] == -1) continue;
	  
	  //set the index at which is located the second wf
	  indexWF2 = _hitMap[icry]->at(ipv2);

	  //check if the hits come from the same RO
	  if(RecoCaloDigis->at(indexWF1).ROid() == RecoCaloDigis->at(indexWF2).ROid())        continue;
              
	  double tdiff = fabs(RecoCaloDigis->at(indexWF1).time() - RecoCaloDigis->at(indexWF2).time());

	  //check if the time difference of the pulses is compatible
	  //is time enough or we should add a check over the amplitudes? 
	  if(tdiff < _time4Merge){
		
	    //RecoCaloDigiCollection    crystalHitRecoDigis;
	    //	    crystalHitRecoDigis.push_back(RecoCaloDigis->at(indexWF1));
	    //	    crystalHitRecoDigis.push_back(RecoCaloDigis->at(indexWF2));
	    std::vector<RecoCaloDigiPtr>    crystalHitRecoDigis;
	    crystalHitRecoDigis.push_back(art::Ptr<RecoCaloDigi>(RecoCaloDigisHandle, indexWF1));
	    crystalHitRecoDigis.push_back(art::Ptr<RecoCaloDigi>(RecoCaloDigisHandle, indexWF2));


	    CaloCrystalHit caloCrystalHit(icry, crystalHitRecoDigis);
	    CaloHits.push_back(caloCrystalHit);
		
	    if (_diagLevel > 0){
	      _hist._hTimeMerged->Fill(caloCrystalHit.time());
	      _hist._hEdepMerged->Fill(caloCrystalHit.energyDep());
	    }
	    
	    //set these RO-pulses as used
	    pairingVec[ipv1] = -1;
	    pairingVec[ipv2] = -1;
	    paired           = 1;
		
	  }
	      
	}
	    
	if(paired == 0) {//FIXME
	     
	  std::vector<RecoCaloDigiPtr>    crystalHitRecoDigis;
	  //	  RecoCaloDigiCollection    crystalHitRecoDigis;
	  crystalHitRecoDigis.push_back(art::Ptr<RecoCaloDigi>(RecoCaloDigisHandle, indexWF1));
	  //	  crystalHitRecoDigis.push_back(RecoCaloDigis->at(indexWF1));

	  CaloCrystalHit caloCrystalHit(icry, crystalHitRecoDigis);

	  CaloHits.push_back(caloCrystalHit);

	  //set the RO-pulse as used	      
	  pairingVec[ipv1] = -1;
	}
	    
      }  //    end of the loop over the RO-pulses
    } //end loop over the crystals
	    	
  } 
  
  
  int  CaloCrystalHitsFromCaloHits::pairRO(int ROid){
   
    if(ROid/2 == 0) return ROid+1;
    else return ROid-1;
    
  }

  
}

using mu2e::CaloCrystalHitsFromCaloHits;
DEFINE_ART_MODULE(CaloCrystalHitsFromCaloHits);
