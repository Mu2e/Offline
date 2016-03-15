//
// An EDProducer Module that reads CaloHit objects and turns them into
// RecoCaloCrystalHit objects, collection
//
// $Id:  $
// $Author:  $
// $Date:  $
//
// Original author G. Pezzullo & L. Morescalchi
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
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"


#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloDigiPacked.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/RecoCaloDigiCollection.hh"
#include "MCDataProducts/inc/CaloDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloDigiMC.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepMC.hh"
#include "MCDataProducts/inc/CaloDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloDigiMC.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/PDGCode.hh"

#include "HitMakers/inc/WaveformProcess.hh"
#include "HitMakers/inc/FitWaveformProcess.hh"
#include "HitMakers/inc/FitPolWaveformProcess.hh"
#include "HitMakers/inc/RawWaveformProcess.hh"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGaussQ.h"


using namespace std;
using art::Event;


namespace mu2e {

  struct Hist_ {
    //some histograms for the occupancy and
    //estimation of the total data flux
    TH1F*  _hNSamples[5];
    TH1F*  _hNHits   [5];   
    TH1F*  _hWFLength[5];
    TH1F*  _hEMC;

  };

  class CaloHitsFromCaloDigis : public art::EDProducer {

  public:
    enum processorStrategy {
      LogNormalFit = 0,
      PolFit       = 1,

      RawProcessor = 2
    };

    explicit CaloHitsFromCaloDigis(fhicl::ParameterSet const& pset) :

      // Parameters

      _diagLevel                  (pset.get<int>                ("diagLevel")),
      _debugLevel                 (pset.get<int>                ("debugLevel")),
      _digiSampling               (pset.get<double>             ("digiSampling")),       // ns
      _caloCalibNoise             (pset.get<double>             ("caloCalibNoise")),
      _toff                       (pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
      _mbbuffer                   (pset.get<double>             ("TimeFoldingBuffer")),  // ns
      _blindTime                  (pset.get<double>             ("blindTime" )),         // ns
      _caloDigiModuleLabel        (pset.get<std::string>        ("caloDigiModuleLabel", "CaloDigisFromStepPointMCs")),
      _stepPoints                 (pset.get<string>             ("calorimeterStepPoints","calorimeter")),
      _caloShowerStepMCModuleLabel(pset.get<std::string>        ("caloShowerStepMCModuleLabel")), 
      _caloShowerMCName           (pset.get<std::string>        ("caloShowerMCName")),		  
      _caloROShowerMCName         (pset.get<std::string>        ("caloROShowerMCName")),
      _fillCaloDigiMC             (pset.get<int>                ("fillCaloDigiMC")),
      _processorStrategy          (pset.get<int>                ("processorStrategy"))
    {

      switch (_processorStrategy){
      case LogNormalFit: default:
	_waveformProcessor = new FitWaveformProcess   (pset.get<fhicl::ParameterSet>("FitWaveformProcessor",
										     fhicl::ParameterSet()));
	break;
      case PolFit:
	_waveformProcessor = new FitPolWaveformProcess(pset.get<fhicl::ParameterSet>("FitPolWaveformProcessor",
										     fhicl::ParameterSet()));
	break;
      case RawProcessor:
	_waveformProcessor = new RawWaveformProcess   (pset.get<fhicl::ParameterSet>("RawProcessor",
										     fhicl::ParameterSet()));
	break;
      }


      if (_diagLevel > 0) {
	_fillCaloDigiMC = 1;
      }
      
      // Tell the framework what we make.
      produces<CaloDigiCollection>();
      produces<CaloDigiMCCollection>();
      produces<RecoCaloDigiCollection>();
    }
    virtual ~CaloHitsFromCaloDigis() { }


    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(art::Run &);
    void produce( art::Event& e);

  private:
    int                        _diagLevel;
    int                        _debugLevel;

    double                     _digiSampling;

    double                     _caloCalibNoise;

    SimParticleTimeOffset      _toff;     // time offset smearing
    double                     _mbbuffer;
    double                     _mbtime;
    double                     _blindTime;

    string                     _caloDigiModuleLabel; // Name of the module that made the calo hits.
			       
    std::string                _stepPoints;

    std::string                 _caloShowerStepMCModuleLabel;   
    std::string                 _caloShowerMCName;   
    std::string                 _caloROShowerMCName;   

    int                         _fillCaloDigiMC;
    int                         _processorStrategy;

    int                         _nHits[5], _nSamples[5];
    double                      _wfWithPileUp;
    		                
    WaveformProcess*            _waveformProcessor;
    

    //some diagnostic histograms
    Hist_                       _hist;
    

    const CaloShowerStepMCCollection*        _caloShowerStepCol;
    
    void              fillCaloDigiMC      (double T0, double TEnd, int CrystalId,
					   CaloDigiMC            &CaloDigiMC,
					   CaloDigiMCCollection  &CaloDigiMCColl,
					   int                   &Found);


    void              makeRecoCaloDigis   (CaloDigiCollection    &    CaloDigiColl,
					   CaloDigiMCCollection  &    CaloDigiMCColl,
					   RecoCaloDigiCollection&    RecoCaloHits, 
					   const CaloDigiPacked*      CaloDigis   );   

    void              processWaveform(int &Index, 
				      const std::vector<int>   *CaloFromDigi, 
				      CaloDigiCollection    &   CaloDigiColl,
				      CaloDigiMCCollection  &   CaloDigiMCColl,
				      RecoCaloDigiCollection& RecoCaloHits);
  
    void              unfoldHitTime  (const CaloShowerStepMC*  h,  double & HitTimme);

    const Calorimeter*   _calorimeter; // cached pointer to the calorimeter geometry
    
  };

  void CaloHitsFromCaloDigis::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _waveformProcessor->book();

    _hist._hEMC       = tfs->make<TH1F>("hEMC",
					"CaloDigMC energy; E [MeV]",
					400, 0, 200);
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
    
    
  }
    


  //------------------------------------------------------------------------------------------
  void   CaloHitsFromCaloDigis::fillCaloDigiMC      (double T0, double TEnd, int CrystalId,
						     CaloDigiMC            &CaloDigiMC,
						     CaloDigiMCCollection  &CaloDigiMCColl,
						     int    &Found){
    int         nShowerSteps  = _caloShowerStepCol->size();
    double      hitTime(0);

    const CaloShowerStepMC* h;

  

    hitTime          = -1;

    //the digitization include a time offset between the MC time
    // and the waveform time. FIX ME
    double      timeOffset = 50;
    T0          -= timeOffset;
      
    for (int i=0; i<nShowerSteps; ++i){
      h =  &_caloShowerStepCol->at(i);
	
      // if (h->volumeId() == 249 ){
      //   h->print(std::cout);
      // }
      if (h->volumeId() != CrystalId)                continue;

      unfoldHitTime(h, hitTime);

      //hits created before _blindTime have not been used, just skip them
      if (hitTime < _blindTime - 100)                      continue;

      if ( (hitTime > T0) && (hitTime < TEnd)){
	CaloDigiMC.addCaloShower(h, hitTime);
	Found = 1;
      }else{
	if ( _debugLevel > 0){
	  printf("[CaloHitsFromCaloDigis::processWaveform] hitTime = %.3f T0 = %5.3f TEnd = %5.3f\n",
		 hitTime, T0, TEnd); 
	}
      }
    }

  
    //store CaloDigMC
    CaloDigiMCColl.push_back(CaloDigiMC);
    
  
  }


  //--------------------------------------------------------------------------------
  
  void CaloHitsFromCaloDigis::unfoldHitTime  (const CaloShowerStepMC* Hit,  double & HitTime){
    double hitTimeUnfolded = _toff.totalTimeOffset(Hit->simParticle()) + Hit->time();
      
    HitTime                = fmod(hitTimeUnfolded, _mbtime);

      if (_debugLevel>0){
	double      totalTimeOffset = _toff.totalTimeOffset(Hit->simParticle()) ;
	Hit->print(cout);
	printf("[CaloHitsFromCaloDigis::unfoldHitTime] totalTimeOffSet = %4.2f\n", totalTimeOffset);
      }
      
      if (HitTime < _mbbuffer) {
	if (HitTime+_mbtime > _blindTime) {
	  HitTime = HitTime + _mbtime;
	}
      }
      else {
	if (HitTime > (_mbtime - _mbbuffer)) {
	  if (HitTime - _mbtime > _blindTime) {
	    HitTime =   HitTime - _mbtime;
	  }
	}
      }
      
      if (_debugLevel>0){
	printf("[CaloHitsFromCaloDigis::unfoldHitTime] HitTime = %4.2f\n", HitTime);
      }

    }

  
  void  CaloHitsFromCaloDigis::endJob(){}

//-----------------------------------------------------------------------------

  void CaloHitsFromCaloDigis::beginRun(art::Run &){
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
  }

  void CaloHitsFromCaloDigis::produce(art::Event& event) {
    if ( _debugLevel > 2 ) {
      printf("[CaloHitsFromCaloDigis::produce] event %i\n", event.id().event());
    }

    //get the time offset
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;
    _toff.updateMap(event);

 
    if (_diagLevel > 0){
      // for (int i=0; i<3500; ++i){
      // 	_waveLengthMap[i] = 0.;
      // }

      //reset some helping variables
      for (int i=0; i<5;++i){
	_nHits   [i] = 0;
	_nSamples[i] = 0;
      }
    }
    
    // Check that calorimeter geometry description exists
    art::ServiceHandle<GeometryService> geom;    
    if( !(geom->hasElement<Calorimeter>()) )                return;


    //Get handles to calorimeter RO (aka APD) collection
    art::Handle<CaloDigiPacked> caloDigisHandle;
    event.getByLabel(_caloDigiModuleLabel, caloDigisHandle);
    if ( !caloDigisHandle.isValid()) return;
    const CaloDigiPacked* caloDigis = caloDigisHandle.product();   //FIX
     
    art::Handle<CaloShowerStepMCCollection> caloShowerStepMCHandle,caloROShowerStepMCHandle;

    //event.getByLabel(_caloShowerStepMCModuleLabel, _caloROShowerMCName, caloROShowerStepMCHandle);  
    if (    event.getByLabel(_caloShowerStepMCModuleLabel, _caloShowerMCName,   caloShowerStepMCHandle)){
      _caloShowerStepCol = caloShowerStepMCHandle.product();
    }else {
      _caloShowerStepCol = 0;
    }

  
    //Create a new RecoCaloCrystalHit collection and fill it
    unique_ptr<CaloDigiCollection>     caloDigiColl  (new CaloDigiCollection    );
    unique_ptr<CaloDigiMCCollection>   caloDigiMCColl(new CaloDigiMCCollection  );
    unique_ptr<RecoCaloDigiCollection> recoHitColl   (new RecoCaloDigiCollection);
  
    makeRecoCaloDigis(*caloDigiColl, *caloDigiMCColl, *recoHitColl, caloDigis);

    if (_diagLevel > 0){
      for (int i=0; i<int(caloDigiMCColl->size()); ++i){
	_hist._hEMC ->Fill(caloDigiMCColl->at(i).totalEDep());
      }
    }


    if ( _debugLevel > 2 ) {
      printf("[CaloHitsFromCaloDigis::produce] produced RecoCrystalHits ");
      printf(": caloDigiColl size  = %i", int(caloDigiColl->size()));
      printf(", recoHitColl  size  = %i \n", int(recoHitColl->size()));
    }

    //    event.put(std::move(caloHits));
    event.put(std::move(caloDigiColl));
    event.put(std::move(caloDigiMCColl));
    event.put(std::move(recoHitColl));

    return;
  }

    
  void CaloHitsFromCaloDigis::makeRecoCaloDigis(CaloDigiCollection    &    CaloDigiColl,
						CaloDigiMCCollection  &    CaloDigiMCColl,
						RecoCaloDigiCollection&    RecoCaloHits, 
						const CaloDigiPacked*      CaloDigis    ){
    
    //get the array which holds all the digitized waveforms
    std::vector<int>        tmpObj       = CaloDigis->output();
    std::vector<int>*       caloFromDigi = &tmpObj;

    //total lenght of the digitized waveforms
    int                cfdDim         = caloFromDigi->size();

    //index used for looping in the array.
    //It starts form one because the element 0 represents the size of the array
    int                stringIndex(1);

    if (_debugLevel > 2){
      printf("[CaloHitsFromCaloDigis::processWaveform]      charge    MC-time     tStamp-0    tStamp-1     tStamp-2    tStamp-3    tStamp-4    tStamp-5\n");
    }
    while ( stringIndex < cfdDim ){
      processWaveform(stringIndex, caloFromDigi, CaloDigiColl, CaloDigiMCColl, RecoCaloHits);
    }

  }
  
  void CaloHitsFromCaloDigis::processWaveform(int&                      Index, 
					      const std::vector<int>*   CaloFromDigi,
					      CaloDigiCollection    &   CaloDigiColl,
					      CaloDigiMCCollection  &   CaloDigiMCColl,
					      RecoCaloDigiCollection&   RecoCaloHits){

    ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

    //get the RO id
    int        roId      = CaloFromDigi->at(Index);
    
    //get the lenght of the digitized hit
    int        hitLenght = CaloFromDigi->at(Index+1);
    
    //instantiate some help parameters
    int        wfSize, wfOffset;
    int        content;
    
    std::vector<int>       waveform;

    double     adc2MeV;     

    //check the lenght: if it is 2 that means that it is empty
    if ( hitLenght != 2) {//goto  PROCESS_END;
    
      //lenght of the waveform
      wfSize   = CaloFromDigi->at(Index + 3) - 2;
    
      wfOffset = Index + 4;

      
      for (int i=0; i<wfSize; ++i)
	{
	  content = CaloFromDigi->at(wfOffset+i);
	  waveform.push_back(content);
	}

      //time of the first digitized timestamp
      int         t0     = CaloFromDigi->at(Index + 2);

      //index of the new  CaloDigi int he output collection 
      int         vecIndex  = CaloDigiColl.size();

      //store the unpacked Calo-Hit
      CaloDigiColl.push_back(CaloDigi(waveform, roId, wfSize, t0, vecIndex));
    
      adc2MeV    = calorimeterCalibrations->ADC2MeV(roId);
    
      _waveformProcessor->processWaveform(adc2MeV, CaloDigiColl.at(CaloDigiColl.size() - 1), RecoCaloHits);


       //fill MC info
      if (_fillCaloDigiMC == 1){
	double      tEnd          = t0 + waveform.size()*_digiSampling;
	int         crystalId     = _calorimeter->crystalByRO(roId);
	int         found(0);
	CaloDigiMC  caloDigiMC;

	caloDigiMC.init();
   	fillCaloDigiMC(double(t0), tEnd, crystalId, caloDigiMC, CaloDigiMCColl, found);

	if ( (found == 0) && ( _debugLevel > 0) ){
	  printf("[CaloHitsFromCaloDigis::processWaveform] no CaloShowerMC found: T0 = %i hitLenght = %d CrystalId = %d\n", 
		 t0, hitLenght, crystalId);
	  CaloDigiColl.at(CaloDigiColl.size()-1).print();
	}

	if (_diagLevel > 0) {
	  _waveformProcessor->fillDiag(&caloDigiMC, &RecoCaloHits);
	}
      }
      
   
    
    }//  PROCESS_END:;
    Index += hitLenght;
  }

}

using mu2e::CaloHitsFromCaloDigis;
DEFINE_ART_MODULE(CaloHitsFromCaloDigis);
