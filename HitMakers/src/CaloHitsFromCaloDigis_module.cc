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

    explicit CaloHitsFromCaloDigis(fhicl::ParameterSet const& pset) :

      // Parameters

      _diagLevel                  (pset.get<int>                ("diagLevel")),
      _debugMode                  (pset.get<int>                ("debugMode")),
      _digiSampling               (pset.get<double>             ("digiSampling")),    // ns
      _acquisitionLenght          (pset.get<double>             ("acquisitionLenght")),    // ns
      _caloCalibNoise             (pset.get<double>             ("caloCalibNoise")),
      _g4ModuleLabel              (pset.get<string>             ("g4ModuleLabel")),
      _toff                       (pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
      _caloDigiModuleLabel        (pset.get<std::string>        ("caloDigiModuleLabel", "CaloDigisFromStepPointMCs")),
      _stepPoints                 (pset.get<string>             ("calorimeterStepPoints","calorimeter")),
      _caloShowerStepMCModuleLabel(pset.get<std::string>        ("caloShowerStepMCModuleLabel")), 
      _caloShowerMCName           (pset.get<std::string>        ("caloShowerMCName")),		  
      _caloROShowerMCName         (pset.get<std::string>        ("caloROShowerMCName"))	  
      //      _DAQTimeThreshold           (pset.get<double>             ("DAQTimeThreshold")) //ns
    {
      //      _wave_point_error = 1.;

      _waveformProcessor = new FitWaveformProcess(pset.get<fhicl::ParameterSet>("FitWaveformProcessor",fhicl::ParameterSet()));

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

    static Double_t logn(Double_t *x, Double_t *par);

    
  private:
    int    _diagLevel;
    int    _debugMode;

    double _digiSampling;
    double _acquisitionLenght;

    double _caloCalibNoise;

    string _g4ModuleLabel;  // Name of the module that made the input hits.
    SimParticleTimeOffset _toff;     // time offset smearing

    string _caloDigiModuleLabel; // Name of the module that made the calo hits.

    std::string _stepPoints;

    std::string            _caloShowerStepMCModuleLabel;   
    std::string            _caloShowerMCName;   
    std::string            _caloROShowerMCName;   

    //    double            _DAQTimeThreshold;
    int               _nHits[5], _nSamples[5];
    double            _wfWithPileUp;
    
    WaveformProcess*  _waveformProcessor;
    

    //some diagnostic histograms
    Hist_              _hist;
    
    //    double            _waveLengthMap[3500];

    const CaloShowerStepMCCollection*        _caloShowerStepCol;
    //----------------------------------------------------------------------//
    void              makeRecoCaloDigis   (CaloDigiCollection    &    CaloDigiColl,
					   CaloDigiMCCollection  &    CaloDigiMCColl,
					   RecoCaloDigiCollection&    RecoCaloHits, 
					   const CaloDigiPacked*      CaloDigis   );   

    void              processWaveform(int &Index, 
				      const std::vector<int>   *CaloFromDigi, 
				      CaloDigiCollection    &   CaloDigiColl,
				      CaloDigiMCCollection  &   CaloDigiMCColl,
				      RecoCaloDigiCollection& RecoCaloHits);
  
    //    void              fillDiag(unique_ptr<CaloDigiCollection> CaloDigis, unique_ptr<RecoCaloDigiCollection> RecoHits);  

    //    double            _eHitMax;

    const Calorimeter*   _calorimeter; // cached pointer to the calorimeter geometry

  };



  //formula from: Particle Detectors, C. Grupen , B. Shwartz
  Double_t CaloHitsFromCaloDigis::logn(Double_t *x, Double_t *par) {
    Double_t Epeak, sigma, eta, norm;
    Double_t Aterm;
    Double_t logterms0,s0;
    Double_t logn,logterm;
    Double_t expterm;
    Double_t pigreco=3.14159265;
      
    eta = par[0];
    sigma = par[1];
    Epeak = par[2];
    norm = par[3];

    logterms0 = eta*2.35/2+sqrt(1+pow((eta*2.35/2),2));
    s0 = (2/2.35)*log(logterms0);
    
    Aterm = eta/(sqrt(2*pigreco)*sigma*s0);

    logterm = 1-(eta/sigma)*(x[0]-Epeak);  

      
    if(logterm<0){
      logterm = 0.0001;
    }
    expterm = log(logterm)/s0;
    expterm = -0.5*pow(expterm,2);

    logn = norm*Aterm *exp(expterm);      
    return logn;
  }  // fine function

  void CaloHitsFromCaloDigis::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _waveformProcessor->book();

    _hist._hEMC       = tfs->make<TH1F>("hEMC%i",
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

  void  CaloHitsFromCaloDigis::endJob(){}

//-----------------------------------------------------------------------------

  void CaloHitsFromCaloDigis::beginRun(art::Run &){
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
  }

  void CaloHitsFromCaloDigis::produce(art::Event& event) {
    if ( _debugMode > 0 ) {
      printf("[CaloHitsFromCaloDigis::produce] event %i\n", event.id().event());
    }

    //get the time offset
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
    //   fillDiag(caloDigiColl, recoHitColl);
      
    // }

    if ( _debugMode > 0 ) {
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
    std::vector<int>   tmpObj  = CaloDigis->output();
    std::vector<int>* caloFromDigi = &tmpObj;

    //total lenght of the digitized waveforms
    int                cfdDim         = caloFromDigi->size();

    //index used for looping in the array.
    //It starts form one because the element 0 represents the size of the array
    int                stringIndex(1);

    if (_debugMode > 0){
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
    //    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    //get the RO id
    int        roId      = CaloFromDigi->at(Index);
    
    //get the lenght of the digitized hit
    int        hitLenght = CaloFromDigi->at(Index+1);
    
    //instantiate some help parameters
    int        wfSize, wfOffset;
    int        content;
    
    std::vector<int>       waveform;


    const CaloShowerStepMC* h;
    double      tEnd;         
    int         nShowerSteps; 
    double      hitTime(0);  
    int         crystalId;    
    double      adc2MeV;     

    //the digitization include a time offset between the MC time
    // and the waveform time. FIX ME
    double      timeOffset = 30;

    //check the lenght: if it is 2 that means that it is empty
    if ( hitLenght != 2) {//goto  PROCESS_END;
    
      //lenght of the waveform
      wfSize   = CaloFromDigi->at(Index + 3) - 2;
    
      wfOffset = Index + 4;

      //fill the histogram
      for (int i=0; i<wfSize; ++i){
	content = CaloFromDigi->at(wfOffset+i);

	waveform.push_back(content);
      }

      //time of the first digitized timestamp
      int         t0     = CaloFromDigi->at(Index + 2);
   
      //store the unpacked Calo-Hit
      CaloDigiColl.push_back(CaloDigi(waveform, roId, hitLenght, t0));
    
      //fill MC info
      tEnd          = t0 + waveform.size()*_digiSampling;
      nShowerSteps  = _caloShowerStepCol->size();
      crystalId     = _calorimeter->crystalByRO(roId);


      CaloDigiMC caloDigiMC;
      caloDigiMC.init();

      int        found=0;
      for (int i=0; i<nShowerSteps; ++i){
	h =  &_caloShowerStepCol->at(i);
	if (h->volumeId() != crystalId)                continue;

	hitTime   = h->time() + _toff.totalTimeOffset(h->simParticle())  + timeOffset;
      
	//FIX ME, there is an offset between data and MC
	if ( (hitTime > t0) && (hitTime < tEnd)){
	  caloDigiMC.addCaloShower(h);
	  found = 1;
	}
      }

      if (found == 0){
	printf("[] no CaloShowerMC found: t0 = %i hitTime = %3.2f \n", 
	       t0, hitTime);
      }
      //store CaloDigMC
      CaloDigiMCColl.push_back(caloDigiMC);
    
      adc2MeV    = calorimeterCalibrations->ADC2MeV(roId);
    
      _waveformProcessor->processWaveform(adc2MeV, CaloDigiColl.at(CaloDigiColl.size() - 1), RecoCaloHits);


      if (_diagLevel > 0) {
	_waveformProcessor->fillDiag(&caloDigiMC, &RecoCaloHits);
      }
    
    }//  PROCESS_END:;
    Index += hitLenght;
  }

}

using mu2e::CaloHitsFromCaloDigis;
DEFINE_ART_MODULE(CaloHitsFromCaloDigis);
