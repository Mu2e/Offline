//
// An EDProducer Module that reads Caloshower and produces CaloDigiPacked
//
// Original author G. Pezzullo
//


// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>

// ROOT includes
#include "TSpline.h"
//#include "CalPatRec/inc/THackData.hh"

#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"
#include "TArc.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TFile.h"
#include "TVector2.h"
#include "TMarker.h"
#include "TCanvas.h"

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
//#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/sort_functors.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "HitMakers/inc/CaloCrystalMCUtil.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloDigiPacked.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepMC.hh"
#include "MCDataProducts/inc/CaloShowerCollection.hh"
#include "MCDataProducts/inc/CaloShower.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoisson.h"



class THackData;

using namespace std;

namespace mu2e {

  struct Hist_t {
    TH1F*  _hEdep;
    TH1F*  _hTime;
    
    //some histograms for the occupancy and
    //estimation of the total data flux
    TH1F*  _hNSamples;
    TH1F*  _hNHits;   
    TH1F*  _hNSamplesVsR;
    TH1F*  _hNHitsVsR;   
    TH1F*  _hWFLength;

    TH2F*  _hWFLengthVsAmp;
  };

  //--------------------------------------------------------------------
  //
  //
  class CaloDigisFromCaloShower : public art::EDProducer {
  
  public:
   
    // First vector is list of crystal steps, associated with particular readout element.
    // Second vector is list of readout steps, associated with particular readout element.

    explicit CaloDigisFromCaloShower(fhicl::ParameterSet const& pset) :

      // Parameters
      _diagLevel                  (pset.get<int>        ("diagLevel" )),		  
      _debugLevel                 (pset.get<int>        ("debugLevel")),		  
      _wfInput                    (pset.get<int>        ("wfInput"   )),	   	  
      _blindTime                  (pset.get<double>     ("blindTime" )),         // ns
      _crystal_refractive_index   (pset.get<double>     ("crystal_refractive_index")),         // 
      _engine                     ( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
      _randGauss                  ( _engine ),
      _randPoisson                ( _engine ),
      _caloShowerModuleLabel      (pset.get<std::string>("caloShowerModuleLabel")), 
      _g4ModuleLabel              (pset.get<string>     ("g4ModuleLabel")),	  
      _toff                       (pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
      _mbbuffer                   (pset.get<double>     ("TimeFoldingBuffer")),  // ns
      _addNoise                   (pset.get<int>       ("addNoise")),           //flag for adding or not Gaussian noise
      _addLightPropagation         (pset.get<int>       ("addLightPropagation")),  
      _noise                      (pset.get<double>     ("noise"   )),           // mV 
      _thresholdVoltage           (pset.get<double>     ("thresholdVoltage"  )), // mV 
      _thresholdAmplitude         (pset.get<double>     ("thresholdAmplitude")), //mV
      _DAQTimeThreshold           (pset.get<double>     ("DAQTimeThreshold"  )), //ns  
      _energyScale                (pset.get<double>     ("energyScale" )),       // mV/MeV
      _digiSampling               (pset.get<double>     ("digiSampling")),       // ns
      _nBits                      (pset.get<int>        ("nBits"       )),       // number of digitizer bits
      _dynamicRange               (pset.get<double>     ("dynamicRange")),       // mV
      _acquisitionEndtime          (pset.get<double>    ("acquisitionEndTime")),  // ns
      _bufferAfter                (pset.get<int>        ("bufferAfter"      )),  //# timestamps	  
      _bufferBefore               (pset.get<int>        ("bufferBefore"     )),	 //# timestamps 
      _waveformOffset             (pset.get<double>     ("waveformOffset"   )),   // ns
      _pulseIntegralSteps         (pset.get<int>        ("pulseIntegralSteps"))	 //# integral steps
    {  
      

      // Tell the framework what we make.
      produces<CaloDigiPacked>();    

      _integralPrecision    = _digiSampling/_pulseIntegralSteps;
 
    }
    
    virtual ~CaloDigisFromCaloShower() { }
    virtual void beginRun(art::Run& );
    virtual void beginJob();
    virtual void endJob();

    void produce( art::Event& e);

  private:

 
    // typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
    // typedef art::Ptr<StepPointMC>   StepPtr;
    // typedef std::vector<StepPtr >   StepPtrs;
    // typedef art::Ptr<SimParticle>   SimPtr;
    // typedef std::vector<SimPtr >    SimPtrs;
    // typedef std::map<int,StepPtrs > HitMap;

	 
    int                     _diagLevel;  
    int                     _debugLevel;  
    int                     _wfInput;

    double                  _blindTime;
    double                  _crystal_refractive_index;

    CLHEP::HepRandomEngine& _engine;
    CLHEP::RandGaussQ       _randGauss;
    CLHEP::RandPoisson      _randPoisson;


    std::string            _caloShowerModuleLabel;   
    std::string            _g4ModuleLabel;   // Name of the module that made these hits.

    SimParticleTimeOffset  _toff;     // time offset smearing
    double                 _mbtime;   // period of 1 microbunch
    double                 _mbbuffer; // buffer on that for ghost hits (wrapping)

    int                    _addNoise;
    int                    _addLightPropagation;
    double                 _noise;
    double                 _thresholdVoltage;         
    int                    _thresholdAmplitude;
    double                 _DAQTimeThreshold;

    double                 _energyScale;
    double                 _digiSampling;
    int                    _nBits;
    double                 _dynamicRange;
    double                 _acquisitionEndtime;
    int                    _bufferAfter;
    int                    _bufferBefore;
    double                 _waveformOffset;    //offset needed for getting the first non null point of the input waveform
    int                    _pulseIntegralSteps;
    double                 _integralPrecision;

    //some diagnostic histograms
    Hist_t                _hist;
    THackData*            _hackData;
    TFolder*              _folder;
 
    int                   _nHits, _nSamples;
      

    void   addNoiseToWaveforms    ();
   
    void   createOutput           (CaloDigiPacked& CaloDigis);

    void   makeDigitization       (art::Handle<CaloShowerCollection>const& CaloShowersHandle,
				   CaloDigiPacked&);


 
    void   fillWaveforms             (art::Handle<CaloShowerCollection>const& crystalStepsHandles);

    void   digitizeWaveform          (int ROId, double Time, double Edep);
    
    void   includeLightPropagation   (double& HitTime, CLHEP::Hep3Vector& PosInCrystalFrame);

    void   readoutResponse           (double Edep, double Time, int CrId);
    
    int                  _nWaveforms;

    TH1F*                _pshape;

    std::vector<int> *   _waveforms[4000];

    std::vector<double>  _pulseDigitized[100];

    int                  _ROFilled[4000];

    const Calorimeter*   _calorimeter; // cached pointer to the calorimeter geometry

    //    double convolutionFunction(double time, double tau, double sigma);
  };

  void   CaloDigisFromCaloShower::includeLightPropagation(double& HitTime, 
							    CLHEP::Hep3Vector& PosInCrystalFrame){
    
    double      path          = PosInCrystalFrame.z();
    double      vLight        = CLHEP::c_light/_crystal_refractive_index;

    HitTime    =  HitTime + path/vLight;
  }

  
  //We want the output in the following data format:
  // nTotWords - nWords_roID - roiID - nWord_roID_ihit - time_roID_ihit - Data_roID_ihit - ... 
  void CaloDigisFromCaloShower::createOutput(CaloDigiPacked& CaloDigis){
    
    int nTotWords = 1;
    std::vector<int> caloDigiOutput; 
      
    double   ADCTomV = _dynamicRange/pow(2.,_nBits);

    for (int it=0; it<_nWaveforms; ++it){

      if (_ROFilled [it] == 0 )        continue;

      std::vector<int>* itWave = _waveforms[it];
      std::vector<int>  output; 
      std::vector< std::vector<int> > MCoutput;

      int         nRoWords = 2;
     
      int         length = 0;
      int         nPeaks = 0;


      //some helping variables
      int         nHitWords(-1), t0(-1), waveformLength, crId;
      double      roRadius, crX, crY;

      // Loop on time stamps of each waveform
      
      int  nSamples = itWave->size();
      for (int timeSample=0; timeSample < nSamples; ++timeSample){
	int         waveContent = itWave->at(timeSample);
	double      funcValue   = waveContent*ADCTomV;

	if (_debugLevel > 1){
	  if (waveContent > 0){
	    printf("wfContent (%4i, %9.3f) \n", waveContent, funcValue);
	  }
	}
	
	if (funcValue > _thresholdVoltage) {
	  
	  int     sampleStart = timeSample;
	  int     sampleStop  = timeSample;
	  int     sampleMax   = -1;
	  // The storing of the peak starts #bufferBefore entries before the first time stamp over the threshold
	  
	  if (sampleStart > _bufferBefore){
	    sampleStart -=  _bufferBefore;
	  }else {
	    sampleStart = 0;
	  }

	  while( itWave->at(sampleStop)*ADCTomV > _thresholdVoltage && 
		 sampleStop < nSamples-1) {
	   
	    waveContent = itWave->at(sampleStop);
	    
	    sampleStop++;

	    if (waveContent > sampleMax) {
	      sampleMax = waveContent;
	    }

	  }

	  if (nSamples - sampleStop > _bufferAfter){
	    sampleStop += _bufferAfter;
	  }else {
	    sampleStop  = nSamples;
	  }

	  
	  //look the amplitude of the waveform and decide to store it or not
	  if (sampleMax*ADCTomV < _thresholdAmplitude)        goto	  NEXT_WAVEFORM;

	  for (int i=sampleStart; i<sampleStop; ++i){
	    output.push_back(itWave->at(i));
	  }

	  // Waveform starting time
       	  nHitWords = 2 + sampleStop - sampleStart;

	  //2015-12-15 G. Pezzullo added the following change
	  if ( nHitWords > 2){
	    t0  = sampleStart*_digiSampling;
       
	    output.insert(output.begin()+length+nPeaks*2, nHitWords);
	    output.insert(output.begin()+length+nPeaks*2, t0);
	    
	    length += sampleStop - sampleStart;
	    nRoWords += nHitWords;
	    nPeaks++;
	  }

	  if (_diagLevel > 0){
	    if (t0 > _DAQTimeThreshold){
	      ++_nHits;

	      waveformLength = sampleStop - sampleStart;	    

	      _nSamples += waveformLength;

	      crId      = _calorimeter->crystalByRO(it);
	      crX       = _calorimeter->crystal(crId).localPosition().x();
	      crY       = _calorimeter->crystal(crId).localPosition().y();
	      roRadius  = sqrt(crX*crX + crY*crY);

	      _hist._hNSamplesVsR   ->Fill(roRadius, waveformLength);
	      _hist._hNHitsVsR      ->Fill(roRadius);
	      _hist._hWFLength      ->Fill(waveformLength);
	      _hist._hWFLengthVsAmp ->Fill(waveformLength, sampleMax);
	    }
	  }
	  
	NEXT_WAVEFORM:;
	  /////increment timeSample
	  timeSample = sampleStop;


	}

      }
    
	

      int iRO = it;

      output.insert(output.begin(), nRoWords);
      output.insert(output.begin(), iRO);
	
      nTotWords += nRoWords;
      
      caloDigiOutput.insert(caloDigiOutput.end(), output.begin(), output.end());
    }
    
    caloDigiOutput.insert(caloDigiOutput.begin(), nTotWords);
    CaloDigis = CaloDigiPacked(caloDigiOutput); 


    if (_diagLevel > 0){
      _hist._hNSamples   ->Fill(_nSamples);
      _hist._hNHits      ->Fill(_nHits);
    }
	  


    if (_debugLevel > 20){
      int     size = caloDigiOutput.size();
      int     content;
      printf("[CaloDigisFromCaloShower::createOutput] caloDigiOutput\n");
      for (int i=0; i<size; ++i){
	content =  caloDigiOutput.at(i);
	printf("[CaloDigisFromCaloShower::createOutput] %i \n", content);
      }
    }
  }


  //2015-09-14 Gianipez and L. Morscalchi:
  //add noise to the digitized waveforms
  void CaloDigisFromCaloShower::addNoiseToWaveforms(){

    for (int it=0; it<_nWaveforms; ++it){
      int size = (_mbtime - _blindTime - _acquisitionEndtime)/_digiSampling;

      if (_ROFilled [it] == 0 )        continue;

      for (int sample=0; sample<size; ++sample){
	//inlcude Gaussian noise if requested
	double content = _waveforms[it]->at(sample);
	content +=  _randGauss.fire(0, _noise)*pow(2,_nBits)/_dynamicRange; // _noise is in mV, sample is in counts
	if (content < 0) {
	  content = 0;
	}
	_waveforms[it]->at(sample) = content;
      }
    }
  }

  
  //2015-09-14 Gianipez and L. Morscalchi:
  //the following function is supposed to propagate the steppoint energy deposition to the photosensor
  // accoridng to the optical transportat. Ray tracing is not implemented yet, so the energy deposition 
  // is evenly splitted btween the two photosensors
  void   CaloDigisFromCaloShower::readoutResponse(double Edep, double Time, int CrId){
    
    int ROidBase = _calorimeter->ROBaseByCrystal(CrId);
    int nROs     = _calorimeter->caloGeomInfo().nROPerCrystal();

    for (int i=0; i<nROs; ++i){
      int ROId = ROidBase + i;
      digitizeWaveform(ROId, Time, Edep);
    }
    
  }

  //--------------------------------------------------------------------
  void   CaloDigisFromCaloShower::digitizeWaveform( int ROId, double Time, double Edep){

    int     size         = _waveforms[ROId]->size();

    //keep memory of the stored RO
    _ROFilled [ROId]     = 1;

    //calculate the amplitude of the pulse
    double  pulseAmp     = Edep*_energyScale;
    double  maxADCCounts = pow(2., _nBits);
    
    if (_debugLevel > 10){
      printf("[CaloDigisFromCaloShower::makeCalorimeterHits]   eDep = %9.3f MeV amplitude = %9.3f mV time = %9.3f ns\n", 
	     Edep, pulseAmp, Time);
      printf("[CaloDigisFromCaloShower::makeCalorimeterHits]   timeSample   |   ADC-counts   |    wfAmp    |   signalAmp  \n");
    }

    int          startSample    = Time/_digiSampling;
    int          precisionIndex = (Time/_digiSampling - int(Time/_digiSampling)) / _integralPrecision;

    for (int timeSample = startSample; timeSample < size; ++timeSample){

      // Energy deposition is converted to counts per bin 1 count = 1 V x 2^nBits / dynamicRange 

      double  funcValue = _pulseDigitized[precisionIndex].at(timeSample - startSample);

      //now convert it into ADC counts
      double  wfAmp     = pulseAmp*funcValue;
      int     ADCCounts = wfAmp / _dynamicRange * maxADCCounts;
      _waveforms[ROId]->at(timeSample)      += ADCCounts;//norm*funcValue;

      //check if the ADC saturates
      if ( _waveforms[ROId]->at(timeSample) > maxADCCounts){
	_waveforms[ROId]->at(timeSample) = maxADCCounts;
      }
      
      if (_debugLevel > 10){
	if(wfAmp > 0){
	  printf("[CaloDigisFromCaloShower::makeCalorimeterHits]   %4i   |   %4i   |    %9.3f  | %9.3f \n", 
		 timeSample, ADCCounts,  wfAmp, funcValue);
	}
      }
      
    }
  }
  
  void CaloDigisFromCaloShower::beginRun(art::Run& ){
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
    _nWaveforms  = _calorimeter->nCrystal()*_calorimeter->caloGeomInfo().nROPerCrystal();

  }

  void CaloDigisFromCaloShower::beginJob(){
  
    TFile *f;

//     TFile *f = new TFile("/mu2e/data/users/gianipez/test-waveform-2015.root", "R");
//     f->cd();

    if (_wfInput == 0)
    {
      f = TFile::Open("/mu2e/data/users/gianipez/test-waveforms-2015.root", "R");
      _pshape = (TH1F*) gDirectory->Get("pshape");
      f->Close();
    }
    if (_wfInput == 1)
    {
      f = TFile::Open("/mu2e/data/users/gianipez/CsI-waveform-2016-02-03.root", "R");
      _pshape = (TH1F*) gDirectory->Get("h1");
      f->Close();
    } 
    if (_wfInput == 2)
    {
      f = TFile::Open("/mu2e/data/users/gianipez/CsI-waveform-2016-02-03.root", "R");
      _pshape = (TH1F*) gDirectory->Get("h2"); 
      f->Close();
    }
    if (_wfInput == 3)
    {
      f = TFile::Open("/mu2e/data/users/gianipez/test-CsI-2015-10-06.root");
      _pshape = (TH1F*) gDirectory->Get("CsIPulse");
      f->Close();
    }
    if (_wfInput == 4)
    {
      f = TFile::Open("/mu2e/data/users/gianipez/Mu2eCaloSgnlTempl.root");
      _pshape = (TH1F*) gDirectory->Get("hS");
      f->Close();
    }
    
    double  pulseBinWidth  = _pshape->GetBinWidth  (1);
    double  pulseStartTime = _pshape->GetBinLowEdge(1);
    double  pulseRange     = 1695.;//ns

    int     nTimeStamps    = pulseRange/_digiSampling;
    //    int     startTimeStamp = _waveformOffset/_digiSampling;//pulseBinWidth;
    int     binMin(0), binMax(0);
    int     nBinTimeStamp  = _digiSampling/pulseBinWidth;

    double  pulseNorm      = _pshape->Integral();
    double  pulseIntegral(0);


    
    for (int i=0; i<_pulseIntegralSteps; ++i)
      {
	binMin   = (pulseStartTime + _waveformOffset+ i*_integralPrecision)/pulseBinWidth;
	
	for (int j=0; j<nTimeStamps; ++j)
	  {
	    binMax          = binMin + nBinTimeStamp;
	    pulseIntegral   = _pshape->Integral(binMin, binMax)/pulseNorm;
	    binMin          = binMin + nBinTimeStamp;

	    _pulseDigitized[i].push_back(pulseIntegral);
	  }
      }
    


    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("diag");

    _hist._hEdep          = tfdir.make<TH1F>("hEdep","Hit energy deposition",200,0.,100);
    _hist._hTime          = tfdir.make<TH1F>("hTime","Hit time ", 4000, 0., 2000);

    _hist._hNSamples      = tfdir.make<TH1F>("hNSamples","Numerb of samples / event distribution; nSamples/event [#] ", 
					     100000, 5e3, 5e5);
    _hist._hNHits         = tfdir.make<TH1F>("hNHits","Numerb of hits / event distribution; nHits/event [#] ", 10000, 0., 10000);
    
    _hist._hNSamplesVsR   = tfdir.make<TH1F>("hNSamplesVsR",
					    "Numerb of samples / event distribution vs R [mm]; radius [mm]; nSamples [#]",
					    50, 300, 800);
    _hist._hNHitsVsR      = tfdir.make<TH1F>("hNHitsVsR",
					    "Numerb of hits / event distribution vs R [mm]; radius [mm]; nHits [#]",
					    50, 300, 800);
   
    _hist._hWFLength      = tfdir.make<TH1F>("hWFLength","wavefrom length distribution; waveform length [# ssmples]",
					    1000, 0, 1000);

    _hist._hWFLengthVsAmp = tfdir.make<TH2F>("hWFLengthVsAmp",
					    "wavefrom length vs amplitude distribution; waveform length [# ssmples]; Amp [ADC counts]",
					    1000, 0, 1000, 
					    2000, 0, 2000);
     

    for (int i=0; i<4000; ++i){
      _waveforms[i] = new std::vector<int>();
    }
    
  }

  void CaloDigisFromCaloShower::endJob(){

    for (int i=0; i<4000; ++i)
      {
	delete _waveforms[i];
      }

  }

  void  CaloDigisFromCaloShower::produce(art::Event& event) {

    if ( _diagLevel > 0 ) 
      {
	_nHits    = 0;
	_nSamples = 0;
      }
    
    if (_debugLevel > 0 ) 
      { 
	cout << "CaloDigisFromCaloShower: produce() begin" << endl;
      }

    static int ncalls(0);
    ++ncalls;

    // Check that calorimeter geometry description exists
    art::ServiceHandle<GeometryService> geom;    
    if( !(geom->hasElement<Calorimeter>()) )         return;
   

    //update condition cache
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;
    _toff.updateMap(event);
   
    unique_ptr<CaloDigiPacked>          caloDigis(new CaloDigiPacked);

    art::Handle<CaloShowerCollection> caloShowerHandle;
    event.getByLabel(_caloShowerModuleLabel,   caloShowerHandle);

     
    //clear the holder of the waveforms
    for (int i=0; i<_nWaveforms; ++i){
      _waveforms[i]->clear();
      _ROFilled [i] = 0;
    }
  
    //initizlize the waveforms
    for (int it=0; it<_nWaveforms; ++it){
      
      int size = (_mbtime - _blindTime - _acquisitionEndtime)/_digiSampling;

      for (int sample=0; sample<size; ++sample){
	_waveforms[it]->push_back(0);
      }
    }
    
    makeDigitization(caloShowerHandle, *caloDigis);
    
    // Add the output hit collection to the event
    event.put(std::move(caloDigis));

    if ( _debugLevel > 0 ) cout << "CaloDigisFromCaloShower: produce() end" << endl;

  } 

  void CaloDigisFromCaloShower::makeDigitization (art::Handle<CaloShowerCollection>const& CaloShowerCollHandle,
						  CaloDigiPacked& CaloDigis){
    
    
     // Fill the waveforms
    fillWaveforms( CaloShowerCollHandle);
    
    if (_addNoise == 1){                                        //add noise to the waveforms
      addNoiseToWaveforms();
    }
    
    createOutput(CaloDigis);

  } // end makeCalorimeterHits

  //-----------------------------------------------------------------------------
  void CaloDigisFromCaloShower::fillWaveforms(art::Handle<CaloShowerCollection>const& CaloShowerCollHandle) {
        
    GlobalConstantsHandle<ParticleDataTable> pdt;
  
    // Handle to the conditions service
    ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
  
    CaloShowerCollection const& crystalStepsHandles(*CaloShowerCollHandle);
  
    int    NCompressedHits = int(crystalStepsHandles.size());

    for (int i=0; i < NCompressedHits; ++i){

      CaloShower const& h =  crystalStepsHandles.at(i);
      if ( h.energy() <= 0.0 )         continue;
	double      edep_corr = h.energy();
	int         crystalId = h.crystalId();

	CLHEP::Hep3Vector  posInCrystalFrame  =_calorimeter->crystal(crystalId).localPositionFF();

	double      hitTime         =  h.time();//fmod(hitTimeUnfolded,_mbtime);

	if (_addLightPropagation == 1)
	  {
	    includeLightPropagation(hitTime, posInCrystalFrame);
	  }
	
	if (_debugLevel>0)
	  {
	    printf("[CaloDigisFromCaloShower::fillWaveforms] HitTime = %4.2f\n", hitTime);
	  }
	
	//avoid the use of hits created before _blindTime
	if (hitTime < _blindTime)               continue;
	
	readoutResponse(edep_corr, hitTime, crystalId);

	if (_diagLevel > 0){
	  _hist._hEdep  ->Fill(edep_corr);
	  _hist._hTime  ->Fill(hitTime);
	}
	
      }
    }


} // end namespace mu2e

using mu2e::CaloDigisFromCaloShower;
DEFINE_ART_MODULE(CaloDigisFromCaloShower);



