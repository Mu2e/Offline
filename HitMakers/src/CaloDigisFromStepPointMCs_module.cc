//
// An EDProducer Module that reads StepPointMC objects and turns them into
// CaloHit, CaloHitMCTruth, CaloHitMCPtr, CrystalOnlyHit, CrystalHitMCPtr
// objects.
//
// Original author L. Morescalchi, G. Pezzullo and S. R. Soleti
//
// Notes
// 1) The CrystalOnlyHitCollection is a form of MC truth on a per crystal basis.
//    It represents an idea per crystal response if no readouts were hit directly.
//
// 2) We still have two times the crystal StepPointMC hits saved per hit in a crystal. Note sure we can 
//    keep only a single copy before simulating the whole readout chain, so keep them for now. 


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
#include "messagefacility/MessageLogger/MessageLogger.h"

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
  class CaloDigisFromStepPointMCs : public art::EDProducer {
  private:

 
    typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
    typedef art::Ptr<StepPointMC>   StepPtr;
    typedef std::vector<StepPtr >   StepPtrs;
    typedef art::Ptr<SimParticle>   SimPtr;
    typedef std::vector<SimPtr >    SimPtrs;
    typedef std::map<int,StepPtrs > HitMap;

	 
    int   _diagLevel;  
    int   _debugLevel;  
    int   _wfInput;

    double                _blindTime;
    double                _crystal_refractive_index;
    bool                  _caloLRUcorrection;
    bool                  _caloNonLinCorrection;
    int                   _caloChargeProductionEffects;

    CLHEP::HepRandomEngine& _engine;
    CLHEP::RandGaussQ       _randGauss;
    CLHEP::RandPoisson      _randPoisson;


    std::string _stepPoints;
    std::string _rostepPoints;
    std::string            _caloShowerStepMCModuleLabel;   
    std::string            _caloShowerMCName;   
    std::string            _caloROShowerMCName;   

    std::string _g4ModuleLabel;   // Name of the module that made these hits.

    SimParticleTimeOffset _toff;     // time offset smearing
    double                _mbtime;   // period of 1 microbunch
    double                _mbbuffer; // buffer on that for ghost hits (wrapping)


    const std::string _messageCategory;

    int                   _addNoise;
    int                   _addLightPropagation;
    double                _noise;
    double                _thresholdVoltage;         
    int                   _thresholdAmplitude;
    double                _DAQTimeThreshold;

    double                _energyScale;
    double                _digiSampling;
    int                   _nBits;
    double                _dynamicRange;
    double                _acquisitionEndtime;
    int                   _bufferAfter;
    int                   _bufferBefore;
    double                _waveformOffset;    //offset needed for getting the first non null point of the input waveform

    //some diagnostic histograms
    Hist_t                _hist;
    THackData*            _hackData;
    TFolder*              _folder;
 
    int                   _nHits, _nSamples;
      
  public:
   
    // First vector is list of crystal steps, associated with particular readout element.
    // Second vector is list of readout steps, associated with particular readout element.

    explicit CaloDigisFromStepPointMCs(fhicl::ParameterSet const& pset) :

      // Parameters
      _diagLevel                  (pset.get<int>        ("diagLevel" )),		  
      _debugLevel                 (pset.get<int>        ("debugLevel")),		  
      _wfInput                    (pset.get<int>        ("wfInput"   )),	   	  
      _blindTime                  (pset.get<double>     ("blindTime" )),         // ns
      _crystal_refractive_index   (pset.get<double>     ("crystal_refractive_index")),         // 
      _caloLRUcorrection          (pset.get<bool>       ("caloLRUcorrection")),	  
      _caloNonLinCorrection       (pset.get<bool>       ("caloNonLinCorrection")),	  
      _caloChargeProductionEffects(pset.get<int>        ("caloChargeProductionEffects")),
      _engine                     ( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
      _randGauss                  ( _engine ),
      _randPoisson                ( _engine ),
      _stepPoints                 (pset.get<string>     ("calorimeterStepPoints")),
      _rostepPoints               (pset.get<string>     ("calorimeterROStepPoints")),
      _caloShowerStepMCModuleLabel(pset.get<std::string>("caloShowerStepMCModuleLabel")), 
      _caloShowerMCName           (pset.get<std::string>("caloShowerMCName")),		  
      _caloROShowerMCName         (pset.get<std::string>("caloROShowerMCName")),	  
      _g4ModuleLabel              (pset.get<string>     ("g4ModuleLabel")),	  
      _toff                       (pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
      _mbbuffer                   (pset.get<double>     ("TimeFoldingBuffer")),  // ns
      _messageCategory            ("CaloDigisFromStepPointMCs"),			  
      _addNoise                    (pset.get<int>       ("addNoise")),           //flag for adding or not Gaussian noise
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
      _waveformOffset             (pset.get<double>     ("waveformOffset"   ))   // ns
    {  
      

      // Tell the framework what we make.
      produces<CaloDigiPacked>();     
    }
    
    virtual ~CaloDigisFromStepPointMCs() { }
    virtual void beginRun(art::Run& );
    virtual void beginJob();
    virtual void endJob();

    void produce( art::Event& e);

  private:
    
    void addNoiseToWaveforms    ();
    
    void addNuclearCounterEffect();

    void createOutput     (CaloDigiPacked& CaloDigis);

    void makeDigitization (art::Handle<CaloShowerStepMCCollection>const& crystalStepsHandles,
			   art::Handle<CaloShowerStepMCCollection>const& readoutStepsHandles,
			   CaloDigiPacked&);


    void nonLinearityCorrection(CaloShowerStepMC const& h, 
				double& energy, 
				int cryId, 
				ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations,
				GlobalConstantsHandle<ParticleDataTable>& pdt );

    void longitudinalResponseUniformityCorrection(double posZ, 
						  double cryhalflength, 
						  double& energy, 
						  int crid,
						  ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);

    void chargeProductionCorrection(double &Edep, int roId, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibration);

    void fillWaveforms        (art::Handle<CaloShowerStepMCCollection>const& crystalStepsHandles);
    void printStepPointDebug(Calorimeter const & cal,StepPointMC const& h, int crid);


    // Print information about the data products found by the selector functions. 
    void printDataProductInfo( HandleVector const& crystalStepsHandles,
			       HandleVector const& readoutStepsHandles );

    void   digitizeWaveform(int ROId, double Time, double Edep);
    
    void   includeLightPropagation(double& HitTime, CLHEP::Hep3Vector& PosInCrystalFrame);

    void   readoutResponse(double Edep, double Time, int CrId, CLHEP::Hep3Vector PosInCrystalFrame);
    
    int                             _nWaveforms;

    double               _wfnorm;
    TH1F*                _pshape;

    std::vector<int> *   _waveforms[4000];

    int                  _ROFilled[4000];

    const Calorimeter*   _calorimeter; // cached pointer to the calorimeter geometry

    //    double convolutionFunction(double time, double tau, double sigma);
  };

  void   CaloDigisFromStepPointMCs::includeLightPropagation(double& HitTime, 
							    CLHEP::Hep3Vector& PosInCrystalFrame){
    
    //    double      crystalLenght = 2.*_calorimeter->caloGeomInfo().crystalHalfLength();
    
    double      path          = PosInCrystalFrame.z();
    double      vLight        = CLHEP::c_light/_crystal_refractive_index;

    HitTime    =  HitTime + path/vLight;
  }

  
  void CaloDigisFromStepPointMCs::chargeProductionCorrection(double &Edep, int RoId, ConditionsHandle<CalorimeterCalibrations>& CalorimeterCalibrations){
  
    double       lightYield = CalorimeterCalibrations->ROpe(RoId);
    double       mean       = Edep*lightYield/CalorimeterCalibrations->ROfano(RoId);
    double       npe        = _randPoisson.fire(mean);//gRandom->PoissonD(mean);//

    Edep       = npe*CalorimeterCalibrations->ROfano(RoId)/lightYield;
    
  }


  //FIX ME
  //2015-09-14 G. Pezzullo, L. Morescalchi: still need to implement effects
  //on the digitization of the nuclear counter effect
  void CaloDigisFromStepPointMCs::addNuclearCounterEffect(){}
  

  //We want the output in the following data format:
  // nTotWords - nWords_roID - roiID - nWord_roID_ihit - time_roID_ihit - Data_roID_ihit - ... 
  void CaloDigisFromStepPointMCs::createOutput(CaloDigiPacked& CaloDigis){
    
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
      printf("[CaloDigisFromStepPointMCs::createOutput] caloDigiOutput\n");
      for (int i=0; i<size; ++i){
	content =  caloDigiOutput.at(i);
	printf("[CaloDigisFromStepPointMCs::createOutput] %i \n", content);
      }
    }
  }


  //2015-09-14 Gianipez and L. Morscalchi:
  //add noise to the digitized waveforms
  void CaloDigisFromStepPointMCs::addNoiseToWaveforms(){

    for (int it=0; it<_nWaveforms; ++it){
      int size = (_mbtime - _acquisitionEndtime)/_digiSampling;

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
  void   CaloDigisFromStepPointMCs::readoutResponse(double Edep, double Time, int CrId, CLHEP::Hep3Vector PosInCrystalFrame){
    
    int ROidBase = _calorimeter->ROBaseByCrystal(CrId);
    int nROs     = _calorimeter->caloGeomInfo().nROPerCrystal();

    for (int i=0; i<nROs; ++i){
      int ROId = ROidBase + i;
      digitizeWaveform(ROId, Time, Edep);
    }
    
  }

  //--------------------------------------------------------------------
  void   CaloDigisFromStepPointMCs::digitizeWaveform( int ROId, double Time, double Edep){

    int size = _waveforms[ROId]->size();

    //keep memory of the stored RO
    _ROFilled [ROId] =1;

    //calculate the amplitude of the pulse
    double  pulseAmp  = (Edep/_calorimeter->caloGeomInfo().nROPerCrystal())*_energyScale;
    double  binWidth  = _pshape->GetBinWidth(1);

    double  maxADCCounts = pow(2.,_nBits);
    
    if (_debugLevel > 10){
      printf("[CaloDigisFromStepPointMCs::makeCalorimeterHits]   eDep = %9.3f MeV amplitude = %9.3f mV time = %9.3f ns\n", Edep, pulseAmp, Time);
      printf("[CaloDigisFromStepPointMCs::makeCalorimeterHits]   timeSample   |   ADC-counts   |    wfAmp    |   signalAmp  \n");
    }
    
    for (int timeSample = 0; timeSample < size; ++timeSample){

      // Energy deposition is converted to counts per bin 1 count = 1 V x 2^nBits / dynamicRange 
      // TODO: add propagation delay and time smearing effects

      //_waveformOffset is a shift that determines the start of the pulse
      double  timeStamp = timeSample*_digiSampling - Time;
      if (timeStamp < 0.)         continue;
      timeStamp        +=  _waveformOffset;

      int     binmin    = (timeStamp - _digiSampling/2.)/binWidth;
      int     binmax    = (timeStamp + _digiSampling/2.)/binWidth;

      double  funcValue = _pshape->Integral(binmin, binmax) / _digiSampling;

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
	  printf("[CaloDigisFromStepPointMCs::makeCalorimeterHits]   %4i   |   %4i   |    %9.3f  | %9.3f \n", timeSample, ADCCounts,  wfAmp, funcValue);
	}
      }
      
    }
  }
  
  void CaloDigisFromStepPointMCs::beginRun(art::Run& ){
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
    _nWaveforms  = _calorimeter->nCrystal()*_calorimeter->caloGeomInfo().nROPerCrystal();

  }

  void CaloDigisFromStepPointMCs::beginJob(){
  
    TFile *f;

//     TFile *f = new TFile("/mu2e/data/users/gianipez/test-waveform-2015.root", "R");
//     f->cd();

    if (_wfInput == 0){
      f = TFile::Open("/mu2e/data/users/gianipez/test-waveforms-2015.root", "R");
      _pshape = (TH1F*) gDirectory->Get("pshape");
      f->Close();
    }else  if (_wfInput == 1){
      f = TFile::Open("/mu2e/data/users/gianipez/test-waveforms-2015.root", "R");
      _pshape = (TH1F*) gDirectory->Get("histo50");
      f->Close();
    }    if (_wfInput == 2){
      f = TFile::Open("/mu2e/data/users/gianipez/test-waveforms-2015.root", "R");
      _pshape = (TH1F*) gDirectory->Get("histo100"); 
      f->Close();
    }    if (_wfInput == 3){
      f = TFile::Open("/mu2e/data/users/gianipez/test-CsI-2015-10-06.root");
      _pshape = (TH1F*) gDirectory->Get("CsIPulse");
      f->Close();
    }


    _wfnorm = _pshape->Integral();


    art::ServiceHandle<art::TFileService> tfs;
    _hist._hEdep   = tfs->make<TH1F>("hEdep","Hit energy deposition",200,0.,100);
    _hist._hTime   = tfs->make<TH1F>("hTime","Hit time ", 4000, 0., 2000);

    _hist._hNSamples   = tfs->make<TH1F>("hNSamples","Numerb of samples / event distribution; nSamples/event [#] ", 100000, 5e3, 5e5);
    _hist._hNHits      = tfs->make<TH1F>("hNHits","Numerb of hits / event distribution; nHits/event [#] ", 10000, 0., 10000);
    
    _hist._hNSamplesVsR   = tfs->make<TH1F>("hNSamplesVsR",
					    "Numerb of samples / event distribution vs R [mm]; radius [mm]; nSamples [#]",
					    50, 300, 800);
    _hist._hNHitsVsR      = tfs->make<TH1F>("hNHitsVsR",
					    "Numerb of hits / event distribution vs R [mm]; radius [mm]; nHits [#]",
					    50, 300, 800);
   
    _hist._hWFLength      = tfs->make<TH1F>("hWFLength","wavefrom length distribution; waveform length [# ssmples]",
					    1000, 0, 1000);

    _hist._hWFLengthVsAmp = tfs->make<TH2F>("hWFLengthVsAmp",
					    "wavefrom length vs amplitude distribution; waveform length [# ssmples]; Amp [ADC counts]",
					    1000, 0, 1000, 
					    2000, 0, 2000);
     

    for (int i=0; i<4000; ++i){
      _waveforms[i] = new std::vector<int>();
    }
    
  }

  void CaloDigisFromStepPointMCs::endJob(){

    for (int i=0; i<4000; ++i){
      delete _waveforms[i];
    }

  }

  void  CaloDigisFromStepPointMCs::produce(art::Event& event) {

    if ( _diagLevel > 0 ) {
      _nHits    = 0;
      _nSamples = 0;
    }
    
    if (_debugLevel > 0 ) {
      cout << "CaloDigisFromStepPointMCs: produce() begin" << endl;
    }

    static int ncalls(0);
    ++ncalls;

    // Check that calorimeter geometry description exists
    art::ServiceHandle<GeometryService> geom;    
    if( !(geom->hasElement<Calorimeter>()) ) return;
   

    //update condition cache
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;
    _toff.updateMap(event);
   
    unique_ptr<CaloDigiPacked>                         caloDigis   (new CaloDigiPacked);

    // These selectors will select data products with the given instance name, and ignore
    // all other fields of the product ID.
    art::ProductInstanceNameSelector getCrystalSteps(_stepPoints);
    art::ProductInstanceNameSelector getReadoutSteps(_rostepPoints);

    // Get the StepPointMCsNew from the event.
    HandleVector           crystalStepsHandles, readoutStepsHandles;
    event.getMany( getCrystalSteps, crystalStepsHandles);
    event.getMany( getReadoutSteps, readoutStepsHandles);


    art::Handle<CaloShowerStepMCCollection> caloShowerStepMCHandle,caloROShowerStepMCHandle;
    event.getByLabel(_caloShowerStepMCModuleLabel, _caloShowerMCName,   caloShowerStepMCHandle);
    event.getByLabel(_caloShowerStepMCModuleLabel, _caloROShowerMCName, caloROShowerStepMCHandle);
     
    //clear the holder of the waveforms
    for (int i=0; i<_nWaveforms; ++i){
      _waveforms[i]->clear();
      _ROFilled [i] = 0;
    }
  
    //initizlize the waveforms
    for (int it=0; it<_nWaveforms; ++it){
      
      int size = (_mbtime - _acquisitionEndtime)/_digiSampling;

      for (int sample=0; sample<size; ++sample){
	_waveforms[it]->push_back(0);
      }
    }
    
    makeDigitization(caloShowerStepMCHandle, caloROShowerStepMCHandle,
		     *caloDigis);
    
    // Add the output hit collection to the event
    event.put(std::move(caloDigis));

    if ( _debugLevel > 0 ) cout << "CaloDigisFromStepPointMCs: produce() end" << endl;

  } 

  void CaloDigisFromStepPointMCs::makeDigitization (art::Handle<CaloShowerStepMCCollection>const& caloShowerStepMCCollHandle,
						    art::Handle<CaloShowerStepMCCollection>const& caloROShowerStepMCCollHandle,
						    CaloDigiPacked& CaloDigis){
    
    
     // Fill the waveforms
    fillWaveforms( caloShowerStepMCCollHandle);
    
    if (_addNoise == 1){                                        //add noise to the waveforms
      addNoiseToWaveforms();
    }
    
    createOutput(CaloDigis);

  } // end makeCalorimeterHits

  //-----------------------------------------------------------------------------
  void CaloDigisFromStepPointMCs::fillWaveforms(art::Handle<CaloShowerStepMCCollection>const& caloShowerStepMCCollHandle) {
    
    double      first_sp_time = -1;
    double      cryhalflength = _calorimeter->caloGeomInfo().crystalHalfLength();
    
    GlobalConstantsHandle<ParticleDataTable> pdt;
  
    // Handle to the conditions service
    ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
  
    CaloShowerStepMCCollection const& crystalStepsHandles(*caloShowerStepMCCollHandle);
  
    int    NCompressedHits = int(crystalStepsHandles.size());

    for (int i=0; i < NCompressedHits; ++i){

      CaloShowerStepMC const& h =  crystalStepsHandles.at(i);
      if ( h.energy() <= 0.0 )         continue;
	double      edep_corr = h.energy();
	int         crid      = h.volumeId();

	CLHEP::Hep3Vector  posInCrystalFrame  = h.position();

	//should probably go to another class when the code for these corrections grows larger
	if (_caloNonLinCorrection && h.simParticle().isNonnull()) {
	  nonLinearityCorrection(h, edep_corr, crid, calorimeterCalibrations, pdt);
	}

	if (_caloLRUcorrection) {
	  // Calculate correction for edep
	  double posZ = posInCrystalFrame.z();
	  
	  longitudinalResponseUniformityCorrection(posZ, cryhalflength, edep_corr,crid, calorimeterCalibrations);
	}
       
	if (_caloChargeProductionEffects == 1) {
	  chargeProductionCorrection(edep_corr, crid, calorimeterCalibrations);
	}
	
	// time folding and Adding ghost hits to properly treat boundary conditions with folding, see docdb-3425
	double      hitTimeUnfolded = _toff.totalTimeOffset(h.simParticle()) + h.time();
	double      hitTime         = fmod(hitTimeUnfolded,_mbtime);

	if (_debugLevel>0){
	  double      totalTimeOffset = _toff.totalTimeOffset(h.simParticle()) ;
	  h.print(cout);
	  printf("[CaloDigisFromStepPointMCs::fillWaveforms] i = %i totalTimeOffSet = %4.2f\n", i, totalTimeOffset);
	}
	
	if (_addLightPropagation == 1){
	  includeLightPropagation(hitTime, posInCrystalFrame);
	}

	if (hitTime < _mbbuffer) {
	  if (hitTime+_mbtime > _blindTime) {
	    hitTime = hitTime + _mbtime;
	  }
	}
	else {
	  if (hitTime > (_mbtime-_mbbuffer)) {
	    if (hitTime-_mbtime > _blindTime) {
	      hitTime =   hitTime - _mbtime;
	    }
	  }
	}

	if (_debugLevel>0){
	  printf("[CaloDigisFromStepPointMCs::fillWaveforms] HitTime = %4.2f\n", hitTime);
	}
	
	//avoid the use of hits created before _blindTime
	if (hitTime < _blindTime)               continue;
	
	readoutResponse(edep_corr, hitTime, crid, posInCrystalFrame);

	if (_diagLevel > 0){
	  _hist._hEdep  ->Fill(edep_corr);
	  _hist._hTime  ->Fill(hitTime);
	}
	
        if(first_sp_time < 0.) first_sp_time = hitTime;
      }
    }
    
    //    _hackData->fData[20] = first_sp_time;

    
  //  } 


  //-----------------------------------------------------------------------------
  void CaloDigisFromStepPointMCs::printDataProductInfo( HandleVector const& crystalStepsHandles,
							HandleVector const& readoutStepsHandles )
  {
    mf::LogInfo log(_messageCategory);
    log << "MakeCaloReadoutHit::produce will use StepPointMCsNew from: \n";
    for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end();
	  i != e; ++i ){
      art::Provenance const& prov(*(i->provenance()));
      log  << "   " << prov.branchName() << "\n";
    }
    log << "\nMakeCaloReadoutHit::produce will use StepPointMCsNew from: \n";
    for ( HandleVector::const_iterator i=readoutStepsHandles.begin(), e=readoutStepsHandles.end();
	  i != e; ++i ){
      art::Provenance const& prov(*(i->provenance()));
      log  << "   " << prov.branchName() << "\n";
    }
  } 
  
  
  //-----------------------------------------------------------------------------
  void CaloDigisFromStepPointMCs::nonLinearityCorrection(CaloShowerStepMC const& h, 
							 double& energy, int cryId, 
							 ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations, 
							 GlobalConstantsHandle<ParticleDataTable>& pdt  )
  { 
    double edep_save(energy);
    double MeV2keV = 1000.0;
      
    int particleCode = std::abs(h.simParticle()->pdgId());
    if (particleCode!= 11 && particleCode != 22 ) return;
    
    const HepPDT::ParticleData& data = pdt->particle(particleCode).ref();
    double mass = data.mass().value();
    double kinetic_energy = h.energy() - mass;//std::sqrt(h.momentum().mag2() + mass*mass) - mass;
    if (kinetic_energy > 1.0) return;
    
    //the formula used returns positive values if tmpEnergy>(approximately) 3keV, see Mu2e doc 1748-v1
    kinetic_energy = MeV2keV*kinetic_energy;
    double trackKine = kinetic_energy;
    kinetic_energy = std::log10(kinetic_energy);
    kinetic_energy *= _randGauss.fire(calorimeterCalibrations->LINpar2(cryId), calorimeterCalibrations->LINpar2Err(cryId) );
    kinetic_energy += _randGauss.fire(calorimeterCalibrations->LINpar1(cryId), calorimeterCalibrations->LINpar1Err(cryId) );
    kinetic_energy = std::log10(kinetic_energy);
    kinetic_energy *= _randGauss.fire(calorimeterCalibrations->LINpar0(cryId), calorimeterCalibrations->LINpar0Err(cryId) );
    kinetic_energy /= _randGauss.fire(calorimeterCalibrations->LINpar3(cryId), calorimeterCalibrations->LINpar3Err(cryId) );
    kinetic_energy /= std::log10(trackKine);
    
    if (kinetic_energy>0) energy *= kinetic_energy;		
    
    if (_debugLevel > 2) {
      std::cout<<"************************** BEFORE / AFTER NON-LINEARITY EFFECT-> edep_corr = "
	       << edep_save<<"  /  "<<energy<<std::endl
	       <<", energyKin = "
	       << trackKine 
	       << ", mass = "<< mass
	//	       << ", momentum.mag2() = "<< h.momentum().mag2()
	       <<std::endl;
    }
  }
  
  //-----------------------------------------------------------------------------
  void CaloDigisFromStepPointMCs::longitudinalResponseUniformityCorrection(double posZ, 
									   double cryhalflength, 
									   double& energy, 
									   int crid, 
									   ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations )
  {
    double edep_save(energy);
    
    posZ *= -_randGauss.fire(calorimeterCalibrations->LRUpar0(crid), calorimeterCalibrations->LRUpar0Err(crid) );
    posZ += 1.0;
    energy *= posZ;
    
    if (_debugLevel > 2) { 
      std::cout<<"***************BEFORE /  AFTER LRU EFFECT-> edep_corr = "
	       << edep_save<<"  /  "<<energy
	       << std::endl;	  
    }
  }


  //-----------------------------------------------------------------------------
  void CaloDigisFromStepPointMCs::printStepPointDebug(Calorimeter const & cal,StepPointMC const& h, int crid)
  {
    CLHEP::Hep3Vector testpos = h.position();
    std::cout<<"Reco "<<crid<<"   "<<h.volumeId()<<std::endl;
    std::cout<<"reco position Mu2e    "<<testpos<<std::endl;
    std::cout<<"reco position disk    "<<_calorimeter->toSectionFrame(_calorimeter->crystal(crid).sectionId(),testpos)<<std::endl;
    std::cout<<"reco position diskFF  "<<_calorimeter->toSectionFrameFF(_calorimeter->crystal(crid).sectionId(),testpos)<<std::endl;
    std::cout<<"reco position local   "<<_calorimeter->toCrystalFrame(crid,testpos)<<std::endl;
    std::cout<<"reco position disk    "<<_calorimeter->fromSectionFrame(_calorimeter->crystal(crid).sectionId(),_calorimeter->toSectionFrame(_calorimeter->crystal(crid).sectionId(),testpos))<<std::endl;
    std::cout<<"reco position diskFF  "<<_calorimeter->fromSectionFrameFF(_calorimeter->crystal(crid).sectionId(),_calorimeter->toSectionFrameFF(_calorimeter->crystal(crid).sectionId(),testpos))<<std::endl;
    std::cout<<"reco position local   "<<_calorimeter->fromCrystalFrame(crid,_calorimeter->toCrystalFrame(crid,testpos))<<std::endl;

    std::cout<<"reco Crystal orig     "<<_calorimeter->crystal(crid).position()<<std::endl;
    std::cout<<"reco Crystal orig sec "<<_calorimeter->crystal(crid).localPosition()<<std::endl;
    std::cout<<"Is inside I           "<<_calorimeter->isInsideCalorimeter(testpos)<<std::endl;
    std::cout<<"Is inside II          "<<_calorimeter->isInsideCalorimeter(testpos+CLHEP::Hep3Vector(0,-1,0))<<std::endl;
    std::cout<<"Is inside III         "<<_calorimeter->isInsideCalorimeter(testpos+CLHEP::Hep3Vector(0,0,-1))<<std::endl;
  }

} // end namespace mu2e

using mu2e::CaloDigisFromStepPointMCs;
DEFINE_ART_MODULE(CaloDigisFromStepPointMCs);



