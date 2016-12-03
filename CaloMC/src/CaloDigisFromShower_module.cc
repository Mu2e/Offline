//
// An EDProducer Module that reads Caloshower and produces CaloDigiPacked
//
// Fill the readout response from the deosited energy. Photo-statistic fluctuations are calculated independently 
// for each sensor and must be done at this stage.
//
// The output is split between the different digitization boards
//

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

#include "CaloMC/src/CaloPulseShape.cc" 
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/CaloShowerStepROCollection.hh"
#include "RecoDataProducts/inc/CaloDigiPackedCollection.hh"
#include "SeedService/inc/SeedService.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH2F.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>




namespace mu2e {


  class CaloDigisFromShower : public art::EDProducer {
  
    public:

      explicit CaloDigisFromShower(fhicl::ParameterSet const& pset) :

        // Parameters
        caloShowerModuleLabel_ (pset.get<std::string>("caloShowerModuleLabel")), 
        addNoise_              (pset.get<bool>       ("addNoise")),           //flag for adding or not Gaussian noise
        noise_                 (pset.get<double>     ("noise")),           // mV 
        thresholdVoltage_      (pset.get<double>     ("thresholdVoltage")), // mV 
        thresholdAmplitude_    (pset.get<double>     ("thresholdAmplitude")), //mV
        energyScale_           (pset.get<double>     ("energyScale")),       // mV/MeV
        digiSampling_          (pset.get<double>     ("digiSampling")),       // ns
        nBits_                 (pset.get<int>        ("nBits")),       // number of digitizer bits
        dynamicRange_          (pset.get<double>     ("dynamicRange")),       // mV
        endTimeBuffer_         (pset.get<double>     ("endTimeBuffer")),  // ns
        bufferDigi_            (pset.get<int>        ("bufferDigi")),  //# timestamps          
        pulseIntegralSteps_    (pset.get<int>        ("pulseIntegralSteps")),         //# integral steps
        diagLevel_             (pset.get<int>        ("diagLevel",0)),
        engine_                (createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
        randGauss_             (engine_),
	pulseShape_(CaloPulseShape(digiSampling_,pulseIntegralSteps_))
      {  

           produces<CaloDigiPackedCollection>();    
  
           maxADCCounts_ =  1 << nBits_;
           ADCTomV_      = dynamicRange_/maxADCCounts_;
           mVToADC_      = maxADCCounts_/dynamicRange_;

      }

      virtual ~CaloDigisFromShower() { }

      virtual void produce(art::Event& e);
      virtual void beginJob();
      virtual void beginRun(art::Run& aRun);




    private:

       std::string             caloShowerModuleLabel_;   

       double                  blindTime_;
       double                  mbtime_;

       bool                    addNoise_;
       double                  noise_;

       double                  thresholdVoltage_;         
       int                     thresholdAmplitude_;
       double                  energyScale_;
       double                  digiSampling_;
       int                     nBits_;
       double                  dynamicRange_;
       double                  endTimeBuffer_;
       int                     bufferDigi_;
       int                     pulseIntegralSteps_;

       int                     diagLevel_;
       CLHEP::HepRandomEngine& engine_;
       CLHEP::RandGaussQ       randGauss_;
       CaloPulseShape          pulseShape_;

       int                     maxADCCounts_;
       double                  ADCTomV_;
       double                  mVToADC_;
       const Calorimeter*      calorimeter_; 


       std::vector< std::vector<double> > pulseDigitized_;
       std::vector< std::vector<double> > waveforms_;


       //some diagnostic histograms
       TH1F*  hEdep_;
       TH1F*  hTime_;
       TH1F*  hNSamples_VsIro;
       TH1F*  hWFLength_;
       TH2F*  hWFLength_VsAmp;


       void   resetWaveforms();
       void   makeDigitization(const CaloShowerStepROCollection& caloShowerStepROs, CaloDigiPackedCollection&);
       void   fillWaveforms(const CaloShowerStepROCollection& caloShowerStepROs);
       void   readoutResponse(int ROID, double energyCorr, double time);
       void   buildOutputDigi(CaloDigiPackedCollection& caloDigis);
       void   fillOutoutRO(int iRO, int& nTotWords, std::vector<int>& caloDigiOutputs);
       int    findROCardIdx(int ROID);
       

  };


  
  //-----------------------------------------------------------------------------
  void CaloDigisFromShower::beginJob()
  {
       if ( diagLevel_ > 2)
       {
           art::ServiceHandle<art::TFileService> tfs;
           hEdep_          = tfs->make<TH1F>("hEdep","Hit energy deposition",            100,    0.,     50);
           hTime_          = tfs->make<TH1F>("hTime","Hit time ",                       4000,    0.,   2000);
           hNSamples_VsIro = tfs->make<TH1F>("hNSamplesVsIro","Number of samples /ro",  2000,    0,   2000);
           hWFLength_      = tfs->make<TH1F>("hWFLength","wavefrom length",              100,    0,    100);
           hWFLength_VsAmp = tfs->make<TH2F>("hWFLengthVsAmp","wavefrom length vs amp", 1000,    0, 100,2000, 0, 2000);
       }        
  }

  //-----------------------------------------------------------------------------
  void CaloDigisFromShower::beginRun(art::Run& aRun)
  {
      pulseShape_.buildShapes();
      if ( diagLevel_ > 3) pulseShape_.printShape();
  }



  //---------------------------------------------------------
  void CaloDigisFromShower::produce(art::Event& event)
  {   
      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigisFromShower::produce] begin" << std::endl;

      //update condition cache
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime_ = accPar->deBuncherPeriod;

      art::Handle<CaloShowerStepROCollection> caloShowerStepROHandle;
      event.getByLabel(caloShowerModuleLabel_, caloShowerStepROHandle);
      const CaloShowerStepROCollection& caloShowerStepROs(*caloShowerStepROHandle);


      std::unique_ptr<CaloDigiPackedCollection> caloDigis(new CaloDigiPackedCollection);    
      makeDigitization(caloShowerStepROs, *caloDigis);    
            
      event.put(std::move(caloDigis));

      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigisFromShower::produce] end" << std::endl;
  } 

  
  
  //-------------------------------------------------------------------------------------------------------------
  void CaloDigisFromShower::makeDigitization(const CaloShowerStepROCollection& caloShowerStepROs,CaloDigiPackedCollection& caloDigis)
  {
       mu2e::GeomHandle<mu2e::Calorimeter> ch;
       calorimeter_ = ch.get();                       

       resetWaveforms();
       fillWaveforms(caloShowerStepROs);
       buildOutputDigi(caloDigis);
  } 



  //-------------------------------------------------
  void CaloDigisFromShower::resetWaveforms()
  {
      unsigned int nWaveforms   = calorimeter_->nCrystal()*calorimeter_->caloGeomInfo().nROPerCrystal();
      unsigned int waveformSize = (mbtime_ - blindTime_ + endTimeBuffer_) / digiSampling_;  

      if (waveforms_.size() < nWaveforms || (waveforms_.size() && waveforms_[0].size() < waveformSize))
      {
          waveforms_.clear();
          for (unsigned int i=0;i< nWaveforms;++i) waveforms_.push_back(std::vector<double>(waveformSize,0.0)); 
      } 
      
      
      //if we want noise, fill all waveforms with noise, noise_ is in mV, sample is in counts, noise must be positive
      //otherwise simply reset the waveform to zero 
      
      if (addNoise_)
	for (unsigned int i=0; i<nWaveforms; ++i) std::generate(waveforms_[i].begin(),waveforms_[i].end(),
	                                                       [&]() {return  std::max(0.0,randGauss_.fire(0.0,noise_)*mVToADC_);});
      else 
      	for (unsigned int i=0; i<nWaveforms; ++i) std::fill(waveforms_[i].begin(),waveforms_[i].end(),0);
              
  }


  //----------------------------------------------------------------------------------------------------------
  void CaloDigisFromShower::fillWaveforms(const CaloShowerStepROCollection& caloShowerStepROs)
  {          
      double totalEdep(0);
      
      for (const auto& caloShowerStepRO : caloShowerStepROs)
      {       
          readoutResponse(caloShowerStepRO.ROID(),caloShowerStepRO.energy(),caloShowerStepRO.time());          
          if (diagLevel_ > 0) totalEdep +=  caloShowerStepRO.energy();
      }
            
      if (diagLevel_ > 1) std::cout<<"[CaloDigisFromShower::fillWaveforms] total energy processed "
                                   <<totalEdep/calorimeter_->caloGeomInfo().nROPerCrystal()<<std::endl;      
   }


  //-----------------------------------------------------------------------------------------
  void CaloDigisFromShower::readoutResponse(int ROID, double energyCorr, double time)
  {          
      double                     pulseAmp       = energyCorr*energyScale_;          
      double                     timeCorr       = time - blindTime_;          
      int                        startSample    = timeCorr/digiSampling_;
      int                        precisionIndex = (timeCorr/digiSampling_ - int(timeCorr/digiSampling_))*pulseIntegralSteps_;
      std::vector<double>&       waveform       = waveforms_.at(ROID);
      const std::vector<double>& pulse          = pulseShape_.pulseDigitized(precisionIndex);
      int                        stopSample     = std::min(startSample+pulse.size(), waveform.size());
     
      for (int timeSample = startSample; timeSample < stopSample; ++timeSample)
      {
          int     funcIdx   = timeSample - startSample;
          double  funcValue = pulse.at(funcIdx);
          double  wfAmp     = pulseAmp*funcValue;
          double  ADCCounts = wfAmp*mVToADC_;

          waveform.at(timeSample) += ADCCounts;  
          if ( waveform.at(timeSample) > maxADCCounts_) waveform.at(timeSample) = maxADCCounts_; 
      }
    
  }



  
 


  //----------------------------------------------------------------------------
  void CaloDigisFromShower::buildOutputDigi(CaloDigiPackedCollection& CaloDigis)
  {
        
      unsigned int Noutput(80);
      std::vector<int> nTotWords(Noutput,0);
      std::vector<std::vector<int>> caloDigiOutputs; 
      for (unsigned int i=0;i< Noutput;++i) caloDigiOutputs.push_back(std::vector<int>(1,0)); 



      for (unsigned int iRO=0; iRO<waveforms_.size(); ++iRO)
      {
          int outputIdx = findROCardIdx(iRO);
	  fillOutoutRO(iRO,nTotWords[outputIdx], caloDigiOutputs[outputIdx]);
      }

      for (unsigned int i=0;i< Noutput;++i) 
      {           
	  std::vector<int>& caloDigiOutput = caloDigiOutputs[i];
	  
	  caloDigiOutput[0] = nTotWords[i];          
	  if (nTotWords[i]>2) CaloDigis.emplace_back(CaloDigiPacked(caloDigiOutput));

	  if (diagLevel_ > 3)
	  {
              printf("[CaloDigisFromShower::buildOutputDigi] caloDigiOutput\n");
	      for (const auto& digi : caloDigiOutput) std::cout<<digi<<" "; std::cout<<std::endl;
	  }
      }

  }


  //--------------------------------------------------------------------
  // The output is written in the following data format:
  // nTotWords - nWords_roID - roiID - nWord_roID_ihit - time_roID_ihit - Data_roID_ihit - ...
  //  
  void CaloDigisFromShower::fillOutoutRO(int iRO, int& nTotWords, std::vector<int>& caloDigiOutput)
  {  
       int nRoWords = 2;     
       std::vector<int>  output; 
       std::vector<double>  &itWave = waveforms_.at(iRO);

       if (diagLevel_ > 4)
       {
           std::cout<<"CaloDigisFromShower::fillOutoutRO] Waveform content for readout "<<iRO<<std::endl; 
           for (const auto  &v : itWave) std::cout<<v<<" "; 
           std::cout<<std::endl;

           if (diagLevel_ > 5) std::cout<<"wfContent content (timesample: waveContent, funcValue)"<<std::endl;
       }

       output.push_back(0); //placeholder for number of words
       output.push_back(0); //placeholder for iRO

       int waveSize = itWave.size();
       int timeSample(0);
       while (timeSample < waveSize)
       {
           double waveContent = itWave.at(timeSample);
           double funcValue   = waveContent*ADCTomV_;

           if (diagLevel_ > 5 && waveContent > 0) printf("wfContent (%4i:  %4i, %9.3f) \n", timeSample, int(waveContent), funcValue);                    
           if (funcValue < thresholdVoltage_) {++timeSample; continue;}


           // find the starting / stopping point of the peak              
           // the stopping point is the first value below the threshold _and_ the buffer is also below the threshold
           int sampleStart = std::max(timeSample - bufferDigi_,0);
           int sampleStop  = timeSample;
           for (; sampleStop < waveSize; ++sampleStop)
           {
               int sampleCheck    = std::min(sampleStop+bufferDigi_+1,waveSize-1);
               double waveOverBuffer = *std::max_element(&itWave.at(sampleStop),&itWave.at(sampleCheck));
               if (waveOverBuffer*ADCTomV_ < thresholdVoltage_) break;
           }             
           sampleStop = std::min(sampleStop + bufferDigi_, waveSize-1);


           //forward the scanning time
           timeSample = sampleStop+1;


           //check if peak is acceptable and digitize
           if (sampleStop == sampleStart)                continue;

           double sampleMax = *std::max_element(&itWave.at(sampleStart),&itWave.at(sampleStop));
           if (sampleMax*ADCTomV_ < thresholdAmplitude_) continue;


           int nHitWords = 2 + sampleStop - sampleStart +1;
           int t0        = int(sampleStart*digiSampling_);

           output.emplace_back(nHitWords);
           output.emplace_back(t0 + blindTime_);
           for (int i=sampleStart; i<=sampleStop; ++i) output.push_back(int(itWave.at(i)));

           nRoWords += nHitWords;

           if (diagLevel_ > 2 && t0 > blindTime_)
           {
                hNSamples_VsIro->Fill(iRO, sampleStop - sampleStart);
                hWFLength_->Fill(sampleStop - sampleStart);
                hWFLength_VsAmp->Fill(sampleStop - sampleStart, sampleMax);
           }
       }

       output[0] = nRoWords;
       output[1] = iRO;

       if (diagLevel_ > 3 && output.size() > 2)
       {
           std::cout<<"CaloDigisFromShower::fillOutoutRO] Readout "<<iRO<<std::endl;
           for (const auto& val : output) std::cout<<val<<" ";std::cout<<std::endl;                   
       }        

       //comment if you want to keep all readouts
       if (nRoWords==2) return; 
       
       nTotWords += nRoWords;
       caloDigiOutput.insert(caloDigiOutput.end(), output.begin(), output.end());       
  }




  //----------------------------------------------
  int CaloDigisFromShower::findROCardIdx(int iRO)
  {
      return iRO / 40;
  }



}

using mu2e::CaloDigisFromShower;
DEFINE_ART_MODULE(CaloDigisFromShower);
