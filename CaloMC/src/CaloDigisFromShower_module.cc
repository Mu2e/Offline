//
// An EDProducer Module that reads Caloshower and produces CaloDigiPacked
//
// Fill the readout response from the deosited energy. Photo-statistic fluctuations are calculated independently 
// for each sensor and must be done at this stage.
//
// Original author B. Echenard
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
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerCollection.hh"
#include "RecoDataProducts/inc/CaloDigiPacked.hh"
#include "SeedService/inc/SeedService.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoisson.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>

#include "TH2F.h"
#include "TFile.h"



namespace mu2e {


  class CaloDigisFromShower : public art::EDProducer {
  
    public:

      explicit CaloDigisFromShower(fhicl::ParameterSet const& pset) :

        // Parameters
        caloShowerModuleLabel_ (pset.get<std::string>("caloShowerModuleLabel")), 
        blindTime_             (pset.get<double>     ("blindTime")),         // ns
        caloPEStatCorrection_  (pset.get<bool>       ("caloPEStatcorrection")),
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
        randPoisson_           (engine_),
        randGauss_             (engine_),
	pulseShape_(CaloPulseShape(digiSampling_,pulseIntegralSteps_))
      {  

           produces<CaloDigiPacked>();    
  
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

       bool                    caloPEStatCorrection_;
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
       CLHEP::RandPoisson      randPoisson_;
       CLHEP::RandGaussQ       randGauss_;
       CaloPulseShape          pulseShape_;

       int                     maxADCCounts_;
       double                  ADCTomV_;
       double                  mVToADC_;
       const Calorimeter*      calorimeter_; 


       std::vector< std::vector<double> > pulseDigitized_;
       std::vector< std::vector<int> >    waveforms_;


       //some diagnostic histograms
       TH1F*  hEdep_;
       TH1F*  hTime_;
       TH2F*  hPECorr_;
       TH1F*  hPECorr2_;
       TH1F*  hNSamples_;
       TH1F*  hNHits_;   
       TH1F*  hNSamples_VsIro;
       TH1F*  hWFLength_;
       TH2F*  hWFLength_VsAmp;


       void   resetWaveforms();
       void   makeDigitization(const art::Handle<CaloShowerCollection>& CaloShowersHandle, CaloDigiPacked&);
       void   fillWaveforms(const art::Handle<CaloShowerCollection>& crystalStepsHandles);
       void   readoutResponse(double Edep, double Time, int CrId, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);
       void   buildOutputDigi(CaloDigiPacked& caloDigis);
       double photoStatisticsCorrection(int crystalId, double edep, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);

  };




  //------------------------------------------
  // Calculate the normalized waveform content for different starting time
  // The t0 is set when the pulse crosses the 1e-3 level, fill all previou bins with zero
  // instead of shifting the waveform forward for each time step and integrating between the digitizer time boundaries, 
  // we can shift the start time backward and integrate over a fixed time to get the same effect
  // Assuming T0 = 180, set all previous bins to zero. 
  // Now integrating from 180-185, 185-190,... with a waveform shifted at 182 is the same as integrating from 178-183, 183-188,....
  // (a drawing helps) and the second solution is simpler to implement (see comment in text)
  
  void CaloDigisFromShower::beginJob()
  {
       art::ServiceHandle<art::TFileService> tfs;
       hEdep_          = tfs->make<TH1F>("hEdep","Hit energy deposition",            100,    0.,     50);
       hTime_          = tfs->make<TH1F>("hTime","Hit time ",                       4000,    0.,   2000);
       hPECorr_        = tfs->make<TH2F>("hPECorr", "PE stat correction",            100,    0., 50, 100,    0,  50);
       hPECorr2_       = tfs->make<TH1F>("hPECorr2","PE stat correction",            100,  -0.5,    0.5);
       hNSamples_      = tfs->make<TH1F>("hNSamples","Numerb of samples / evt",    10000,   5e3,    5e3);
       hNHits_         = tfs->make<TH1F>("hNHits","Numer of hits / evt",           10000,     0., 10000);
       hNSamples_VsIro = tfs->make<TH1F>("hNSamplesVsIro","Number of samples /ro",  2000,     0,   2000);
       hWFLength_      = tfs->make<TH1F>("hWFLength","wavefrom length",              100,     0,    100);
       hWFLength_VsAmp = tfs->make<TH2F>("hWFLengthVsAmp","wavefrom length vs amp", 1000,     0, 100,2000, 0, 2000);        
  }

  //-----------------------------------------------------------------------------
  void CaloDigisFromShower::beginRun(art::Run& aRun)
  {
      pulseShape_.buildShapes(); 
  }



  //---------------------------------------------------------
  void CaloDigisFromShower::produce(art::Event& event)
  {   
      if ( diagLevel_ > 0 ) std::cout << "[CaloDigisFromShower:] produce() begin" << std::endl;

      //update condition cache
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime_ = accPar->deBuncherPeriod;

      art::Handle<CaloShowerCollection> caloShowerHandle;
      event.getByLabel(caloShowerModuleLabel_, caloShowerHandle);

      std::unique_ptr<CaloDigiPacked> caloDigis(new CaloDigiPacked);    
      makeDigitization(caloShowerHandle, *caloDigis);    
      event.put(std::move(caloDigis));

      if ( diagLevel_ > 0 ) std::cout << "[CaloDigisFromShower:] produce() end" << std::endl;
  } 

  
  
  //-------------------------------------------------------------------------------------------------------------
  void CaloDigisFromShower::makeDigitization(const art::Handle<CaloShowerCollection>& CaloShowerCollHandle,CaloDigiPacked& caloDigis)
  {
       mu2e::GeomHandle<mu2e::Calorimeter> ch;
       calorimeter_ = ch.get();                       

       resetWaveforms();
       fillWaveforms(CaloShowerCollHandle);
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
          for (unsigned int i=0;i< nWaveforms;++i) waveforms_.push_back(std::vector<int>(waveformSize,0)); 
      } 
      
      //if we want noise, fill all waveforms with noise, noise_ is in mV, sample is in counts, noise must be positive
      //otherwise simply reset the waveform to zero 

      if (addNoise_)
	for (unsigned int i=0; i<nWaveforms; ++i) std::generate(waveforms_[i].begin(),waveforms_[i].end(),[&]() {return  int(std::max(0.0,randGauss_.fire(0.0,noise_)*mVToADC_));});
      else 
      	for (unsigned int i=0; i<nWaveforms; ++i) std::fill(waveforms_[i].begin(),waveforms_[i].end(),0);
              
  }


  //----------------------------------------------------------------------------------------------------------
  void CaloDigisFromShower::fillWaveforms(const art::Handle<CaloShowerCollection>& CaloShowerCollHandle)
  {
          
      double totalEdep(0);
      
      GlobalConstantsHandle<ParticleDataTable> pdt;
      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
      const CaloShowerCollection& caloShowers(*CaloShowerCollHandle);

      
      for (unsigned int i=0; i < caloShowers.size(); ++i)
      {
           const CaloShower& shower =  caloShowers.at(i);
           if ( shower.energy() < 1e-4 ) continue;
           
           double energyDep = shower.energy();
           int    crystalId = shower.crystalId();
           double hitTime   = shower.time();
                                 
           if (diagLevel_ > 3) printf("[CaloDigisFromShower::fillWaveforms] HitTime = %4.2f  Energy = %4.2f Cry = %4i\n", 
                                      hitTime, energyDep, crystalId);
           
           if (hitTime < blindTime_) continue;
           if (hitTime > mbtime_ )   continue;
                        
           readoutResponse(energyDep, hitTime, crystalId, calorimeterCalibrations);

           if (diagLevel_ > 0) totalEdep += energyDep; 
           if (diagLevel_ > 2)
           {                   
               hEdep_->Fill(energyDep);
               hTime_->Fill(hitTime);
           }
       }
       
       if (diagLevel_ > 1) std::cout<<"[CaloDigisFromShower] total energy processed "<<totalEdep<<std::endl;
       
   }


  //--------------------------------------------------------------------
  // 2015-09-14 Gianipez and L. Morscalchi:
  // Propagate the steppoint energy deposition to the photosensor accoridng to optical transportation
  // Ray tracing is not implemented yet, so the energy deposition is evenly splitted between the two photosensors
  //
  void CaloDigisFromShower::readoutResponse(double energyDep, double time, int crystalId, 
                                            ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {          

      int ROidBase = calorimeter_->ROBaseByCrystal(crystalId);
      int nROs     = calorimeter_->caloGeomInfo().nROPerCrystal();

      for (int i=0; i<nROs; ++i)
      {
          int ROId = ROidBase + i;

          double energyCorr(energyDep);
          if (caloPEStatCorrection_)
                 energyCorr = photoStatisticsCorrection(crystalId, energyDep, calorimeterCalibrations);          
                    
          double                     pulseAmp       = energyCorr*energyScale_;          
          double                     timeCorr       = time - blindTime_;          
          int                        startSample    = timeCorr/digiSampling_;
          int                        precisionIndex = (timeCorr/digiSampling_ - int(timeCorr/digiSampling_))*pulseIntegralSteps_;
          std::vector<int>&          waveform       = waveforms_.at(ROId);
          const std::vector<double>& pulse          = pulseShape_.pulseDigitized(precisionIndex);
          int                        stopSample     = std::min(startSample+pulse.size(), waveform.size());

          if (diagLevel_ > 1 )
          {
               printf("[CaloDigisFromShower::makeCalorimeterHits]   Readout %4i   eDep = %9.3f MeV   amplitude = %9.3f mV   time = %9.3f ns\n", 
                       ROId, energyCorr, pulseAmp, time);
               
               if (diagLevel_ > 3 )
                   printf("[CaloDigisFromShower::makeCalorimeterHits]   timeSample   |   ADC-counts   |    wfAmp    |   signalAmp  | Total ADC \n");
               
               if (diagLevel_ > 2)
               {
                   hPECorr_->Fill(energyCorr,energyDep);
                   hPECorr2_->Fill((energyCorr-energyDep)/sqrt(energyDep));
               }               
          }


          for (int timeSample = startSample; timeSample < stopSample; ++timeSample)
          {
              int     funcIdx   = timeSample - startSample;
              double  funcValue = pulse.at(funcIdx);
              double  wfAmp     = pulseAmp*funcValue;
              int     ADCCounts = int(wfAmp*mVToADC_);

              waveform.at(timeSample) += ADCCounts;  
              if ( waveform.at(timeSample) > maxADCCounts_) waveform.at(timeSample) = maxADCCounts_; 

              if (diagLevel_ > 3 && wfAmp > 1e-6 )
              {
                  printf("[CaloDigisFromShower::makeCalorimeterHits]   %4i   |   %4i   |    %9.3f  | %9.3f | %4i \n", 
                         timeSample, ADCCounts,  wfAmp, funcValue, waveform.at(timeSample));
              }
          }

      }           
    
  }



  
 


  //--------------------------------------------------------------------
  // The output is written in the following data format:
  // nTotWords - nWords_roID - roiID - nWord_roID_ihit - time_roID_ihit - Data_roID_ihit - ...
  //  
  void CaloDigisFromShower::buildOutputDigi(CaloDigiPacked& CaloDigis)
  {
  
      int nHits(0), nSamples(0), nTotWords(1);

      std::vector<int> caloDigiOutput; 
      caloDigiOutput.push_back(0);       //reserve nTotword word, will be filled at the end


      for (unsigned int iRO=0; iRO<waveforms_.size(); ++iRO)
      {
           int nRoWords = 2;     
           std::vector<int>  output; 
           std::vector<int>  &itWave = waveforms_.at(iRO);


           if (diagLevel_ > 1)
           {
               std::cout<<"CaloDigisFromShower::buildOutputDigi] Waveform content for readout "<<iRO<<std::endl; 
               for (const auto  &v : itWave) std::cout<<v<<" "; 
               std::cout<<std::endl;

               if (diagLevel_ > 3) std::cout<<"wfContent content (timesample: waveContent, funcValue)"<<std::endl;
           }


           output.push_back(0); //placeholder for number of words
           output.push_back(0); //placeholder for iRO

           int waveSize = itWave.size();
           int timeSample(0);
           while (timeSample < waveSize)
           {
               int    waveContent = itWave.at(timeSample);
               double funcValue   = waveContent*ADCTomV_;

               if (diagLevel_ > 3 && waveContent > 0) printf("wfContent (%4i:  %4i, %9.3f) \n", timeSample, waveContent, funcValue);                    
               if (funcValue < thresholdVoltage_) {++timeSample; continue;}


               // find the starting / stopping point of the peak              
               // the stopping point is the first value below the threshold _and_ the buffer is also below the threshold
               int sampleStart = std::max(timeSample - bufferDigi_,0);
               int sampleStop  = timeSample;
               for (; sampleStop < waveSize; ++sampleStop)
               {
                   int sampleCheck    = std::min(sampleStop+bufferDigi_+1,waveSize-1);
                   int waveOverBuffer = *std::max_element(&itWave.at(sampleStop),&itWave.at(sampleCheck));
                   if (waveOverBuffer*ADCTomV_ < thresholdVoltage_) break;
               }             
               sampleStop = std::min(sampleStop + bufferDigi_, waveSize-1);


               //forward the scanning time
               timeSample = sampleStop+1;


               //check if peak is acceptable and digitize
               if (sampleStop == sampleStart)                continue;

               int sampleMax = *std::max_element(&itWave.at(sampleStart),&itWave.at(sampleStop));
               if (sampleMax*ADCTomV_ < thresholdAmplitude_) continue;


               int nHitWords = 2 + sampleStop - sampleStart +1;
               int t0        = int(sampleStart*digiSampling_);

               output.emplace_back(nHitWords);
               output.emplace_back(t0 + blindTime_);
               for (int i=sampleStart; i<=sampleStop; ++i) output.push_back(itWave.at(i));

               nRoWords += nHitWords;

               if (diagLevel_ > 2 && t0 > blindTime_)
               {
                    ++nHits;
                    nSamples += (sampleStop - sampleStart);

                    hNSamples_VsIro->Fill(iRO, sampleStop - sampleStart);
                    hWFLength_->Fill(sampleStop - sampleStart);
                    hWFLength_VsAmp->Fill(sampleStop - sampleStart, sampleMax);
               }
           }

           nTotWords += nRoWords;

           output[0] = nRoWords;
           output[1] = iRO;

           caloDigiOutput.insert(caloDigiOutput.end(), output.begin(), output.end());

           if (diagLevel_ > 3)
           {
               std::cout<<"CaloDigisFromShower::buildOutputDigi] Readout "<<iRO<<std::endl;
               for (const auto& val : output) std::cout<<val<<" ";std::cout<<std::endl;                   
           }        
       }


      caloDigiOutput[0] = nTotWords; 



      if (diagLevel_  > 2)
      {
          hNSamples_->Fill(nSamples);
          hNHits_->Fill(nHits);
      }
      if (diagLevel_ > 3)
      {
          printf("[CaloDigisFromShower::buildOutputDigi] caloDigiOutput\n");
          for (const auto& digi : caloDigiOutput) std::cout<<digi<<" "; std::cout<<std::endl;
      }

      
      
      CaloDigis = std::move(CaloDigiPacked(caloDigiOutput)); 
  }



  //-----------------------------------------------------------------------------
  double CaloDigisFromShower::photoStatisticsCorrection(int crystalId, double edep, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  { 
      double lightYield = randPoisson_.fire(edep*calorimeterCalibrations->peMeV(crystalId));
      return lightYield/calorimeterCalibrations->peMeV(crystalId);
  }


}

using mu2e::CaloDigisFromShower;
DEFINE_ART_MODULE(CaloDigisFromShower);
