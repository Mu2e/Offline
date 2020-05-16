//
// An EDProducer Module that reads Caloshower and produces CaloDigiPacked
//
// Fill the readout response from the deosited energy. Photo-statistic fluctuations are calculated independently
// for each sensor and must be done at this stage.
//
// The output is split between the different digitization boards
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "CaloMC/inc/CaloPulseShape.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/CaloShowerStepROCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "SeedService/inc/SeedService.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandPoissonQ.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>


namespace mu2e {


  class CaloDigiFromShower : public art::EDProducer 
  {
     public:
         
         struct Config 
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             
             fhicl::Atom<art::InputTag> caloShowerCollection { Name("caloShowerROCollection"), Comment("CaloShowerRO collection name") }; 
             fhicl::Atom<double>        blindTime            { Name("blindTime"),              Comment("Microbunch blind time") }; 
             fhicl::Atom<bool>          addNoise             { Name("addNoise"),               Comment("Add noise to waveform") }; 
             fhicl::Atom<double>        noise                { Name("noise"),                  Comment("Noise level - ADC equivalent") }; 
             fhicl::Atom<double>        thresholdVoltage     { Name("thresholdVoltage"),       Comment("Threshold voltage for sigitizing") }; 
             fhicl::Atom<int>           thresholdAmplitude   { Name("thresholdAmplitude"),     Comment("Threshold amplitude for keeping hit") }; 
             fhicl::Atom<double>        digiSampling         { Name("digiSampling"),           Comment("Digitization time sampling") }; 
             fhicl::Atom<int>           nBits                { Name("nBits"),                  Comment("ADC Number of bits") }; 
             fhicl::Atom<double>        dynamicRange         { Name("dynamicRange"),           Comment("ADC dynamic range") }; 
             fhicl::Atom<double>        endTimeBuffer        { Name("endTimeBuffer"),          Comment("Number of extra timestamps after end of pulse") }; 
             fhicl::Atom<int>           bufferDigi           { Name("bufferDigi"),             Comment("Number of timeStamps for the buffer digi") }; 
             fhicl::Atom<int>           pulseIntegralSteps   { Name("pulseIntegralSteps"),     Comment("Numer of time sub-division for CaloPulseChape") }; 
             fhicl::Atom<int>           nBinsPeak            { Name("nBinsPeak"),              Comment("XXXXXXXXX") }; 
             fhicl::Atom<int>           diagLevel            { Name("diagLevel"),              Comment("Diag Level"),0 };
         };
         

         explicit CaloDigiFromShower(const art::EDProducer::Table<Config>& config) :
            EDProducer{config},
            caloShowerToken_{consumes<CaloShowerStepROCollection>(config().caloShowerCollection())},
            blindTime_         (config().blindTime()),
            addNoise_          (config().addNoise()),
            noise_             (config().noise()),
            thresholdVoltage_  (config().thresholdVoltage()),
            thresholdAmplitude_(config().thresholdAmplitude()),
            digiSampling_      (config().digiSampling()),
            maxADCCounts_      (1 << config().nBits()),
            dynamicRange_      (config().dynamicRange()),
            endTimeBuffer_     (config().endTimeBuffer()),
            bufferDigi_        (config().bufferDigi()),
            pulseIntegralSteps_(config().pulseIntegralSteps()),
            nBinsPeak_         (config().nBinsPeak()),
            diagLevel_         (config().diagLevel()),
            engine_            (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
            randPoisson_       (engine_),
            pulseShape_        (CaloPulseShape(digiSampling_,pulseIntegralSteps_))
         {
             produces<CaloDigiCollection>();

             //NEEDS TO BE REFINED
             ADCTomV_      = dynamicRange_/float(maxADCCounts_);
         }

    private:
       void produce(art::Event& e) override;
       void beginRun(art::Run& aRun) override;
       
       void generateNoise(std::vector<std::vector<double>>& waveforms);
       void makeDigitization(const CaloShowerStepROCollection& caloShowerStepROs, CaloDigiCollection&);
       void fillWaveforms(std::vector<std::vector<double>>& waveforms, const CaloShowerStepROCollection& caloShowerStepROs,
                          const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);
       void buildOutputDigi(std::vector<std::vector<double>>& waveforms, CaloDigiCollection& caloDigiColl);
       void diag0(int iRO,std::vector<double>& itWave );
       void diag1(int iRO, double time, std::vector<int>& wf );

       
       const art::ProductToken<CaloShowerStepROCollection> caloShowerToken_;
       double                  blindTime_;
       double                  mbtime_;
       bool                    addNoise_;
       double                  noise_;
       double                  thresholdVoltage_;
       int                     thresholdAmplitude_;
       double                  digiSampling_;
       int                     maxADCCounts_;
       double                  dynamicRange_;
       double                  endTimeBuffer_;
       int                     bufferDigi_;
       int                     pulseIntegralSteps_;
       size_t  		       nBinsPeak_;
       int                     diagLevel_;
       CLHEP::HepRandomEngine& engine_;
       CLHEP::RandPoissonQ     randPoisson_;
       CaloPulseShape          pulseShape_;

       double                  MeVToADC_;
       double                  ADCTomV_;
       const Calorimeter*      calorimeter_;
  };


  //-----------------------------------------------------------------------------
  void CaloDigiFromShower::beginRun(art::Run& aRun)
  {
      pulseShape_.buildShapes();
      //if ( diagLevel_ > 3) pulseShape_.printShape();
  }



  //---------------------------------------------------------
  void CaloDigiFromShower::produce(art::Event& event)
  {
      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiFromShower::produce] begin" << std::endl;

      //update condition cache
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime_ = accPar->deBuncherPeriod;
      
      auto caloShowerStepHandle     = event.getValidHandle(caloShowerToken_);
      const auto& caloShowerStepROs = *caloShowerStepHandle;
      auto caloDigiColl             = std::make_unique<CaloDigiCollection>();
      
      makeDigitization(caloShowerStepROs, *caloDigiColl);

      event.put(std::move(caloDigiColl));

      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiFromShower::produce] end" << std::endl;
  }

  
  //-----------------------------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::makeDigitization(const CaloShowerStepROCollection& caloShowerStepROs,CaloDigiCollection& caloDigiColl)
  {
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      calorimeter_ = ch.get();

      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
      
      unsigned nWaveforms   = calorimeter_->nCrystal()*calorimeter_->caloInfo().nROPerCrystal();
      unsigned waveformSize = (mbtime_ - blindTime_ + endTimeBuffer_) / digiSampling_;     
      std::vector<std::vector<double>> waveforms(nWaveforms,std::vector<double>(waveformSize,0.0));
      
      if (addNoise_) generateNoise(waveforms);
      fillWaveforms(waveforms, caloShowerStepROs,calorimeterCalibrations);
      buildOutputDigi(waveforms, caloDigiColl);
  }


  //-------------------------------------------------
  void CaloDigiFromShower::generateNoise(std::vector<std::vector<double>>& waveforms)
  {
     for (auto& waveform : waveforms)
     {
        std::generate(waveform.begin(),waveform.end(),[this]{return std::max(randPoisson_.fire(noise_)-noise_,0.0);});
     }
  }


  //----------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::fillWaveforms(std::vector<std::vector<double>>& waveforms, const CaloShowerStepROCollection& caloShowerStepROs,
                                         const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
      unsigned totalPE(0);
      
      for (const auto& caloShowerStepRO : caloShowerStepROs)
      {
          int ROID = caloShowerStepRO.ROID();
          auto& waveform = waveforms.at(ROID);
          for (const float PEtime : caloShowerStepRO.PETime())
          {        
              float       time           = PEtime - blindTime_;         
              int         startSample    = time/digiSampling_;
              const auto& pulse          = pulseShape_.digitizedPulse(time);
              int         stopSample     = std::min(startSample+pulse.size(), waveform.size());

              for (int timeSample = startSample; timeSample < stopSample; ++timeSample) 
                 waveform.at(timeSample) += pulse.at(timeSample - startSample)*calorimeterCalibrations->peMeV(ROID)*MeVToADC_;
          }
          totalPE += caloShowerStepRO.NPE();
      }
      
      //for (auto& waveform : waveforms)
      //   std::transform(waveform.begin(),waveform.end(),waveform.begin(),[this](double val){return std::min(val,double(maxADCCounts_));});
      
      if (diagLevel_ > 1) std::cout<<"[CaloDigiFromShower::fillWaveforms] total PE processed "<<totalPE<<std::endl;
  }




  //----------------------------------------------------------------------------
  void CaloDigiFromShower::buildOutputDigi(std::vector<std::vector<double>>& waveforms, CaloDigiCollection& caloDigiColl)
  {

    for (unsigned int iRO=0; iRO<waveforms.size(); ++iRO)
      {

        if (diagLevel_ > 5) std::cout<<"wfContent content (timesample: waveContent, funcValue)"<<std::endl;
        if (diagLevel_ > 4) diag0(iRO,waveforms.at(iRO));


        std::vector<double>& itWave = waveforms.at(iRO);

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
                int sampleCheck = std::min(sampleStop+bufferDigi_+1,waveSize-1);
                double waveOverBuffer = *std::max_element(&itWave.at(sampleStop),&itWave.at(sampleCheck));
                if (waveOverBuffer*ADCTomV_ < thresholdVoltage_) break;
            }
            
            sampleStop = std::min(sampleStop + bufferDigi_, waveSize-1);

            timeSample = sampleStop+1;  //forward the scanning time


            if (sampleStop == sampleStart) continue;  //check if peak is acceptable and digitize

		
		float wfInt(0);
		size_t peakP(0);
		for(size_t i =sampleStart; i<sampleStop-nBinsPeak_;++i){

			float sum(0);
			for(size_t j=0; j< nBinsPeak_;++j){
				sum+=itWave.at(i+j);
			}
			if(sum>wfInt){
				wfInt = sum;
				peakP = i+nBinsPeak_/2;
			}


		}
			
		double sampleMax = itWave.at(peakP);
		peakP = peakP - sampleStart;
            	if (sampleMax*ADCTomV_ < thresholdAmplitude_) continue;


            int t0 = int(sampleStart*digiSampling_+ blindTime_);
            std::vector<int> wf;
            for (int i=sampleStart; i<=sampleStop; ++i) wf.push_back(int(itWave.at(i)));

            caloDigiColl.emplace_back( CaloDigi(iRO,t0,wf, peakP) );

            if (diagLevel_ > 4) diag1(iRO,t0,wf);
          }
      }
  }




  void CaloDigiFromShower::diag0(int iRO,std::vector<double>& itWave )
  {
     if (*std::max_element(itWave.begin(),itWave.end())<1) return;
     std::cout<<"CaloDigiFromShower::fillOutoutRO] Waveform content for readout "<<iRO<<std::endl;
     for (const auto  &v : itWave) std::cout<<v<<" ";
     std::cout<<std::endl;
  }

  void CaloDigiFromShower::diag1(int iRO, double time, std::vector<int>& wf )
  {
     std::cout<<"Created caloDigi with roID = "<<iRO<<"  time="<<time<<" and content ";
     for (const auto  &v : wf) std::cout<<v<<" ";
     std::cout<<std::endl;
  }

}

DEFINE_ART_MODULE(mu2e::CaloDigiFromShower);
