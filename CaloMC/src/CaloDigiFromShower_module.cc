//
// An EDProducer Module that reads CaloshowerStepROs and produces the readout waveform.
// Photo-statistic fluctuations are calculated independently for each sensor.
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
#include <numeric>


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
             fhicl::Atom<double>        noise                { Name("noise"),                  Comment("Noise level - PE equivalent") }; 
             fhicl::Atom<double>        digiSampling         { Name("digiSampling"),           Comment("Digitization time sampling") }; 
             fhicl::Atom<int>           nBits                { Name("nBits"),                  Comment("ADC Number of bits") }; 
             fhicl::Atom<int>           nBinsPeak            { Name("nBinsPeak"),              Comment("Window size for finding local maximum to digitize wf") }; 
             fhicl::Atom<double>        MeVToADC             { Name("MeVToADC"),               Comment("MeV to ADC conversion factor") }; 
             fhicl::Atom<int>           minPeakADC           { Name("minPeakADC"),             Comment("Minimum ADC hits of local peak to digitize") }; 
             fhicl::Atom<double>        endTimeBuffer        { Name("endTimeBuffer"),          Comment("Number of extra timestamps after end of pulse") }; 
             fhicl::Atom<int>           bufferDigi           { Name("bufferDigi"),             Comment("Number of timeStamps for the buffer digi") }; 
             fhicl::Atom<int>           pulseIntegralSteps   { Name("pulseIntegralSteps"),     Comment("Numer of time sub-division for CaloPulseChape") }; 
             fhicl::Atom<int>           diagLevel            { Name("diagLevel"),              Comment("Diag Level"),0 };
         };
         
         explicit CaloDigiFromShower(const art::EDProducer::Table<Config>& config) :
            EDProducer{config},
            caloShowerToken_{consumes<CaloShowerStepROCollection>(config().caloShowerCollection())},
            blindTime_         (config().blindTime()),
            addNoise_          (config().addNoise()),
            noise_             (config().noise()),
            digiSampling_      (config().digiSampling()),
            bufferDigi_        (config().bufferDigi()),
            nBinsPeak_         (config().nBinsPeak()),
            minPeakADC_        (config().minPeakADC()),
            MeVToADC_          (config().MeVToADC()),
            maxADCCounts_      (1 << config().nBits()),
            endTimeBuffer_     (config().endTimeBuffer()),
            pulseIntegralSteps_(config().pulseIntegralSteps()),
            diagLevel_         (config().diagLevel()),
            engine_            (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
            randPoisson_       (engine_),
            pulseShape_        (CaloPulseShape(digiSampling_,pulseIntegralSteps_,false))
         {
             produces<CaloDigiCollection>();
         }

    private:
       void produce(art::Event& e)   override;
       void beginRun(art::Run& aRun) override;
       
       void generateNoise(std::vector<std::vector<double>>& waveforms, const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);
       void makeDigitization(const CaloShowerStepROCollection& caloShowerStepROs, CaloDigiCollection&);
       void fillWaveforms(std::vector<std::vector<double>>& waveforms, const CaloShowerStepROCollection& caloShowerStepROs,
                          const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);
       void buildOutputDigi(std::vector<std::vector<double>>& waveforms, CaloDigiCollection& caloDigiColl);
       void diag0(int iRO, const std::vector<double>& wf);
       void diag1(int iRO, double time, size_t peakP, const std::vector<int>& wf);
       void diag2(const CaloDigiCollection& caloDigiColl);

       
       const art::ProductToken<CaloShowerStepROCollection> caloShowerToken_;
       double                  blindTime_;
       double                  mbtime_;
       bool                    addNoise_;
       double                  noise_;
       double                  digiSampling_;
       int                     bufferDigi_;
       int  		       nBinsPeak_;
       int  		       minPeakADC_;
       double                  MeVToADC_;
       int                     maxADCCounts_;
       double                  endTimeBuffer_;
       int                     pulseIntegralSteps_;
       int                     diagLevel_;
       CLHEP::HepRandomEngine& engine_;
       CLHEP::RandPoissonQ     randPoisson_;
       CaloPulseShape          pulseShape_;
       const Calorimeter*      calorimeter_;
  };


  //-----------------------------------------------------------------------------
  void CaloDigiFromShower::beginRun(art::Run& aRun)
  {
      pulseShape_.buildShapes();
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
      
      if (addNoise_) generateNoise(waveforms,calorimeterCalibrations);
      fillWaveforms(waveforms, caloShowerStepROs,calorimeterCalibrations);
      buildOutputDigi(waveforms, caloDigiColl);
      
      if (diagLevel_ > 1) diag2(caloDigiColl);
  }


  //-------------------------------------------------
  void CaloDigiFromShower::generateNoise(std::vector<std::vector<double>>& waveforms, const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
     for (unsigned i=0; i<waveforms.size(); ++i)
     {
        auto& waveform = waveforms[i];
        double scaleFactor  = MeVToADC_/calorimeterCalibrations->peMeV(i);
        auto noiseGenerator = [this,scaleFactor]{return std::max( (randPoisson_.fire(noise_)-noise_),0.0)*scaleFactor;};
        
        std::generate(waveform.begin(),waveform.end(),noiseGenerator);
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
              {
                 waveform.at(timeSample) += pulse.at(timeSample - startSample)/calorimeterCalibrations->peMeV(ROID)*MeVToADC_;
              }
          }
          totalPE += caloShowerStepRO.NPE();
      }
            
      if (diagLevel_ > 1) std::cout<<"[CaloDigiFromShower::fillWaveforms] total PE processed "<<totalPE<<std::endl;
  }



  //----------------------------------------------------------------------------
  void CaloDigiFromShower::buildOutputDigi(std::vector<std::vector<double>>& waveforms, CaloDigiCollection& caloDigiColl)
  {
    float totEdepEq(0);
    for (unsigned int iRO=0; iRO<waveforms.size(); ++iRO)
    {
         const auto& wf = waveforms[iRO];
         int wfSize = wf.size();

         if (diagLevel_ > 2) diag0(iRO,wf);

         int timeSample(nBinsPeak_);
         while (timeSample < wfSize)
         {
              // find the local maximum over a window of size peakWindow_ and above threshold value
              if (wf[timeSample] < minPeakADC_) {++timeSample; continue;}
              if (std::max_element(&wf[timeSample]-nBinsPeak_,&wf[timeSample]+nBinsPeak_) != &wf[timeSample]) {++timeSample; continue;}

              // find the starting / stopping point of the peak
              // the stopping point is the first value below the threshold

              int sampleStart = std::max(timeSample - bufferDigi_,0);
              int sampleStop  = timeSample;
              while (sampleStop < wfSize)
              {
                  if (wf[sampleStop]<minPeakADC_) break;
                  ++sampleStop;
              }

              timeSample = sampleStop+1; //forward the scanning time                       
              if (diagLevel_ > 2) std::cout<<"[CaloDigiFromShower] found peak with startSample="<<sampleStart<<"  stopSample="<<sampleStop<<std::endl;
   
              //build the digi, find t0, PeakPosition, wf with saturation
              int t0 = int(sampleStart*digiSampling_+ blindTime_);

	      double wfsumMax(0);
	      size_t peakP(0);
	      for (int i=sampleStart; i<sampleStop-nBinsPeak_; ++i)
              {
                  float sum = std::accumulate(&wf[i],&wf[i+nBinsPeak_],0.0);
                  if (sum>wfsumMax){wfsumMax=sum; peakP = i-sampleStart+nBinsPeak_/2;}
	      }
              
              std::vector<int> wfsample;
              for (int i=sampleStart; i<=sampleStop; ++i) wfsample.push_back(std::min(int(wf[i]),maxADCCounts_));

              // make the CaloDigi
              caloDigiColl.emplace_back(CaloDigi(iRO,t0,wfsample,peakP) );

              if (diagLevel_ > 2) diag1(iRO,t0, peakP, wfsample);
              if (iRO%2==0) totEdepEq += float(wfsample[peakP])/MeVToADC_;
          }
      }
      if (diagLevel_ >0) std::cout<<"[CaloDigiFromShower] Total energy equivalent digitized "<<totEdepEq<<std::endl;
  }




  void CaloDigiFromShower::diag0(int iRO, const std::vector<double>& wf )
  {
     if (*std::max_element(wf.begin(),wf.end())<1) return;
     std::cout<<"CaloDigiFromShower::fillOutoutRO] Waveform content for readout "<<iRO<<std::endl;
     for (unsigned i=0;i<wf.size();++i) {if (i%10==0 && i>0) std::cout<<"- "; std::cout<<int(wf[i])<<" ";}
     std::cout<<std::endl;
  }

  void CaloDigiFromShower::diag1(int iRO, double time, size_t peakP, const std::vector<int>& wf)
  {
     std::cout<<"Created caloDigi with roID = "<<iRO<<"  time="<<time<<" peak="<<peakP<<"  and content ";
     for (const auto  &v : wf) std::cout<<v<<" ";
     std::cout<<std::endl;
  }

  void CaloDigiFromShower::diag2(const CaloDigiCollection& caloDigiColl)
  {      
     std::map<int,double> enerMap;
     for (const auto& digi : caloDigiColl) enerMap[digi.roId()] += digi.waveform().at(digi.peakpos())/MeVToADC_;
     std::cout<<"[CaloDigiFromShower] energy equivalent per RoID"<<std::endl;
     for (auto& kv : enerMap) std::cout<<" roID: "<<kv.first<<"   Ener: "<<kv.second<<std::endl;
  }






}

DEFINE_ART_MODULE(mu2e::CaloDigiFromShower);
