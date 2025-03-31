//
// Simulate the readout waveform for each sensors from CaloShowerROs.
// Individual photo-electrons are generated for each readout, including photo-statistic fluctuations
// Simulate digitization procedure and produce CaloDigis.
//
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "Offline/Mu2eUtilities/inc/CaloPulseShape.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloMC/inc/CaloNoiseSimGenerator.hh"
#include "Offline/CaloMC/inc/CaloWFExtractor.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/ConditionsService/inc/CalorimeterCalibrations.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/MCDataProducts/inc/CaloShowerRO.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"

#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TGraph.h"

#include <iostream>
#include <string>
#include <vector>
#include <numeric>


namespace mu2e {


  class CaloDigiMaker : public art::EDProducer
  {
     public:
         struct Config
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             using CNG     = mu2e::CaloNoiseSimGenerator::Config;
             fhicl::Table<CNG>          noise_gen_conf       { Name("NoiseGenerator"),         Comment("Noise generator config") };
             fhicl::Atom<art::InputTag> caloShowerCollection { Name("caloShowerROCollection"), Comment("CaloShowerRO collection name") };
             fhicl::Atom<art::InputTag> ewMarkerTag          { Name("eventWindowMarker"),      Comment("EventWindowMarker producer") };
             fhicl::Atom<art::InputTag> pbtmcTag             { Name("protonBunchTimeMC"),      Comment("ProtonBunchTimeMC producer") };
             fhicl::Atom<double>        digitizationStart    { Name("digitizationStart"),      Comment("Start of digitization window relative to nominal pb time") };
             fhicl::Atom<double>        digitizationEnd      { Name("digitizationEnd"),        Comment("End of digitization window relative to nominal pb time")};
             fhicl::Atom<bool>          addNoise             { Name("addNoise"),               Comment("Add noise to waveform") };
             fhicl::Atom<bool>          addRandomNoise       { Name("addRandomNoise"),         Comment("Add random salt and pepper noise") };
             fhicl::Atom<double>        digiSampling         { Name("digiSampling"),           Comment("Digitization time sampling") };
             fhicl::Atom<int>           nBits                { Name("nBits"),                  Comment("ADC Number of bits") };
             fhicl::Atom<unsigned>      nBinsPeak            { Name("nBinsPeak"),              Comment("Window size for finding local maximum to digitize wf") };
             fhicl::Atom<int>           minPeakADC           { Name("minPeakADC"),             Comment("Minimum ADC hits of local peak to digitize") };
             fhicl::Atom<unsigned>      bufferDigi           { Name("bufferDigi"),             Comment("Number of timeStamps for the buffer digi") };
             fhicl::Atom<int>           diagLevel            { Name("diagLevel"),              Comment("Diag Level"),0 };
         };

         explicit CaloDigiMaker(const art::EDProducer::Table<Config>& config) :
            EDProducer{config},
            caloShowerToken_{consumes<CaloShowerROCollection>(config().caloShowerCollection())},
            ewMarkerTag_       (config().ewMarkerTag()),
            pbtmcTag_          (config().pbtmcTag()),
            digitizationStart_ (config().digitizationStart()),
            digitizationEnd_   (config().digitizationEnd()),
            digiSampling_      (config().digiSampling()),
            bufferDigi_        (config().bufferDigi()),
            startTimeBuffer_   (config().digiSampling()*config().bufferDigi()),
            maxADCCounts_      (1 << config().nBits()),
            pulseShape_        (CaloPulseShape(config().digiSampling())),
            wfExtractor_       (config().bufferDigi(),config().nBinsPeak(),config().minPeakADC(),config().bufferDigi()),
            engine_            (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
            addNoise_          (config().addNoise()),
            noiseGenerator_    (config().noise_gen_conf(), engine_, 0),
            addRandomNoise_    (config().addRandomNoise()),
            diagLevel_         (config().diagLevel())
         {
             consumes<EventWindowMarker>(ewMarkerTag_);
             consumes<ProtonBunchTimeMC>(pbtmcTag_);

             produces<CaloDigiCollection>();

             //check that StartTimeBuffer is shorter than BlindTime_
             if (startTimeBuffer_ > digitizationStart_) throw cet::exception("CATEGORY")<< "CaloDigiMaker: digitizationStart is too small to accommodate start time buffer";
          }

         void produce(art::Event& e)   override;
         void beginRun(art::Run& aRun) override;


    private:

       void makeDigitization  (const CaloShowerROCollection&, CaloDigiCollection&, const EventWindowMarker&, const ProtonBunchTimeMC&);
       bool fillROHits        (unsigned iRO, std::vector<double>& waveform, const CaloShowerROCollection&,
                               const ConditionsHandle<CalorimeterCalibrations>&, const ProtonBunchTimeMC&);
       void generateSpotNoise (std::vector<double>& waveform, unsigned iRO, const ConditionsHandle<CalorimeterCalibrations>&);
       void buildOutputDigi   (unsigned iRO, std::vector<double>& waveform, int pedestal, CaloDigiCollection&);
       void diag0             (unsigned, const std::vector<int>&);
       void diag1             (unsigned, double, size_t, const std::vector<int>&, int);
       void plotWF            (const std::vector<int>& waveform,    const std::string& pname, int pedestal);
       void plotWF            (const std::vector<double>& waveform, const std::string& pname, int pedestal);

       const art::ProductToken<CaloShowerROCollection> caloShowerToken_;
       art::InputTag           ewMarkerTag_;
       art::InputTag           pbtmcTag_;
       double                  digitizationStart_;
       double                  digitizationEnd_;
       double                  timeFromProtonsToDRMarker_;
       double                  digiSampling_;
       unsigned                bufferDigi_;
       double                  startTimeBuffer_;
       int                     maxADCCounts_;
       CaloPulseShape          pulseShape_;
       CaloWFExtractor         wfExtractor_;
       CLHEP::HepRandomEngine& engine_;
       bool                    addNoise_;
       CaloNoiseSimGenerator   noiseGenerator_;
       bool                    addRandomNoise_;
       const Calorimeter*      calorimeter_;
       int                     diagLevel_;
  };


  //-----------------------------------------------------------------------------
  void CaloDigiMaker::beginRun(art::Run& aRun)
  {
      pulseShape_.buildShapes();
      if (addNoise_) noiseGenerator_.initialize(wfExtractor_);
  }



  //---------------------------------------------------------
  void CaloDigiMaker::produce(art::Event& event)
  {

      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiMaker::produce] begin" << std::endl;

      //get Event window and bunch timing info
      art::Handle<EventWindowMarker> ewMarkerHandle;
      event.getByLabel(ewMarkerTag_, ewMarkerHandle);
      const EventWindowMarker& ewMarker(*ewMarkerHandle);

      art::Handle<ProtonBunchTimeMC> pbtmcHandle;
      event.getByLabel(pbtmcTag_, pbtmcHandle);
      const ProtonBunchTimeMC& pbtmc(*pbtmcHandle);

      ProditionsHandle<EventTiming> eventTimingHandle;
      const EventTiming &eventTiming = eventTimingHandle.get(event.id());
      timeFromProtonsToDRMarker_ = eventTiming.timeFromProtonsToDRMarker();

      auto caloShowerStepHandle = event.getValidHandle(caloShowerToken_);
      const auto& CaloShowerROs = *caloShowerStepHandle;

      auto caloDigiColl = std::make_unique<CaloDigiCollection>();
      makeDigitization(CaloShowerROs, *caloDigiColl,ewMarker, pbtmc);
      event.put(std::move(caloDigiColl));

      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiMaker::produce] end" << std::endl;
  }


  //-----------------------------------------------------------------------------------------------------------------------------
  // Note: DigitizationStart include the fixed delay from timeFromProtonsToDRMarker, need to subtract it to be in the digitizer frame
  void CaloDigiMaker::makeDigitization(const CaloShowerROCollection& CaloShowerROs, CaloDigiCollection& caloDigiColl,
                                       const EventWindowMarker& ewMarker, const ProtonBunchTimeMC& pbtmc)
  {
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      calorimeter_ = ch.get();

      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

      if (calorimeter_->nCrystals()<1 || calorimeter_->caloInfo().getInt("nSiPMPerCrystal")<1) return;
      int waveformSize = (digitizationEnd_ - digitizationStart_ + startTimeBuffer_) / digiSampling_;
      if (ewMarker.spillType() != EventWindowMarker::SpillType::onspill)
      {
        waveformSize = (ewMarker.eventLength() - digitizationStart_ + startTimeBuffer_) / digiSampling_;
      }

      int nWaveforms   = calorimeter_->nCrystals()*calorimeter_->caloInfo().getInt("nSiPMPerCrystal");
      if (waveformSize<1) throw cet::exception("Rethrow")<< "[CaloMC/CaloDigiMaker] digitization size too short " << std::endl;
      bool resetWaveform(false);
      std::vector<double> waveform(waveformSize,0.0);

      for (int iRO=0;iRO<nWaveforms;++iRO)
      {
          if (resetWaveform) std::fill(waveform.begin(), waveform.end(), 0.0);
          bool isEmpty = fillROHits(iRO, waveform, CaloShowerROs, calorimeterCalibrations, pbtmc);
          resetWaveform = (addRandomNoise_ || !isEmpty);

          // if we add random noise, then we need to scan all waveforms. Otherwise we can skip empty waveforms
          if (addRandomNoise_) {
            if (!isEmpty) generateSpotNoise(waveform, iRO, calorimeterCalibrations);
            noiseGenerator_.addSaltAndPepper(waveform);
            buildOutputDigi(iRO, waveform, noiseGenerator_.pedestal(), caloDigiColl);
          }
          else
          {
            if (isEmpty) continue;
            if (addNoise_) generateSpotNoise(waveform, iRO, calorimeterCalibrations);
            buildOutputDigi(iRO, waveform, noiseGenerator_.pedestal(), caloDigiColl);
          }
     }
  }


  //--------------------------------------------------------------------------
  bool CaloDigiMaker::fillROHits(unsigned iRO, std::vector<double>& waveform, const CaloShowerROCollection& CaloShowerROs,
                                 const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations, const ProtonBunchTimeMC& pbtmc)
  {
      bool  isEmpty(true);
      float scaleFactor = calorimeterCalibrations->MeV2ADC(iRO)/calorimeterCalibrations->peMeV(iRO);

      for (const auto& CaloShowerRO : CaloShowerROs)
      {
          unsigned SiPMID = CaloShowerRO.SiPMID();
          if (SiPMID != iRO) continue;

          isEmpty = false;
          for (const auto PEtime : CaloShowerRO.PETime())
          {
              //PE time is given in DR frame, we need to subtract the event window start and the digi Start time
              float       time           = PEtime + pbtmc.pbtime_- digitizationStart_ + timeFromProtonsToDRMarker_ + startTimeBuffer_;
              unsigned    startSample    = std::max(0u,unsigned(time/digiSampling_));
              const auto& pulse          = pulseShape_.digitizedPulse(time);
              unsigned    stopSample     = std::min(startSample+pulse.size(), waveform.size());

              for (size_t timeSample = startSample; timeSample < stopSample; ++timeSample)
                 waveform.at(timeSample) += pulse.at(timeSample - startSample)*scaleFactor;
          }
      }
      return isEmpty;
  }


  //----------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::generateSpotNoise(std::vector<double>& waveform, unsigned iRO,
                                        const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
       double minAmplitude = 0.1*calorimeterCalibrations->MeV2ADC(iRO);

       size_t timeSample(0);
       std::vector<size_t> hitStarts{}, hitStops{};
       hitStarts.reserve(16);hitStops.reserve(16);

       // First, find the ranges in the waveform with non-zero bins.
       while (timeSample < waveform.size())
       {
           if (waveform[timeSample] < minAmplitude) {++timeSample; continue;}

           size_t sampleStart = (timeSample > bufferDigi_) ? timeSample - bufferDigi_ : 0;
           size_t sampleStop(timeSample);
           while (sampleStop < waveform.size() && waveform[sampleStop] > minAmplitude) ++sampleStop;

           hitStarts.push_back(sampleStart);
           hitStops.push_back(sampleStop);
           timeSample = sampleStop+1;
       }

       // ranges might overlap and need to be concatenated if this is the case
       size_t iprev(0),ic(1);
       while (ic < hitStarts.size())
       {
           if (hitStops[iprev] >= hitStarts[ic]) {hitStops[iprev]=hitStops[ic]; hitStarts[ic]=hitStops[ic]=waveform.size()+1;}
           else                                  {iprev = ic;}
           ++ic;
       }
       auto pred = [&waveform](const auto a) {return a>waveform.size();};
       hitStarts.erase(std::remove_if(hitStarts.begin(),hitStarts.end(),pred),hitStarts.end());
       hitStops.erase( std::remove_if(hitStops.begin(), hitStops.end(), pred),hitStops.end());

       //Now take a random part of the noise waveform and add it to the waveform content
       for (size_t ihit=0; ihit<hitStarts.size(); ++ihit)
       {
          noiseGenerator_.addSampleNoise(waveform,hitStarts[ihit],hitStops[ihit]-hitStarts[ihit]);
       }
  }


  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::buildOutputDigi(unsigned iRO, std::vector<double>& waveform, int pedestal, CaloDigiCollection& caloDigiColl)
  {
       // round the waveform into non-null integers and apply maxADC cut
       std::vector<int> wf(waveform.size(),0);
       for (unsigned i=0; i<waveform.size(); ++i) {if (waveform[i] > pedestal) wf[i] = std::min(maxADCCounts_, int(waveform[i] - pedestal));}
       if (diagLevel_ > 2) diag0(iRO, wf);

       //extract hits start / stop times
       std::vector<size_t> hitStarts, hitStops;
       hitStarts.reserve(16);hitStops.reserve(16);
       wfExtractor_.extract(wf,hitStarts,hitStops);

       // Build digi for concatenated hits
       for (size_t ihit=0;ihit<hitStarts.size();++ihit)
       {
           size_t sampleStart = hitStarts[ihit];
           size_t sampleStop  = hitStops[ihit];
           size_t t0          = size_t(sampleStart*digiSampling_ + digitizationStart_ - timeFromProtonsToDRMarker_ - startTimeBuffer_);
           //t0 is given in the "digitizer time frame"

           std::vector<int> wfsample{};
           wfsample.reserve(sampleStop-sampleStart);
           for (size_t i=sampleStart; i<sampleStop; ++i) wfsample.push_back(std::min(int(waveform[i]),maxADCCounts_));

           // only consider hits above digitizationStart
           size_t peakPosition(0u);
           for (auto i = 0u; i<wfsample.size();++i) {
              if (t0+i*digiSampling_+timeFromProtonsToDRMarker_ >= digitizationStart_ && wfsample[i]>=wfsample[peakPosition]) peakPosition=i;
           }
           if (diagLevel_ >2) std::cout<<"[CaloDigiMaker] Start=" << sampleStart << " Stop=" << sampleStop
                                       << " peak in position " << peakPosition << std::endl;

           caloDigiColl.emplace_back(CaloDigi(iRO,t0,wfsample,peakPosition) );
           if (diagLevel_ > 2) diag1(iRO, t0, peakPosition, wfsample, pedestal);
       }
  }





  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::diag0(unsigned iSiPM, const std::vector<int>& wf)
  {
      if (*std::max_element(wf.begin(),wf.end())<1) return;
      std::cout<<"CaloDigiMaker::fillOutoutRO] Waveform content for readout "<<iSiPM<<std::endl;
      for (size_t i=0;i<wf.size();++i) {if (i%10==0 && i>0) std::cout<<"- "; std::cout<<wf[i]<<" ";}
      std::cout<<std::endl;
  }
  void CaloDigiMaker::diag1(unsigned iSiPM, double time, size_t peakP, const std::vector<int>& wf, int ped)
  {
      std::cout<<"Created caloDigi with SiPM = "<<iSiPM<<"  t0="<<time<<" peak="<<peakP<<" and content ";
      for (const auto  &v : wf) std::cout<<v-ped<<" ";
      std::cout<<std::endl;
  }




  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::plotWF(const std::vector<int>& waveform, const std::string& pname, int pedestal)
  {
      double startTime = digitizationStart_ - timeFromProtonsToDRMarker_ - startTimeBuffer_;
      TH1F h("h","Waveform",waveform.size(),startTime,waveform.size()*digiSampling_+startTime);
      for (size_t i=1;i<=waveform.size();++i) h.SetBinContent(i,waveform[i-1]);
      TLine line;
      line.SetLineStyle(2);
      TLine line2;
      line2.SetLineStyle(3);

      gStyle->SetOptStat(0);
      TCanvas c1("c1","c1");
      h.Draw();
      line.DrawLine(startTime,pedestal,waveform.size()*digiSampling_+startTime,pedestal);
      line2.DrawLine(startTime,pedestal+16,waveform.size()*digiSampling_+startTime,pedestal+16);
      c1.SaveAs(pname.c_str());
  }

  void CaloDigiMaker::plotWF(const std::vector<double>& waveform, const std::string& pname, int pedestal)
  {
      std::vector<int> v;
      for (size_t i=0;i<waveform.size();++i) v.push_back(int(waveform[i]));
      plotWF(v,pname, pedestal);
  }


}

DEFINE_ART_MODULE(mu2e::CaloDigiMaker)
