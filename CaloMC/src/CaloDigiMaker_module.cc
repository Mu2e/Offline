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

#include "Offline/CaloConditions/inc/CalSimParams.hh"
#include "Offline/Mu2eUtilities/inc/CaloPulseUtil.hh"
#include "Offline/Mu2eUtilities/inc/CaloNoiseUtil.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/MCDataProducts/inc/CaloShowerRO.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"

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
             using CNG     = CaloNoiseUtil::Config;
             using CPG     = CaloPulseUtil::Config;
             fhicl::Table<CNG>          noiseCache           { Name("noiseCache"),             Comment("Noise cache maker config") };
             fhicl::Table<CPG>          pulseCache           { Name("pulseCache"),             Comment("Pulse cache maker config") };
             fhicl::Atom<art::InputTag> caloShowerCollection { Name("caloShowerROCollection"), Comment("CaloShowerRO collection name") };
             fhicl::Atom<art::InputTag> ewMarkerTag          { Name("eventWindowMarker"),      Comment("EventWindowMarker producer") };
             fhicl::Atom<art::InputTag> pbtmcTag             { Name("protonBunchTimeMC"),      Comment("ProtonBunchTimeMC producer") };
             fhicl::Atom<float>         digitizationStart    { Name("digitizationStart"),      Comment("Start of digitization window relative to nominal pb time") };
             fhicl::Atom<float>         digitizationEnd      { Name("digitizationEnd"),        Comment("End of digitization window relative to nominal pb time")};
             fhicl::Atom<bool>          addNoise             { Name("addNoise"),               Comment("Add noise to waveform") };
             fhicl::Atom<float>         digiSampling         { Name("digiSampling"),           Comment("Digitization time sampling") };
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
            maxADCCounts_      ((1 << config().nBits()) -1),
            pulseCache_        (CaloPulseUtil(config().pulseCache())),
            nBinsPeak_         (config().nBinsPeak()),
            minPeakADC_        (config().minPeakADC()),
            engine_            (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
            addNoise_          (config().addNoise()),
            noiseSampler_      (config().noiseCache(), engine_),
            diagLevel_         (config().diagLevel())
         {
             consumes<EventWindowMarker>(ewMarkerTag_);
             consumes<ProtonBunchTimeMC>(pbtmcTag_);

             produces<CaloDigiCollection>();

             //check that StartTimeBuffer is shorter than BlindTime_
             if (startTimeBuffer_ > digitizationStart_) throw cet::exception("CALODIGIMAKER")
               << "digitizationStart too small to accommodate start time buffer";
          }

         void produce(art::Event& e)   override;
         void beginRun(art::Run& aRun) override;


    private:

       void   makeDigitization  (const CaloShowerROCollection&, CaloDigiCollection&, const EventWindowMarker&, const ProtonBunchTimeMC&, const CalSimParams&);
       bool   fillROHits        (unsigned iRO, std::vector<double>& waveform, const CaloShowerROCollection&, const ProtonBunchTimeMC&, const CalSimParams&);
       void   generateSpotNoise (std::vector<double>& waveform, double scaleFactor);
       void   buildOutputDigi   (unsigned iRO, std::vector<double>& waveform, double pedestal, CaloDigiCollection&);
       void   extract           (const std::vector<int>& wf, std::vector<size_t>& starts, std::vector<size_t>& stops) const;
       double readoutScaleFactor(unsigned iRO, const CalSimParams& conds) const;
       void   diag0             (unsigned, const std::vector<int>&);
       void   diag1             (unsigned, double, size_t, const std::vector<int>&, int);

       ProditionsHandle<CalSimParams>                  calCrystalConds_;
       const art::ProductToken<CaloShowerROCollection> caloShowerToken_;
       art::InputTag           ewMarkerTag_;
       art::InputTag           pbtmcTag_;
       float                   digitizationStart_;
       float                   digitizationEnd_;
       float                   timeFromProtonsToDRMarker_;
       float                   digiSampling_;
       unsigned                bufferDigi_;
       float                   startTimeBuffer_;
       int                     maxADCCounts_;
       CaloPulseUtil           pulseCache_;
       unsigned                nBinsPeak_;
       int                     minPeakADC_;
       CLHEP::HepRandomEngine& engine_;
       bool                    addNoise_;
       CaloNoiseUtil           noiseSampler_;
       const Calorimeter*      calorimeter_;
       int                     diagLevel_;
  };


  //-----------------------------------------------------------------------------
  void CaloDigiMaker::beginRun(art::Run& aRun)
  {
      pulseCache_.buildCache();
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

      const auto& calCrystalConds = calCrystalConds_.get(event.id());

      auto caloDigiColl = std::make_unique<CaloDigiCollection>();
      makeDigitization(CaloShowerROs, *caloDigiColl,ewMarker, pbtmc,calCrystalConds);
      event.put(std::move(caloDigiColl));

      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiMaker::produce] end" << std::endl;
  }


  //-----------------------------------------------------------------------------------------------------------------------------
  // Note: DigitizationStart include the fixed delay from timeFromProtonsToDRMarker, need to subtract it to be in the digitizer frame
  void CaloDigiMaker::makeDigitization(const CaloShowerROCollection& CaloShowerROs, CaloDigiCollection& caloDigiColl,
                                       const EventWindowMarker& ewMarker, const ProtonBunchTimeMC& pbtmc, const CalSimParams& calCrystalConds)
  {
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      calorimeter_ = ch.get();

      if (calorimeter_->nCrystals()<1 || calorimeter_->G4Info().get<int>("nSiPMPerCrystal")<1) return;
      int waveformSize = (digitizationEnd_ - digitizationStart_ + startTimeBuffer_) / digiSampling_;
      if (ewMarker.spillType() != EventWindowMarker::SpillType::onspill)
      {
        waveformSize = (ewMarker.eventLength() - digitizationStart_ + startTimeBuffer_) / digiSampling_;
      }

      int nWaveforms   = calorimeter_->nCrystals()*calorimeter_->G4Info().get<int>("nSiPMPerCrystal");
      if (waveformSize<1) throw cet::exception("Rethrow")<< "[CaloMC/CaloDigiMaker] digitization size too short " << std::endl;
      bool resetWaveform(false);
      std::vector<double> waveform(waveformSize,0.0);

      for (int iRO=0;iRO<nWaveforms;++iRO)
      {
          if (resetWaveform) std::fill(waveform.begin(), waveform.end(), 0.0);
          bool isEmpty = fillROHits(iRO, waveform, CaloShowerROs, pbtmc, calCrystalConds );
          resetWaveform = !isEmpty;

          if (isEmpty) continue;
          if (addNoise_) {
            const double scaleFactor = readoutScaleFactor(iRO, calCrystalConds);
            generateSpotNoise(waveform,scaleFactor);
          }
          buildOutputDigi(iRO, waveform, noiseSampler_.pedestal(), caloDigiColl);
     }
  }

  //----------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::generateSpotNoise(std::vector<double>& waveform, double scaleFactor)
  {
       const int NoiseWFID(0); // will need this from proditions later;

       const double minAmplitude(0.01);
       noiseSampler_.prepare(NoiseWFID, scaleFactor);

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
          unsigned istart  = hitStarts[ihit];
          unsigned ilength = hitStops[ihit]-hitStarts[ihit];
          const auto& noiseWF = noiseSampler_.noiseSegment(NoiseWFID,istart,ilength);
          for (unsigned i=0;i<ilength;++i) waveform[istart+i] += noiseWF[i];
       }
  }

  //--------------------------------------------------------------------------
  bool CaloDigiMaker::fillROHits(unsigned iRO, std::vector<double>& waveform, const CaloShowerROCollection& CaloShowerROs,
                                 const ProtonBunchTimeMC& pbtmc, const CalSimParams& calCrystalConds)
  {
      bool isEmpty = true;
      const double scaleFactor = readoutScaleFactor(iRO, calCrystalConds);

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
              const auto& pulse          = pulseCache_.digitizedPulse(time);
              unsigned    stopSample     = std::min(startSample+pulse.size(), waveform.size());

              for (size_t timeSample = startSample; timeSample < stopSample; ++timeSample)
                 waveform.at(timeSample) += pulse.at(timeSample - startSample)*scaleFactor;
          }
      }
      return isEmpty;
  }


  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::buildOutputDigi(unsigned iRO, std::vector<double>& waveform, double pedestal, CaloDigiCollection& caloDigiColl)
  {
       // round the waveform into non-null integers and apply maxADC cut
       std::vector<int> wf(waveform.size(),0);
       for (size_t i=0; i<waveform.size(); ++i) {
         if (waveform[i] > pedestal) wf[i] = std::min(maxADCCounts_, int(waveform[i] - pedestal));
       }
       if (diagLevel_ > 2) diag0(iRO, wf);

       //extract hits start / stop times
       std::vector<size_t> hitStarts, hitStops;
       hitStarts.reserve(16);hitStops.reserve(16);
       extract(wf,hitStarts,hitStops);

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
  void CaloDigiMaker::extract(const std::vector<int>& wf, std::vector<size_t>& starts, std::vector<size_t>& stops) const
  {
       size_t timeSample(nBinsPeak_+bufferDigi_);
       while (timeSample+nBinsPeak_ < wf.size())
       {
           // find starting point
           if (wf[timeSample] < minPeakADC_) {++timeSample; continue;}

           size_t imax(timeSample-nBinsPeak_);
           for (auto i = timeSample-nBinsPeak_; i<=timeSample+nBinsPeak_;++i) {if (wf[i]>wf[imax]) imax=i;}
           if (timeSample != imax)  {++timeSample; continue;}

           // find the starting / stopping point of the peak (stop = first value under threshold)
           size_t sampleStart = (timeSample > bufferDigi_) ? timeSample - bufferDigi_ : 0;
           size_t sampleStop(timeSample);
           ++sampleStop;
           while (sampleStop < wf.size() && wf[sampleStop] >= minPeakADC_) ++sampleStop;

           starts.push_back(sampleStart);
           stops.push_back(sampleStop);

           //fast forward to end of waveform to search for next one
           timeSample = sampleStop+1;
       }


       // Concatenate peaks and remove unused values (flag value to remove past wf.size() since the latter is a legitimate value)
       size_t iprev(0), icurrent(1);
       while (icurrent < starts.size())
       {
           if (stops[iprev] >= starts[icurrent]) {stops[iprev]=stops[icurrent]; starts[icurrent]=stops[icurrent]=wf.size()+1;}
           else                                  {iprev = icurrent;}
           ++icurrent;
       }

       auto pred = [&wf](const auto a) {return a>wf.size();};
       starts.erase(std::remove_if(starts.begin(),starts.end(),pred),starts.end());
       stops.erase(std::remove_if(stops.begin(), stops.end(), pred),stops.end());
  }

  //-------------------------------------------------------------------------------------------------------------------
  double CaloDigiMaker::readoutScaleFactor(unsigned iRO, const CalSimParams& conds) const
  {
      const auto SiPMID    = CaloSiPMId(iRO);
      const auto crystalID = SiPMID.crystal();
      const auto pePerMeV  = conds.pePerMeVs (crystalID).at(SiPMID.SiPMLocalId());
      const auto ADCPerMeV = conds.ADCPerMeVs(crystalID).at(SiPMID.SiPMLocalId());
      return ADCPerMeV / pePerMeV;
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

}

DEFINE_ART_MODULE(mu2e::CaloDigiMaker)
