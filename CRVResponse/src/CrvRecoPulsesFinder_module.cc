//
// A module to extract leading edge times, pulse heights, integrals and number of photons from the CRV waveforms
//
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvRecoPulses.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/CrvParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulse.hh"
#include "RecoDataProducts/inc/CrvRecoPulseFlags.hh"
#include "RecoDataProducts/inc/ProtonBunchTime.hh"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>

#include <TMath.h>

namespace mu2e 
{
  class CrvRecoPulsesFinder : public art::EDProducer 
  {

    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> crvDigiModuleLabel{Name("crvDigiModuleLabel"), Comment("module label for CrvDigis")};
      fhicl::Atom<float> minADCdifference{Name("minADCdifference"), Comment("minimum ADC difference above pedestal to be considered for reconstruction")};  //5.0
      fhicl::Atom<float> defaultBeta{Name("defaultBeta"), Comment("initialization value for fit and default value for invalid fits (regular pulses: 19.0ns, dark counts for calibration: 12.0ns)")};
      fhicl::Atom<float> minBeta{Name("minBeta"), Comment("smallest accepted beta for valid fit [ns]")};  //5.0ns
      fhicl::Atom<float> maxBeta{Name("maxBeta"), Comment("largest accepted beta for valid fit [ns]")}; //50.0ns
      fhicl::Atom<float> maxTimeDifference{Name("maxTimeDifference"), Comment("largest accepted difference between time of largest ADC value and fitted peak [ns]")}; //20.0ns
      fhicl::Atom<float> minPulseHeightRatio{Name("minPulseHeightRatio"), Comment("smallest accepted ratio between largest ADC value and fitted peak")}; //0.7
      fhicl::Atom<float> maxPulseHeightRatio{Name("maxPulseHeightRatio"), Comment("largest accepted ratio between largest ADC value and fitted peak")}; //1.5
      fhicl::Atom<float> LEtimeFactor{Name("LEtimeFactor"), Comment("time of leading edge is peakTime-LEtimeFactor*beta (0.985,1.385,1.587 for a leading edge of 0.5,0.2,0.1 pulse height")};
      fhicl::Atom<bool> allowDoubleGumbel{Name("allowDoubleGumbel"), Comment("tries fitting with two Gumbel functions")};
      fhicl::Atom<float> doubleGumbelThreshold{Name("doubleGumbelThreshold"), Comment("Chi2/#ADCsamples (based on single Gumbel fit) at which a fit with two Gumbel functions should be attempted")};
      fhicl::Atom<float> minPEs{Name("minPEs"), Comment("minimum number of PEs")}; //0
      fhicl::Atom<art::InputTag> protonBunchTimeTag{ Name("protonBunchTimeTag"), Comment("ProtonBunchTime producer"),"EWMProducer" };
    };

    typedef art::EDProducer::Table<Config> Parameters;

    explicit CrvRecoPulsesFinder(const Parameters& config);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    boost::shared_ptr<mu2eCrv::MakeCrvRecoPulses> _makeCrvRecoPulses;

    std::string _crvDigiModuleLabel;
    float       _minPEs;
    float       _digitizationPeriod;
    float       _pedestal;           //100 ADC
    float       _calibrationFactor;  //394.6 ADC*ns/PE
    float       _calibrationFactorPulseHeight;  //11.4 ADC/PE
    art::InputTag _protonBunchTimeTag;
  };


  CrvRecoPulsesFinder::CrvRecoPulsesFinder(const Parameters& conf) :
    art::EDProducer(conf),
    _crvDigiModuleLabel(conf().crvDigiModuleLabel()),
    _minPEs(conf().minPEs()),
    _protonBunchTimeTag(conf().protonBunchTimeTag())
  {
    produces<CrvRecoPulseCollection>();
    _makeCrvRecoPulses=boost::shared_ptr<mu2eCrv::MakeCrvRecoPulses>(new mu2eCrv::MakeCrvRecoPulses(conf().minADCdifference(),
                                                                                                    conf().defaultBeta(),
                                                                                                    conf().minBeta(),
                                                                                                    conf().maxBeta(),
                                                                                                    conf().maxTimeDifference(),
                                                                                                    conf().minPulseHeightRatio(),
                                                                                                    conf().maxPulseHeightRatio(),
                                                                                                    conf().LEtimeFactor(),
                                                                                                    conf().allowDoubleGumbel(),
                                                                                                    conf().doubleGumbelThreshold()));
  }

  void CrvRecoPulsesFinder::beginJob()
  {
  }

  void CrvRecoPulsesFinder::endJob()
  {
  }

  void CrvRecoPulsesFinder::beginRun(art::Run &run)
  {
    mu2e::ConditionsHandle<mu2e::CrvParams> crvPar("ignored");
    _digitizationPeriod = crvPar->digitizationPeriod;
    _pedestal           = crvPar->pedestal;
    _calibrationFactor  = crvPar->calibrationFactor;
    _calibrationFactorPulseHeight  = crvPar->calibrationFactorPulseHeight;
  }

  void CrvRecoPulsesFinder::produce(art::Event& event) 
  {
    std::unique_ptr<CrvRecoPulseCollection> crvRecoPulseCollection(new CrvRecoPulseCollection);

    art::Handle<ProtonBunchTime> protonBunchTime;
    event.getByLabel(_protonBunchTimeTag, protonBunchTime);
    double TDC0time = -protonBunchTime->pbtime_;

    art::Handle<CrvDigiCollection> crvDigiCollection;
    event.getByLabel(_crvDigiModuleLabel,"",crvDigiCollection);

    size_t waveformIndex = 0;
    while(waveformIndex<crvDigiCollection->size())
    {
      const CrvDigi &digi = crvDigiCollection->at(waveformIndex);
      const CRSScintillatorBarIndex &barIndex = digi.GetScintillatorBarIndex();
      int SiPM = digi.GetSiPMNumber();
      unsigned int startTDC = digi.GetStartTDC();
      std::vector<unsigned int> ADCs;
      std::vector<size_t> waveformIndices;
      for(size_t i=0; i<CrvDigi::NSamples; ++i) ADCs.push_back(digi.GetADCs()[i]);
      waveformIndices.push_back(waveformIndex);

      //checking following digis whether they are a continuation of the current digis
      //if that is the case, append the next digis
      while(++waveformIndex<crvDigiCollection->size())
      {
        const CrvDigi &nextDigi = crvDigiCollection->at(waveformIndex);
        if(barIndex!=nextDigi.GetScintillatorBarIndex()) break;
        if(SiPM!=nextDigi.GetSiPMNumber()) break;
        if(startTDC+ADCs.size()!=nextDigi.GetStartTDC()) break;
        for(size_t i=0; i<CrvDigi::NSamples; ++i) ADCs.push_back(nextDigi.GetADCs()[i]);
        waveformIndices.push_back(waveformIndex);
      }

      _makeCrvRecoPulses->SetWaveform(ADCs, startTDC, _digitizationPeriod, _pedestal, _calibrationFactor, _calibrationFactorPulseHeight);

      size_t n = _makeCrvRecoPulses->GetPEs().size();
      for(size_t j=0; j<n; ++j)
      {
        //the TDC times were recorded with respect to the event window start.
        //need to shift the times back to the original time scale (i.e. microbunch time)
        double pulseTime   = _makeCrvRecoPulses->GetPulseTimes().at(j) + TDC0time;
        double LEtime      = _makeCrvRecoPulses->GetLEtimes().at(j) + TDC0time;
        float  PEs         = _makeCrvRecoPulses->GetPEs().at(j);
        float  PEsPulseHeight = _makeCrvRecoPulses->GetPEsPulseHeight().at(j);
        float  pulseHeight = _makeCrvRecoPulses->GetPulseHeights().at(j); 
        float  pulseBeta   = _makeCrvRecoPulses->GetPulseBetas().at(j);
        float  pulseFitChi2= _makeCrvRecoPulses->GetPulseFitChi2s().at(j);

        bool   failedFit              = _makeCrvRecoPulses->GetFailedFits().at(j);
        CrvRecoPulseFlags flags;
        if(failedFit)              flags.set(CrvRecoPulseFlagEnums::failedFit);

        float  PEsNoFit          = _makeCrvRecoPulses->GetPEsNoFit().at(j);
        double pulseTimeNoFit    = _makeCrvRecoPulses->GetPulseTimesNoFit().at(j) + TDC0time;
        double pulseStart        = _makeCrvRecoPulses->GetPulseStarts().at(j) + TDC0time;
        double pulseEnd          = _makeCrvRecoPulses->GetPulseEnds().at(j) + TDC0time;

        if(PEs<_minPEs) continue; 
        crvRecoPulseCollection->emplace_back(PEs, PEsPulseHeight, pulseTime, pulseHeight, pulseBeta, pulseFitChi2, LEtime, flags, 
                                             PEsNoFit, pulseTimeNoFit, pulseStart, pulseEnd,
                                             waveformIndices, barIndex, SiPM);
      }

    }

    event.put(std::move(crvRecoPulseCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvRecoPulsesFinder;
DEFINE_ART_MODULE(CrvRecoPulsesFinder)
