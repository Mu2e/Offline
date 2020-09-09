//
// A module to extract leading edge times, pulse heights, integrals and number of photons from the CRV waveforms
//
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvRecoPulses.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CrvParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulseCollection.hh"

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
      fhicl::Atom<double> minADCdifference{Name("minADCdifference"), Comment("minimum ADC difference above pedestal to be considered for reconstruction")};  //5.0
      fhicl::Atom<double> defaultBeta{Name("defaultBeta"), Comment("initialization value for fit and default value for invalid fits (regular pulses: 19.0ns, dark counts for calibration: 12.0ns)")};
      fhicl::Atom<double> minBeta{Name("minBeta"), Comment("smallest accepted beta for valid fit [ns]")};  //5.0ns
      fhicl::Atom<double> maxBeta{Name("maxBeta"), Comment("largest accepted beta for valid fit [ns]")}; //50.0ns
      fhicl::Atom<double> maxTimeDifference{Name("maxTimeDifference"), Comment("largest accepted difference between time of largest ADC value and fitted peak [ns]")}; //20.0ns
      fhicl::Atom<double> minPulseHeightRatio{Name("minPulseHeightRatio"), Comment("smallest accepted ratio between largest ADC value and fitted peak")}; //0.7
      fhicl::Atom<double> maxPulseHeightRatio{Name("maxPulseHeightRatio"), Comment("largest accepted ratio between largest ADC value and fitted peak")}; //1.5
      fhicl::Atom<double> LEtimeFactor{Name("LEtimeFactor"), Comment("time of leading edge is peakTime-LEtimeFactor*beta (0.985,1.385,1.587 for a leading edge of 0.5,0.2,0.1 pulse height")};
      fhicl::Atom<int> minPEs{Name("minPEs"), Comment("minimum number of PEs")}; //0
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
    int         _minPEs;
    double      _digitizationPeriod;
    double      _pedestal;           //100 ADC
    double      _calibrationFactor;  //394.6 ADC*ns/PE
    double      _calibrationFactorPulseHeight;  //11.4 ADC/PE
    double      _microBunchPeriod;
  };


  CrvRecoPulsesFinder::CrvRecoPulsesFinder(const Parameters& conf) :
    art::EDProducer(conf),
    _crvDigiModuleLabel(conf().crvDigiModuleLabel()),
    _minPEs(conf().minPEs())
  {
    produces<CrvRecoPulseCollection>();
    _makeCrvRecoPulses=boost::shared_ptr<mu2eCrv::MakeCrvRecoPulses>(new mu2eCrv::MakeCrvRecoPulses(conf().minADCdifference(),
                                                                                                    conf().defaultBeta(),
                                                                                                    conf().minBeta(),
                                                                                                    conf().maxBeta(),
                                                                                                    conf().maxTimeDifference(),
                                                                                                    conf().minPulseHeightRatio(),
                                                                                                    conf().maxPulseHeightRatio(),
                                                                                                    conf().LEtimeFactor()));
  }

  void CrvRecoPulsesFinder::beginJob()
  {
  }

  void CrvRecoPulsesFinder::endJob()
  {
  }

  void CrvRecoPulsesFinder::beginRun(art::Run &run)
  {
    mu2e::ConditionsHandle<mu2e::AcceleratorParams> accPar("ignored");
    _microBunchPeriod = accPar->deBuncherPeriod;

    mu2e::ConditionsHandle<mu2e::CrvParams> crvPar("ignored");
    _digitizationPeriod = crvPar->digitizationPeriod;
    _pedestal           = crvPar->pedestal;
    _calibrationFactor  = crvPar->calibrationFactor;
    _calibrationFactorPulseHeight  = crvPar->calibrationFactorPulseHeight;
  }

  void CrvRecoPulsesFinder::produce(art::Event& event) 
  {
    std::unique_ptr<CrvRecoPulseCollection> crvRecoPulseCollection(new CrvRecoPulseCollection);

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
      for(size_t i=0; i<CrvDigi::NSamples; i++) ADCs.push_back(digi.GetADCs()[i]);
      waveformIndices.push_back(waveformIndex);

      //checking following digis whether they are a continuation of the current digis
      //if that is the case, append the next digis
      while(++waveformIndex<crvDigiCollection->size())
      {
        const CrvDigi &nextDigi = crvDigiCollection->at(waveformIndex);
        if(barIndex!=nextDigi.GetScintillatorBarIndex()) break;
        if(SiPM!=nextDigi.GetSiPMNumber()) break;
        if(startTDC+ADCs.size()!=nextDigi.GetStartTDC()) break;
        for(size_t i=0; i<CrvDigi::NSamples; i++) ADCs.push_back(nextDigi.GetADCs()[i]);
        waveformIndices.push_back(waveformIndex);
      }

      _makeCrvRecoPulses->SetWaveform(ADCs, startTDC, _digitizationPeriod, _pedestal, _calibrationFactor, _calibrationFactorPulseHeight);

      unsigned int n = _makeCrvRecoPulses->GetNPulses();
      for(unsigned int j=0; j<n; j++)
      {
        double pulseTime   = _makeCrvRecoPulses->GetPulseTime(j);
        int    PEs         = _makeCrvRecoPulses->GetPEs(j);
        int    PEsPulseHeight = _makeCrvRecoPulses->GetPEsPulseHeight(j);
        double pulseHeight = _makeCrvRecoPulses->GetPulseHeight(j); 
        double pulseBeta   = _makeCrvRecoPulses->GetPulseBeta(j);
        double pulseFitChi2= _makeCrvRecoPulses->GetPulseFitChi2(j);
        double LEtime      = _makeCrvRecoPulses->GetLEtime(j);
//        if(pulseTime<0) continue;
//        if(pulseTime>_microBunchPeriod) continue;
        if(PEs<_minPEs) continue; 
        crvRecoPulseCollection->emplace_back(PEs, PEsPulseHeight, pulseTime, pulseHeight, pulseBeta, pulseFitChi2, LEtime, 
                                                  waveformIndices, barIndex, SiPM);
      }
    }

    event.put(std::move(crvRecoPulseCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvRecoPulsesFinder;
DEFINE_ART_MODULE(CrvRecoPulsesFinder)
