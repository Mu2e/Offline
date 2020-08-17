//
// A module to extract leading edge times, pulse heights, integrals and number of photons from the CRV waveforms
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
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
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>

#include <TMath.h>

namespace mu2e 
{
  class CrvRecoPulsesFinder : public art::EDProducer 
  {

    public:
    explicit CrvRecoPulsesFinder(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    boost::shared_ptr<mu2eCrv::MakeCrvRecoPulses> _makeCrvRecoPulses;

    double      _digitizationPeriod;
    std::string _crvDigiModuleLabel;
    double      _pedestal;           //100 ADC
    double      _calibrationFactor;  //394.6 ADC*ns/PE
    double      _calibrationFactorPulseHeight;  //11.4 ADC/PE
    double      _minADCdifference;
    double      _defaultBeta;
    double      _minBeta;
    double      _maxBeta;
    double      _maxTimeDifference;
    double      _minPulseHeightRatio;
    double      _maxPulseHeightRatio;
    double      _LEtimeFactor;
    int         _minPEs;
    double      _microBunchPeriod;
  };

  CrvRecoPulsesFinder::CrvRecoPulsesFinder(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _crvDigiModuleLabel(pset.get<std::string>("crvDigiModuleLabel")),
    _minADCdifference(pset.get<double>("minADCdifference",5)),//minimum ADC difference above pedestal
                                                              //to be considered for reconstruction
    _defaultBeta(pset.get<double>("defaultBeta",19.0)), //19.0ns for regular pulses
                                                        //12.6ns for dark noise pulses used for calibration
                                                        //used for initialization of fit function
                                                        //and as default value for invalid fits
    _minBeta(pset.get<double>("minBeta",5.0)),          //5.0ns smallest accepted beta for valid fit
    _maxBeta(pset.get<double>("maxBeta",50.0)),         //50.0ns largest accepted beta for valid fit
    _maxTimeDifference(pset.get<double>("maxTimeDifference",20.0)), //20.0ns largest accepted differences between
                                                                    //time of largest ADC value and fitted peak
    _minPulseHeightRatio(pset.get<double>("minPulseHeightRatio",0.7)), //smallest accepted ratio between
                                                                       //largest ADC value and fitted peak
    _maxPulseHeightRatio(pset.get<double>("maxPulseHeightRatio",1.5)), //largest accepted ratio between
                                                                       //largest ADC value and fitted peak
    _LEtimeFactor(pset.get<double>("LEtimeFactor",0.985)), //time of leading edge is peakTime-LEtimeFactor*beta
                                                           //e.g. 0.985 for a leading edge at 50% pulse height
                                                           //e.g. 1.385 for a leading edge at 20% pulse height
                                                           //e.g. 1.587 for a leading edge at 10% pulse height
    _minPEs(pset.get<int>("minPEs",0))                  //0 PEs
  {
    produces<CrvRecoPulseCollection>();
    _makeCrvRecoPulses=boost::shared_ptr<mu2eCrv::MakeCrvRecoPulses>(new mu2eCrv::MakeCrvRecoPulses(_minADCdifference,
                                                                                                    _defaultBeta,
                                                                                                    _minBeta,
                                                                                                    _maxBeta,
                                                                                                    _maxTimeDifference,
                                                                                                    _minPulseHeightRatio,
                                                                                                    _maxPulseHeightRatio,
                                                                                                    _LEtimeFactor));
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
