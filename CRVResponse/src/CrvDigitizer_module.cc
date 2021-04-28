//
// A module to convert waveform voltages into ADC counts
//
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvDigis.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "DataProducts/inc/EventWindowMarker.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CrvParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "DAQConditions/inc/EventTiming.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvDigiMC.hh"
#include "MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"

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
  class CrvDigitizer : public art::EDProducer 
  {

    public:
    explicit CrvDigitizer(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    boost::shared_ptr<mu2eCrv::MakeCrvDigis> _makeCrvDigis;

    double      _digitizationPeriod;
    std::string _crvWaveformsModuleLabel;
    double      _digitizationStart, _digitizationEnd;
    std::string _eventWindowMarkerLabel;
    std::string _protonBunchTimeMCLabel;
    double      _ADCconversionFactor;
    int         _pedestal;
  };

  CrvDigitizer::CrvDigitizer(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel")),
    _digitizationStart(pset.get<double>("digitizationStart")),       //400ns
    _digitizationEnd(pset.get<double>("digitizationEnd")),           //1750ns
    _eventWindowMarkerLabel(pset.get<std::string>("eventWindowMarker","EWMProducer")),
    _protonBunchTimeMCLabel(pset.get<std::string>("protonBunchTimeMC","EWMProducer")),
    _ADCconversionFactor(pset.get<double>("ADCconversionFactor")),
    _pedestal(pset.get<int>("pedestal"))
  {
    produces<CrvDigiCollection>();
    _makeCrvDigis = boost::shared_ptr<mu2eCrv::MakeCrvDigis>(new mu2eCrv::MakeCrvDigis());
  }

  void CrvDigitizer::beginJob()
  {
  }

  void CrvDigitizer::endJob()
  {
  }

  void CrvDigitizer::beginRun(art::Run &run)
  {
    mu2e::ConditionsHandle<mu2e::CrvParams> crvPar("ignored");
    _digitizationPeriod  = crvPar->digitizationPeriod;
  }

  void CrvDigitizer::produce(art::Event& event) 
  {
    art::Handle<EventWindowMarker> eventWindowMarker;
    event.getByLabel(_eventWindowMarkerLabel,"",eventWindowMarker);
    EventWindowMarker::SpillType spillType = eventWindowMarker->spillType();

    art::Handle<ProtonBunchTimeMC> protonBunchTimeMC;
    event.getByLabel(_protonBunchTimeMCLabel, protonBunchTimeMC);
    double TDC0time = -protonBunchTimeMC->pbtime_;

    ProditionsHandle<EventTiming> eventTimingHandle;
    const EventTiming &eventTiming = eventTimingHandle.get(event.id());
    double jitter = TDC0time - eventTiming.timeFromProtonsToDRMarker();

    double digitizationStart=_digitizationStart+jitter;
    double digitizationEnd=_digitizationEnd+jitter;
    if(spillType!=EventWindowMarker::SpillType::onspill)
    {
      double eventWindowLength = eventWindowMarker->eventLength();
      digitizationStart = TDC0time;
      digitizationEnd = digitizationStart + eventWindowLength;
    }

    std::unique_ptr<CrvDigiCollection> crvDigiCollection(new CrvDigiCollection);

    art::Handle<CrvDigiMCCollection> crvDigiMCCollection;
    event.getByLabel(_crvWaveformsModuleLabel,"",crvDigiMCCollection);

    for(CrvDigiMCCollection::const_iterator iter=crvDigiMCCollection->begin(); 
        iter!=crvDigiMCCollection->end(); iter++)
    {
      const CrvDigiMC &crvDigiMC = *iter;
      const CRSScintillatorBarIndex &barIndex = crvDigiMC.GetScintillatorBarIndex();
      const int SiPM = crvDigiMC.GetSiPMNumber();
      const std::array<double,CrvDigiMC::NSamples> &voltages = crvDigiMC.GetVoltages();
      double startTime = crvDigiMC.GetStartTime();

      //waveform get recorded only between digitization start and end
      if(startTime<digitizationStart || startTime>digitizationEnd) continue;

      //start time gets measured with respect to the event window start
      startTime-=TDC0time;
      if(startTime<0) continue; //this shouldn't happen

      std::vector<double> voltageVector;
      for(size_t i=0; i<voltages.size(); i++) voltageVector.push_back(voltages[i]); 

      _makeCrvDigis->SetWaveform(voltageVector,_ADCconversionFactor,_pedestal, startTime, _digitizationPeriod);
      const std::vector<unsigned int> &ADCs = _makeCrvDigis->GetADCs();
      int startTDC = _makeCrvDigis->GetTDC();

      std::array<unsigned int, CrvDigi::NSamples> ADCArray;
      for(size_t i=0; i<ADCs.size(); i++) ADCArray[i]=ADCs[i]; 

      crvDigiCollection->emplace_back(ADCArray, startTDC, barIndex, SiPM);
    }

    event.put(std::move(crvDigiCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvDigitizer;
DEFINE_ART_MODULE(CrvDigitizer)
