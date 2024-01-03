//
// A module to convert waveform voltages into ADC counts
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CRVResponse/inc/MakeCrvDigis.hh"

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

#include "Offline/CRVConditions/inc/CRVDigitizationPeriod.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/MCDataProducts/inc/CrvDigiMC.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
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

    std::string _crvWaveformsModuleLabel;
    double      _ADCconversionFactor;
    int         _pedestal;
  };

  CrvDigitizer::CrvDigitizer(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel")),
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
  }

  void CrvDigitizer::produce(art::Event& event)
  {
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
      double TDC0time = crvDigiMC.GetTDC0Time();

      //start time gets measured with respect to FEB's TDC near the event window start.
      //that's the TDC that gets set to TDC=0 and can be anywhere within: event window start <= t < event window start + digitization period
      //the waveform generator makes sure that the start time - TDC0time is a multiple of the digitization period
      startTime-=TDC0time;

      std::vector<double> voltageVector;
      for(size_t i=0; i<voltages.size(); i++) voltageVector.push_back(voltages[i]);

      _makeCrvDigis->SetWaveform(voltageVector,_ADCconversionFactor,_pedestal, startTime, CRVDigitizationPeriod);
      const std::vector<int16_t> &ADCs = _makeCrvDigis->GetADCs();
      uint16_t startTDC = _makeCrvDigis->GetTDC();

      std::array<int16_t, CrvDigi::NSamples> ADCArray={};
      for(size_t i=0; i<ADCs.size() && i<CrvDigi::NSamples; i++) ADCArray[i]=ADCs[i];

      crvDigiCollection->emplace_back(ADCArray, startTDC, barIndex, SiPM);
    }

    event.put(std::move(crvDigiCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvDigitizer;
DEFINE_ART_MODULE(CrvDigitizer)
