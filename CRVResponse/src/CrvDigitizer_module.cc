//
// A module to convert waveform voltages into ADC counts
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CRVResponse/inc/MakeCrvDigis.hh"

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

#include "Offline/CRVConditions/inc/CRVADCRange.hh"
#include "Offline/CRVConditions/inc/CRVDigitizationPeriod.hh"
#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
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
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config
    {
      fhicl::Atom<std::string> crvWaveformsModuleLabel{Name("crvWaveformsModuleLabel")};
      fhicl::Atom<double>      ADCconversionFactor{Name("ADCconversionFactor")};
      fhicl::Atom<int>         pedestal{Name("pedestal")};
      fhicl::Atom<bool>        simulateNZS{Name("simulateNZS")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit CrvDigitizer(const Parameters& conf);
    void produce(art::Event& e);

    private:
    boost::shared_ptr<mu2eCrv::MakeCrvDigis> _makeCrvDigis;

    std::string _crvWaveformsModuleLabel;
    double      _ADCconversionFactor;
    int         _pedestal;
    bool        _simulateNZS;

    mu2e::ProditionsHandle<mu2e::CRVOrdinal> _channelMap_h;
  };

  CrvDigitizer::CrvDigitizer(const Parameters& conf) :
    art::EDProducer{conf},
    _crvWaveformsModuleLabel(conf().crvWaveformsModuleLabel()),
    _ADCconversionFactor(conf().ADCconversionFactor()),
    _pedestal(conf().pedestal()),
    _simulateNZS(conf().simulateNZS())
  {
    produces<CrvDigiCollection>();
    if(_simulateNZS) produces<CrvDigiCollection>("NZS");
    _makeCrvDigis = boost::shared_ptr<mu2eCrv::MakeCrvDigis>(new mu2eCrv::MakeCrvDigis());
  }

  void CrvDigitizer::produce(art::Event& event)
  {
    std::array<std::unique_ptr<CrvDigiCollection>,2> crvDigiCollection = {std::unique_ptr<CrvDigiCollection>(new CrvDigiCollection),   //ZS
                                                                          std::unique_ptr<CrvDigiCollection>(new CrvDigiCollection)};  //NZS
    std::array<art::Handle<CrvDigiMCCollection>,2> crvDigiMCCollection;
    event.getByLabel(_crvWaveformsModuleLabel,"",crvDigiMCCollection[0]);
    if(_simulateNZS) event.getByLabel(_crvWaveformsModuleLabel,"NZS",crvDigiMCCollection[1]);

    auto const& channelMap = _channelMap_h.get(event.id());

    for(int i=0; i<(_simulateNZS?2:1); ++i)
    {
      for(CrvDigiMCCollection::const_iterator iter=crvDigiMCCollection[i]->begin(); iter!=crvDigiMCCollection[i]->end(); iter++)
      {
        const CrvDigiMC &crvDigiMC = *iter;
        const CRSScintillatorBarIndex &barIndex = crvDigiMC.GetScintillatorBarIndex();
        const int SiPM = crvDigiMC.GetSiPMNumber();

        uint16_t offlineChannel = barIndex.asUint()*CRVId::nChanPerBar + SiPM;
        CRVROC   onlineChannel  = channelMap.online(offlineChannel);
        uint8_t  roc            = onlineChannel.ROC();
        uint8_t  feb            = onlineChannel.FEB();
        uint8_t  febChannel     = onlineChannel.FEBchannel();

        const std::vector<double> &voltages = crvDigiMC.GetVoltages();
        double startTime = crvDigiMC.GetStartTime();
        double TDC0time = crvDigiMC.GetTDC0Time();
        double NZS = crvDigiMC.IsNZS();

        //start time gets measured with respect to FEB's TDC near the event window start.
        //that's the TDC that gets set to TDC=0 and can be anywhere within: event window start <= t < event window start + digitization period
        //the waveform generator makes sure that the start time - TDC0time is a multiple of the digitization period
        startTime-=TDC0time;

        _makeCrvDigis->SetWaveform(voltages,_ADCconversionFactor,_pedestal, startTime, CRVDigitizationPeriod, CRVMinADC, CRVMaxADC);
        const std::vector<int16_t> &ADCs = _makeCrvDigis->GetADCs();
        uint16_t startTDC = _makeCrvDigis->GetTDC();

        crvDigiCollection[i]->emplace_back(ADCs, startTDC, NZS, false, barIndex, SiPM, roc, feb, febChannel);
      }

      event.put(std::move(crvDigiCollection[i]),(i==1?"NZS":""));
    }
  } // end produce

} // end namespace mu2e

using mu2e::CrvDigitizer;
DEFINE_ART_MODULE(CrvDigitizer)
