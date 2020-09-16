//
// A module to convert waveform voltages into ADC counts
//
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvDigis.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CrvParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
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
    mu2e::ConditionsHandle<mu2e::CrvParams> crvPar("ignored");
    _digitizationPeriod  = crvPar->digitizationPeriod;
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
      const double startTime = crvDigiMC.GetStartTime();

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
