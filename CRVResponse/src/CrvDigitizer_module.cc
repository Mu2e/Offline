//
// A module to convert waveform voltages into ADC counts
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvDigis.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
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

    std::string _crvWaveformsModuleLabel;
    double      _ADCconversionFactor;
    int         _pedestal;
  };

  CrvDigitizer::CrvDigitizer(fhicl::ParameterSet const& pset) :
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel")),
    _ADCconversionFactor(pset.get<double>("ADCconversionFactor")),
    _pedestal(pset.get<int>("Pedestal"))
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
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CrvDigiMC &crvDigiMC = iter->second;
      double digitizationPrecision = crvDigiMC.GetDigitizationPrecision(); //ns //inverse of TDC rate

      CrvDigi &crvDigi = (*crvDigiCollection)[barIndex];
      crvDigi.SetDigitizationPrecision(digitizationPrecision);
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        std::vector<CrvDigi::CrvSingleWaveform> &singleWaveforms = crvDigi.GetSingleWaveforms(SiPM);
        const std::vector<CrvDigiMC::CrvSingleWaveform> &singleWaveformsMC = crvDigiMC.GetSingleWaveforms(SiPM);
        for(size_t i=0; i<singleWaveformsMC.size(); i++)
        {
          //start new single waveform
          CrvDigi::CrvSingleWaveform singleWaveform;

          const std::vector<double> &voltages = singleWaveformsMC[i]._voltages;
          double startTime = singleWaveformsMC[i]._startTime; 
          _makeCrvDigis->SetWaveform(voltages,_ADCconversionFactor,_pedestal, startTime, digitizationPrecision);
          singleWaveform._ADCs = _makeCrvDigis->GetADCs();
          singleWaveform._startTDC = _makeCrvDigis->GetTDC();

          singleWaveforms.push_back(singleWaveform);
        }

      }
    }

    event.put(std::move(crvDigiCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvDigitizer;
DEFINE_ART_MODULE(CrvDigitizer)
