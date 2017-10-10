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
#include "MCDataProducts/inc/CrvWaveformsCollection.hh"
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

    art::Handle<CrvWaveformsCollection> crvWaveformsCollection;
    event.getByLabel(_crvWaveformsModuleLabel,"",crvWaveformsCollection);

    for(CrvWaveformsCollection::const_iterator iter=crvWaveformsCollection->begin(); 
        iter!=crvWaveformsCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CrvWaveforms &crvWaveforms = iter->second;
      double digitizationPrecision = crvWaveforms.GetDigitizationPrecision(); //ns

      CrvDigi &crvDigi = (*crvDigiCollection)[barIndex];
      crvDigi.SetDigitizationPrecision(digitizationPrecision);
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        std::vector<CrvDigi::CrvSingleWaveform> &singleWaveformsDigi = crvDigi.GetSingleWaveforms(SiPM);
        const std::vector<CrvWaveforms::CrvSingleWaveform> &singleWaveforms = crvWaveforms.GetSingleWaveforms(SiPM);
        for(size_t i=0; i<singleWaveforms.size(); i++)
        {
          //start new single waveform
          CrvDigi::CrvSingleWaveform singleWaveformDigi;

          singleWaveformDigi._startTime = singleWaveforms[i]._startTime; 

          const std::vector<double> &voltages = singleWaveforms[i]._voltages;
          _makeCrvDigis->SetWaveform(voltages,_ADCconversionFactor,_pedestal);
          singleWaveformDigi._ADCs = _makeCrvDigis->GetADCs();

          singleWaveformsDigi.push_back(singleWaveformDigi);
        }

      }
    }

    event.put(std::move(crvDigiCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvDigitizer;
DEFINE_ART_MODULE(CrvDigitizer)
