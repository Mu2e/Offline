//
// A module to extract number of PEs, arrival times, hit positions, etc. from the CRV waveforms
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/CrvRecoPulseResponse.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CRVWaveformsCollection.hh"
#include "RecoDataProducts/inc/CRVRecoPulsesCollection.hh"

#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>

#include <TRandom3.h>
#include <TMath.h>

namespace mu2e 
{
  class ExtractCrvRecoPulses : public art::EDProducer 
  {

    public:
    explicit ExtractCrvRecoPulses(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void endJob();

    private:
    boost::shared_ptr<CrvRecoPulseResponse> _crvRecoPulseResponse;

    std::string _crvWaveformsModuleLabel;
    double      _integralFactor;
    double      _pulseThreshold;
    double      _leadingEdgeThreshold;
  };

  ExtractCrvRecoPulses::ExtractCrvRecoPulses(fhicl::ParameterSet const& pset) :
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel")),
    _integralFactor(pset.get<double>("integralFactor",51.0)),
    _pulseThreshold(pset.get<double>("pulseThreshold",0.005)),
    _leadingEdgeThreshold(pset.get<double>("leadingEdgeThreshold",0.2))
  {
    produces<CRVRecoPulsesCollection>();
    _crvRecoPulseResponse = boost::shared_ptr<CrvRecoPulseResponse>(new CrvRecoPulseResponse(_pulseThreshold, _leadingEdgeThreshold, _integralFactor));
  }

  void ExtractCrvRecoPulses::beginJob()
  {
  }

  void ExtractCrvRecoPulses::endJob()
  {
  }

  void ExtractCrvRecoPulses::produce(art::Event& event) 
  {
    std::unique_ptr<CRVRecoPulsesCollection> crvRecoPulsesCollection(new CRVRecoPulsesCollection);

    art::Handle<CRVWaveformsCollection> crvWaveformsCollection;
    event.getByLabel(_crvWaveformsModuleLabel,"",crvWaveformsCollection);

    for(CRVWaveformsCollection::const_iterator iter=crvWaveformsCollection->begin(); 
        iter!=crvWaveformsCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CRVWaveforms &crvWaveforms = iter->second;
      double binWidth = crvWaveforms.GetBinWidth(); //ns

      CRVRecoPulses &crvRecoPulses = (*crvRecoPulsesCollection)[barIndex];
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<double> waveform = crvWaveforms.GetWaveform(SiPM);
        double startTime = crvWaveforms.GetStartTime(SiPM);
        _crvRecoPulseResponse->SetWaveform(waveform, startTime, binWidth);

        unsigned int n = _crvRecoPulseResponse->GetNPulses();
        for(unsigned int i=0; i<n; i++)
        {
          CRVRecoPulses::CRVSingleRecoPulse pulse;
          pulse._PEs = _crvRecoPulseResponse->GetPEs(i);
          pulse._leadingEdge = _crvRecoPulseResponse->GetLeadingEdge(i);
          pulse._pulseHeight = _crvRecoPulseResponse->GetPulseHeight(i);
          crvRecoPulses.GetRecoPulses(SiPM).push_back(pulse);
        }
      }
    }

    event.put(std::move(crvRecoPulsesCollection));
  } // end produce

} // end namespace mu2e

using mu2e::ExtractCrvRecoPulses;
DEFINE_ART_MODULE(ExtractCrvRecoPulses)
