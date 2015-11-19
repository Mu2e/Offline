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
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvWaveformsCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"

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

    std::string _crvWaveformsModuleLabel;
    double      _param0;
    double      _param1;
    double      _pulseThreshold;
    double      _leadingEdgeThreshold;
    double      _microBunchPeriod;
    int         _minPEs;
  };

  CrvRecoPulsesFinder::CrvRecoPulsesFinder(fhicl::ParameterSet const& pset) :
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel")),
    _param0(pset.get<double>("param0")),   //needs to be 0
    _param1(pset.get<double>("param1")),   //51.0
    _pulseThreshold(pset.get<double>("pulseThreshold")),  //0.015
    _leadingEdgeThreshold(pset.get<double>("leadingEdgeThreshold")),   //0.2
    _minPEs(pset.get<int>("minPEs"))   //3
  {
    produces<CrvRecoPulsesCollection>();
    _makeCrvRecoPulses = boost::shared_ptr<mu2eCrv::MakeCrvRecoPulses>(new mu2eCrv::MakeCrvRecoPulses(_pulseThreshold, _leadingEdgeThreshold, _param0, _param1));
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
  }

  void CrvRecoPulsesFinder::produce(art::Event& event) 
  {
    std::unique_ptr<CrvRecoPulsesCollection> crvRecoPulsesCollection(new CrvRecoPulsesCollection);

    art::Handle<CrvWaveformsCollection> crvWaveformsCollection;
    event.getByLabel(_crvWaveformsModuleLabel,"",crvWaveformsCollection);

    for(CrvWaveformsCollection::const_iterator iter=crvWaveformsCollection->begin(); 
        iter!=crvWaveformsCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CrvWaveforms &crvWaveforms = iter->second;
      double binWidth = crvWaveforms.GetBinWidth(); //ns

      CrvRecoPulses &crvRecoPulses = (*crvRecoPulsesCollection)[barIndex];
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<double> waveform = crvWaveforms.GetWaveform(SiPM);
        double startTime = crvWaveforms.GetStartTime(SiPM);
        _makeCrvRecoPulses->SetWaveform(waveform, startTime, binWidth);

        unsigned int n = _makeCrvRecoPulses->GetNPulses();
        for(unsigned int i=0; i<n; i++)
        {
          double time=_makeCrvRecoPulses->GetLeadingEdge(i);
          if(time<0) continue;
          if(time>_microBunchPeriod) continue;
          int PEs = _makeCrvRecoPulses->GetPEs(i);
          double height = _makeCrvRecoPulses->GetPulseHeight(i);
          if(PEs<_minPEs) continue; 
          crvRecoPulses.GetRecoPulses(SiPM).emplace_back(PEs, time,height);
        }
      }
    }

    event.put(std::move(crvRecoPulsesCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvRecoPulsesFinder;
DEFINE_ART_MODULE(CrvRecoPulsesFinder)
