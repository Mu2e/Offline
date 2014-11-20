//
// A module to extract number of PEs, arrival times, hit positions, etc. from the CRV waveforms
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

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
    std::string _crvWaveformsModuleLabel;
    double      _integralFactor;
    double      _threshold;
  };

  ExtractCrvRecoPulses::ExtractCrvRecoPulses(fhicl::ParameterSet const& pset) :
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel")),
    _integralFactor(pset.get<double>("integralFactor",51.0)),
    _threshold(pset.get<double>("threshold",0.005))
  {
    produces<CRVRecoPulsesCollection>();
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
        unsigned int nBins = waveform.size();
        double time = crvWaveforms.GetStartTime(SiPM);
        for(unsigned bin=0; bin<nBins; bin++, time+=binWidth)
        {
          double voltage = waveform[bin];
          if(voltage>_threshold)
          {
            double leadingEdge = time;
            double integral = 0;
            double maxVoltage = 0;
            for( ; bin<nBins; bin++, time+=binWidth)
            {
              voltage = waveform[bin];
              if(voltage<_threshold) break;
              integral += voltage;
              if(voltage>maxVoltage) maxVoltage=voltage;
            }

            CRVRecoPulses::CRVSingleRecoPulse pulse;
            pulse._PEs = static_cast<int>(integral*_integralFactor+0.5);
            pulse._leadingEdge = leadingEdge;
            pulse._timeOverThreshold = time - leadingEdge;
            pulse._pulseHeight = maxVoltage;
            crvRecoPulses.GetRecoPulses(SiPM).push_back(pulse);
          }
        }
      }
    }

    event.put(std::move(crvRecoPulsesCollection));
  } // end produce

} // end namespace mu2e

using mu2e::ExtractCrvRecoPulses;
DEFINE_ART_MODULE(ExtractCrvRecoPulses)
