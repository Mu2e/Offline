//
// A module to create CRV waveforms from CRV SiPMResponses
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvWaveforms.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvSiPMResponsesCollection.hh"
#include "MCDataProducts/inc/CrvWaveformsCollection.hh"

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

#include "Randomize.hh"
#include <TMath.h>

namespace mu2e 
{
  class CrvWaveformsGenerator : public art::EDProducer 
  {

    public:
    explicit CrvWaveformsGenerator(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void endJob();

    private:
    std::string _crvSiPMResponsesModuleLabel;
    std::string _singlePEWaveformFileName;

    boost::shared_ptr<MakeCrvWaveforms> _makeCrvWaveforms;

    double  _binWidth;
  };

  CrvWaveformsGenerator::CrvWaveformsGenerator(fhicl::ParameterSet const& pset) :
    _crvSiPMResponsesModuleLabel(pset.get<std::string>("crvSiPMResponsesModuleLabel")),
    _singlePEWaveformFileName(pset.get<std::string>("singlePEWaveformFileName")),
    _binWidth(pset.get<double>("binWidth",12.5))                   //12.5 ns (digitizer sampling rate)
  {
    double singlePEWaveformBinWidth(pset.get<double>("singlePEWaveformBinWidth",1.0));    //1.0 ns
    int nBins(pset.get<int>("singlePEWaveformBins",200));          //200
    ConfigFileLookupPolicy configFile;
    _singlePEWaveformFileName = configFile(_singlePEWaveformFileName);
    _makeCrvWaveforms = boost::shared_ptr<MakeCrvWaveforms>(new MakeCrvWaveforms());
    _makeCrvWaveforms->LoadSinglePEWaveform(_singlePEWaveformFileName, singlePEWaveformBinWidth, nBins);
    produces<CrvWaveformsCollection>();
  }

  void CrvWaveformsGenerator::beginJob()
  {
  }

  void CrvWaveformsGenerator::endJob()
  {
  }

  void CrvWaveformsGenerator::produce(art::Event& event) 
  {
    std::unique_ptr<CrvWaveformsCollection> crvWaveformsCollection(new CrvWaveformsCollection);

    art::Handle<CrvSiPMResponsesCollection> crvSiPMResponsesCollection;
    event.getByLabel(_crvSiPMResponsesModuleLabel,"",crvSiPMResponsesCollection);

    for(CrvSiPMResponsesCollection::const_iterator iter=crvSiPMResponsesCollection->begin(); 
        iter!=crvSiPMResponsesCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CrvSiPMResponses &siPMResponses = iter->second;

      double startTime = siPMResponses.GetFirstSiPMResponseTime();
      startTime-=G4UniformRand()*_binWidth;

      CrvWaveforms &crvWaveforms = (*crvWaveformsCollection)[barIndex];
      crvWaveforms.SetBinWidth(_binWidth);
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<CrvSiPMResponses::CrvSingleSiPMResponse> &timesAndCharges = siPMResponses.GetSiPMResponses(SiPM);
        std::vector<double> times, charges;
        for(unsigned int i=0; i<timesAndCharges.size(); i++)
        {
          times.push_back(timesAndCharges[i]._time);
          charges.push_back(timesAndCharges[i]._charge);
        }
        std::vector<double> &waveform = crvWaveforms.GetWaveform(SiPM);
        _makeCrvWaveforms->MakeWaveform(times, charges, waveform, startTime, _binWidth);
        crvWaveforms.SetStartTime(SiPM,startTime);
      }
    }

    event.put(std::move(crvWaveformsCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvWaveformsGenerator;
DEFINE_ART_MODULE(CrvWaveformsGenerator)
