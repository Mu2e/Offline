//
// A module to create CRV waveforms from CRV PEs
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/CrvWaveformResponse.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CRVPEsCollection.hh"
#include "MCDataProducts/inc/CRVWaveformsCollection.hh"

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
  class MakeCrvWaveforms : public art::EDProducer 
  {

    public:
    explicit MakeCrvWaveforms(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void endJob();

    private:
    std::string _crvPEsModuleLabel;
    std::string _singlePEWaveformFileName;

    boost::shared_ptr<CrvWaveformResponse> _crvWaveformResponse;

    double  _binWidth;
  };

  MakeCrvWaveforms::MakeCrvWaveforms(fhicl::ParameterSet const& pset) :
    _crvPEsModuleLabel(pset.get<std::string>("crvPEsModuleLabel")),
    _singlePEWaveformFileName(pset.get<std::string>("singlePEWaveformFileName")),
    _binWidth(pset.get<double>("binWidth",12.5))                   //12.5 ns (digitizer sampling rate)
  {
    double singlePEWaveformBinWidth(pset.get<double>("singlePEWaveformBinWidth",1.0));    //1.0 ns
    int nBins(pset.get<int>("singlePEWaveformBins",200));          //200
    _crvWaveformResponse = boost::shared_ptr<CrvWaveformResponse>(new CrvWaveformResponse());
    _crvWaveformResponse->LoadSinglePEWaveform(_singlePEWaveformFileName, singlePEWaveformBinWidth, nBins);
    produces<CRVWaveformsCollection>();
  }

  void MakeCrvWaveforms::beginJob()
  {
  }

  void MakeCrvWaveforms::endJob()
  {
  }

  void MakeCrvWaveforms::produce(art::Event& event) 
  {
    std::unique_ptr<CRVWaveformsCollection> crvWaveformsCollection(new CRVWaveformsCollection);

    art::Handle<CRVPEsCollection> crvPEsCollection;
    event.getByLabel(_crvPEsModuleLabel,"",crvPEsCollection);

    double startTime=NAN;
    for(CRVPEsCollection::const_iterator iter=crvPEsCollection->begin(); iter!=crvPEsCollection->end(); iter++)
    {
      const CRVPEs &crvPEs = iter->second;
      double t = crvPEs.GetFirstPEtime();
      if(isnan(startTime) || t<startTime) startTime=t;
    }
    startTime-=G4UniformRand()*_binWidth;

    for(CRVPEsCollection::const_iterator iter=crvPEsCollection->begin(); iter!=crvPEsCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CRVPEs &crvPEs = iter->second;

      CRVWaveforms &crvWaveforms = (*crvWaveformsCollection)[barIndex];
      crvWaveforms.SetBinWidth(_binWidth);
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<double> &arrivalTimes = crvPEs.GetPEtimes(SiPM);
        std::vector<double> &waveform = crvWaveforms.GetWaveform(SiPM);
        _crvWaveformResponse->MakeWaveforms(arrivalTimes, waveform, startTime, _binWidth);
        crvWaveforms.SetStartTime(SiPM,startTime);
      }
    }

    event.put(std::move(crvWaveformsCollection));
  } // end produce

} // end namespace mu2e

using mu2e::MakeCrvWaveforms;
DEFINE_ART_MODULE(MakeCrvWaveforms)
