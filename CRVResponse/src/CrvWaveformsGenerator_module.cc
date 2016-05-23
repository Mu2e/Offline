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
#include "SeedService/inc/SeedService.hh"

#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Random/Randomize.h"

#include <string>

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

    boost::shared_ptr<mu2eCrv::MakeCrvWaveforms> _makeCrvWaveforms;

    double                              _digitizationPrecision;
    double                              _FEBtimeSpread;
    double                              _minVoltage;

    CLHEP::RandFlat                     _randFlat;
    CLHEP::RandGaussQ                   _randGaussQ;
    
    std::vector<double> _timeShiftFEBsSide0, _timeShiftFEBsSide1;
  };

  CrvWaveformsGenerator::CrvWaveformsGenerator(fhicl::ParameterSet const& pset) :
    _crvSiPMResponsesModuleLabel(pset.get<std::string>("crvSiPMResponsesModuleLabel")),
    _singlePEWaveformFileName(pset.get<std::string>("singlePEWaveformFileName")),
    _digitizationPrecision(pset.get<double>("digitizationPrecision")),   //12.5 ns
    _FEBtimeSpread(pset.get<double>("FEBtimeSpread")),         //2.0 ns (due to cable lengths differences, etc.)
    _minVoltage(pset.get<double>("minVoltage")),               //0.022V (corresponds to 3.5PE)
    _randFlat(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _randGaussQ(art::ServiceHandle<art::RandomNumberGenerator>()->getEngine())
  {
    double singlePEWaveformPrecision(pset.get<double>("singlePEWaveformPrecision"));    //1.0 ns
    double singlePEWaveformMaxTime(pset.get<double>("singlePEWaveformMaxTime"));        //200
    ConfigFileLookupPolicy configFile;
    _singlePEWaveformFileName = configFile(_singlePEWaveformFileName);
    _makeCrvWaveforms = boost::shared_ptr<mu2eCrv::MakeCrvWaveforms>(new mu2eCrv::MakeCrvWaveforms());
    _makeCrvWaveforms->LoadSinglePEWaveform(_singlePEWaveformFileName, singlePEWaveformPrecision, singlePEWaveformMaxTime);
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

    double samplingPointShift = _randFlat.fire()*_digitizationPrecision;

    if(_timeShiftFEBsSide0.empty())
    {
      GeomHandle<CosmicRayShield> CRS;
      unsigned int nCounters = CRS->getAllCRSScintillatorBars().size();
      unsigned int nFEBs = ceil(nCounters/32.0);
      for(unsigned int i=0; i<nFEBs; i++)    
      {
        _timeShiftFEBsSide0.emplace_back(_randGaussQ.fire(0, _FEBtimeSpread));
        _timeShiftFEBsSide1.emplace_back(_randGaussQ.fire(0, _FEBtimeSpread));
      }
    }

    for(CrvSiPMResponsesCollection::const_iterator iter=crvSiPMResponsesCollection->begin(); 
        iter!=crvSiPMResponsesCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CrvSiPMResponses &siPMResponses = iter->second;

      double firstSiPMResponseTime = siPMResponses.GetFirstSiPMResponseTime();

      unsigned int FEB=barIndex.asUint()/32.0; //assume that the counters are ordered in the correct way, i.e. that all counters beloning to the same FEB are grouped together
      double minTimeShiftFEB = std::min(_timeShiftFEBsSide0[FEB],_timeShiftFEBsSide1[FEB]);
      firstSiPMResponseTime += std::min(minTimeShiftFEB,0.0);

      int startTime = floor(firstSiPMResponseTime / _digitizationPrecision) * _digitizationPrecision;  //start time of the waveform in multiples of the digitization interval (12.5ns)
      startTime -= samplingPointShift;  //random shift of start time (same shift for all FEBs of this event)

      CrvWaveforms &crvWaveforms = (*crvWaveformsCollection)[barIndex];
      crvWaveforms.SetBinWidth(_digitizationPrecision);
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        double timeShift=0;
        if(SiPM%2==0 && FEB<_timeShiftFEBsSide0.size()) timeShift=_timeShiftFEBsSide0[FEB];
        if(SiPM%2==1 && FEB<_timeShiftFEBsSide1.size()) timeShift=_timeShiftFEBsSide1[FEB];

        const std::vector<CrvSiPMResponses::CrvSingleSiPMResponse> &timesAndCharges = siPMResponses.GetSiPMResponses(SiPM);
        std::vector<double> times, charges;
        for(unsigned int i=0; i<timesAndCharges.size(); i++)
        {
          times.push_back(timesAndCharges[i]._time + timeShift);
          charges.push_back(timesAndCharges[i]._charge);
        }

        std::vector<double> &waveform = crvWaveforms.GetWaveform(SiPM);
        _makeCrvWaveforms->MakeWaveform(times, charges, waveform, startTime, _digitizationPrecision);
        crvWaveforms.SetStartTime(SiPM,startTime);

        //zero suppression, i.e. set all waveform digi points to zero which are below the minimum voltage, 
        //if the neighboring digi points are also below the minimum voltage
        for(unsigned int i=0; i<waveform.size(); i++)
        {
          bool keepPoint=false;
          if(waveform[i]>_minVoltage) keepPoint=true;
          else
          {
            if(i<waveform.size()-1)
            {
              if(waveform[i+1]>_minVoltage) keepPoint=true;
            }
          }
          if(!keepPoint) waveform[i]=0;
        }
      }
    }

    event.put(std::move(crvWaveformsCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvWaveformsGenerator;
DEFINE_ART_MODULE(CrvWaveformsGenerator)
