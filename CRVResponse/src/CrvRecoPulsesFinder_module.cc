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
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"

#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <TH1D.h>

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

    std::string _crvDigiModuleLabel;
    double      _calibrationFactor, _pedestal;
    bool        _usePulseArea;  //for PE calculation
    int         _minPEs;
    double      _microBunchPeriod;

//    TH1F        *_hRecoPulses;
  };

  CrvRecoPulsesFinder::CrvRecoPulsesFinder(fhicl::ParameterSet const& pset) :
    _crvDigiModuleLabel(pset.get<std::string>("crvDigiModuleLabel")),
    _calibrationFactor(pset.get<double>("calibrationFactor")),   //0.0056 V/PE
    _pedestal(pset.get<double>("pedestal")),   //0 V
    _usePulseArea(pset.get<bool>("usePulseArea")),   //still false, but will be changed in the future
    _minPEs(pset.get<int>("minPEs"))     //6 PEs
  {
    produces<CrvRecoPulsesCollection>();
    _makeCrvRecoPulses = boost::shared_ptr<mu2eCrv::MakeCrvRecoPulses>(new mu2eCrv::MakeCrvRecoPulses(_calibrationFactor, _pedestal, _usePulseArea));
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

    art::Handle<CrvDigiCollection> crvDigiCollection;
    event.getByLabel(_crvDigiModuleLabel,"",crvDigiCollection);

    for(CrvDigiCollection::const_iterator iter=crvDigiCollection->begin(); 
        iter!=crvDigiCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CrvDigi &crvWaveforms = iter->second;
      double digitizationPrecision = crvWaveforms.GetDigitizationPrecision(); //ns

      CrvRecoPulses &crvRecoPulses = (*crvRecoPulsesCollection)[barIndex];
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        //merge single waveforms together if there is no time gap between them
        const std::vector<CrvDigi::CrvSingleWaveform> &singleWaveforms = crvWaveforms.GetSingleWaveforms(SiPM);
        std::vector<std::vector<unsigned int> > ADCs;
        std::vector<std::vector<unsigned int> > singleWaveformIndices;
        std::vector<unsigned int> startTDCs;
        for(size_t i=0; i<singleWaveforms.size(); i++)
        {
          if(singleWaveforms[i]._ADCs.size()==0) continue;

          bool appendWaveform=false;
          if(!ADCs.empty())
          {
            //difference between the time of the last digitization point and the next start time
            unsigned int lastTDC = startTDCs.back()+(ADCs.back().size()-1);
            unsigned int nextTDC = singleWaveforms[i]._startTDC;
            if(lastTDC+1 == nextTDC) appendWaveform=true;   //the next start time seems to be just 
                                                            //one digitization point (12.5ns) away
                                                            //so that one can assume that the following 
                                                            //single waveform is just a continuation
                                                            //and can be appended.
          }

          if(appendWaveform) 
          {
            ADCs.back().insert(ADCs.back().end(),
                               singleWaveforms[i]._ADCs.begin(),
                               singleWaveforms[i]._ADCs.end());
            singleWaveformIndices.back().insert(singleWaveformIndices.back().end(),singleWaveforms[i]._ADCs.size(),i);
          }
          else
          {
            ADCs.push_back(singleWaveforms[i]._ADCs);
            std::vector<unsigned int> tmp(singleWaveforms[i]._ADCs.size(),i);
            singleWaveformIndices.push_back(tmp);
            startTDCs.push_back(singleWaveforms[i]._startTDC); 
          }
        }

        for(size_t i=0; i<ADCs.size(); i++)
        {
          _makeCrvRecoPulses->SetWaveform(ADCs[i], startTDCs[i], digitizationPrecision);

          unsigned int n = _makeCrvRecoPulses->GetNPulses();
          for(unsigned int j=0; j<n; j++)
          {
            double pulseTime   = _makeCrvRecoPulses->GetPulseTime(j);
            int    PEs         = _makeCrvRecoPulses->GetPEs(j);
            double pulseHeight = _makeCrvRecoPulses->GetPulseHeight(j); 
            double pulseWidth  = _makeCrvRecoPulses->GetPulseWidth(j);
            double pulseFitChi2= _makeCrvRecoPulses->GetPulseFitChi2(j);
            double LEtime      = _makeCrvRecoPulses->GetLEtime(j);
            int    peakBin     = _makeCrvRecoPulses->GetPeakBin(j);
//            if(pulseTime<0) continue;
//            if(pulseTime>_microBunchPeriod) continue;
            if(PEs<_minPEs) continue; 
            int singleWaveformIndex = singleWaveformIndices[i][peakBin];
            crvRecoPulses.GetRecoPulses(SiPM).emplace_back(PEs, pulseTime, pulseHeight, pulseWidth, pulseFitChi2, LEtime, singleWaveformIndex);
          }
        }

      }
    }

    event.put(std::move(crvRecoPulsesCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvRecoPulsesFinder;
DEFINE_ART_MODULE(CrvRecoPulsesFinder)
