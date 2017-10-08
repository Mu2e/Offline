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

#include "canvas/Persistency/Common/Ptr.h"
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
    int                                 _digitizationPoints;
    double                              _FEBtimeSpread;
    double                              _minVoltage;
    double                              _noise;

    CLHEP::RandFlat                     _randFlat;
    CLHEP::RandGaussQ                   _randGaussQ;
    
    std::vector<double> _timeShiftFEBsSide0, _timeShiftFEBsSide1;

    bool SingleWaveformStart(std::vector<double> &fullWaveform, size_t i);
  };

  CrvWaveformsGenerator::CrvWaveformsGenerator(fhicl::ParameterSet const& pset) :
    _crvSiPMResponsesModuleLabel(pset.get<std::string>("crvSiPMResponsesModuleLabel")),
    _singlePEWaveformFileName(pset.get<std::string>("singlePEWaveformFileName")),
    _digitizationPrecision(pset.get<double>("digitizationPrecision")),   //12.5 ns
    _digitizationPoints(pset.get<int>("digitizationPoints")),  //8 points for every single waveform
    _FEBtimeSpread(pset.get<double>("FEBtimeSpread")),         //2.0 ns (due to cable lengths differences, etc.)
    _minVoltage(pset.get<double>("minVoltage")),               //0.022V (corresponds to 3.5PE)
    _noise(pset.get<double>("noise")),
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

    GeomHandle<CosmicRayShield> CRS;
    _timeShiftFEBsSide0.clear();
    _timeShiftFEBsSide1.clear();
    unsigned int nCounters = CRS->getAllCRSScintillatorBars().size();
    unsigned int nFEBs = ceil(nCounters/32.0);
    for(unsigned int i=0; i<nFEBs; i++)    
    {
      _timeShiftFEBsSide0.emplace_back(_randGaussQ.fire(0, _FEBtimeSpread));
      _timeShiftFEBsSide1.emplace_back(_randGaussQ.fire(0, _FEBtimeSpread));
    }

    for(CrvSiPMResponsesCollection::const_iterator iter=crvSiPMResponsesCollection->begin(); 
        iter!=crvSiPMResponsesCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CrvSiPMResponses &siPMResponses = iter->second;

      CrvWaveforms &crvWaveforms = (*crvWaveformsCollection)[barIndex];
      crvWaveforms.SetDigitizationPrecision(_digitizationPrecision);

      unsigned int FEB=barIndex.asUint()/32.0; //assume that the counters are ordered in the correct way, 
                                               //i.e. that all counters beloning to the same FEB are grouped together

      for(int SiPM=0; SiPM<4; SiPM++)
      {
        double firstSiPMResponseTime = siPMResponses.GetFirstSiPMResponseTime(SiPM);
        if(isnan(firstSiPMResponseTime)) continue;

        double timeShiftFEB=0;
        if(SiPM%2==0 && FEB<_timeShiftFEBsSide0.size()) timeShiftFEB=_timeShiftFEBsSide0[FEB];
        if(SiPM%2==1 && FEB<_timeShiftFEBsSide1.size()) timeShiftFEB=_timeShiftFEBsSide1[FEB];

        firstSiPMResponseTime += timeShiftFEB;  //Ok, since all SiPMResponse times of this SiPM will be shifted by the same timeShiftFEB

        double startTime = floor(firstSiPMResponseTime / _digitizationPrecision) * _digitizationPrecision;  //start time of the waveform 
                                                                                                            //in multiples of the 
                                                                                                            //digitization interval (12.5ns)

        startTime -= samplingPointShift;  //random shift of start time (same shift for all FEBs of this event)

        const std::vector<CrvSiPMResponses::CrvSingleSiPMResponse> &timesAndCharges = siPMResponses.GetSiPMResponses(SiPM);
        std::vector<double> times, charges;
        for(size_t i=0; i<timesAndCharges.size(); i++)
        {
          times.push_back(timesAndCharges[i]._time + timeShiftFEB);
          charges.push_back(timesAndCharges[i]._chargeInPEs);   //FIXME: needs to be _charge
        }

        //first create the full waveform
        std::vector<double> fullWaveform;
        _makeCrvWaveforms->MakeWaveform(times, charges, fullWaveform, startTime, _digitizationPrecision);
        _makeCrvWaveforms->AddElectronicNoise(fullWaveform, _noise, _randGaussQ);

        //break the waveform apart into short pieces (_digitizationPoints)
        //and apply the zero suppression, i.e. set all waveform digi points to zero which are below the minimum voltage, 
        //if the neighboring digi points are also below the minimum voltage
        std::vector<CrvWaveforms::CrvSingleWaveform> &singleWaveforms = crvWaveforms.GetSingleWaveforms(SiPM);

        for(size_t i=0; i<fullWaveform.size(); i++)
        {
          if(SingleWaveformStart(fullWaveform, i)) //acts as a zero suppression
          {
            //start new single waveform
            CrvWaveforms::CrvSingleWaveform singleWaveform;
            singleWaveform._startTime=startTime+i*_digitizationPrecision;
            for(int singleWaveformIndex=0; 
                i<fullWaveform.size() && singleWaveformIndex<_digitizationPoints; 
                i++, singleWaveformIndex++)
            {
              singleWaveform._voltages.push_back(fullWaveform[i]);
            }
            i--;
            singleWaveforms.push_back(singleWaveform);
          }
        }
      }
    }

    event.put(std::move(crvWaveformsCollection));
  } // end produce

  bool CrvWaveformsGenerator::SingleWaveformStart(std::vector<double> &fullWaveform, size_t i)
  {
    if(fullWaveform[i]>_minVoltage) return true;  //this point is above the threshold --> start recording

    //record at least two points before and after a point above the zero suppression threshold to help with the peak reconstruction
    if(i+2<fullWaveform.size())
    {
      if(fullWaveform[i+2]>_minVoltage) return true;  //the point following the next point is above the threshold --> start recording
    }

    if(i+1<fullWaveform.size())
    {
      if(fullWaveform[i+1]>_minVoltage) return true;  //the following point is above the threshold --> start recording
    }

    if(i-1>=0)
    {
      if(fullWaveform[i-1]>_minVoltage) return true;  //the previous point was above the threshold --> continue recording
    }

    if(i-2>=0)
    {
      if(fullWaveform[i-2]>_minVoltage) return true;  //the point before the previous point was above the threshold --> continue recording
    }

    return false;
  }

} // end namespace mu2e

using mu2e::CrvWaveformsGenerator;
DEFINE_ART_MODULE(CrvWaveformsGenerator)
