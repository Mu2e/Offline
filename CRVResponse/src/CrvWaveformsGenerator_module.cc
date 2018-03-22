//
// A module to create CRV waveforms from CRV SiPMCharges
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
#include "MCDataProducts/inc/CrvSiPMChargesCollection.hh"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
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
    std::string _crvSiPMChargesModuleLabel;
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
    _crvSiPMChargesModuleLabel(pset.get<std::string>("crvSiPMChargesModuleLabel")),
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
    double singlePEReferenceCharge(pset.get<double>("singlePEReferenceCharge")); //1.8564e-13 C (the charge which was used to generate the above 1PE waveform)
    ConfigFileLookupPolicy configFile;
    _singlePEWaveformFileName = configFile(_singlePEWaveformFileName);
    _makeCrvWaveforms = boost::shared_ptr<mu2eCrv::MakeCrvWaveforms>(new mu2eCrv::MakeCrvWaveforms());
    _makeCrvWaveforms->LoadSinglePEWaveform(_singlePEWaveformFileName, singlePEWaveformPrecision, singlePEWaveformMaxTime, singlePEReferenceCharge);
    produces<CrvDigiMCCollection>();
  }

  void CrvWaveformsGenerator::beginJob()
  {
  }

  void CrvWaveformsGenerator::endJob()
  {
  }

  void CrvWaveformsGenerator::produce(art::Event& event) 
  {
    std::unique_ptr<CrvDigiMCCollection> crvDigiMCCollection(new CrvDigiMCCollection);

    art::Handle<CrvSiPMChargesCollection> crvSiPMChargesCollection;
    event.getByLabel(_crvSiPMChargesModuleLabel,"",crvSiPMChargesCollection);

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

    for(CrvSiPMChargesCollection::const_iterator iter=crvSiPMChargesCollection->begin(); 
        iter!=crvSiPMChargesCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CrvSiPMCharges &siPMCharges = iter->second;

      CrvDigiMC &crvDigiMC = (*crvDigiMCCollection)[barIndex];
      crvDigiMC.SetDigitizationPrecision(_digitizationPrecision);
      bool empty=true;

      unsigned int FEB=barIndex.asUint()/32.0; //assume that the counters are ordered in the correct way, 
                                               //i.e. that all counters beloning to the same FEB are grouped together

      for(int SiPM=0; SiPM<4; SiPM++)
      {
        double firstSiPMChargeTime = siPMCharges.GetFirstSiPMChargeTime(SiPM);
        if(isnan(firstSiPMChargeTime)) continue;

        double timeShiftFEB=0;
        if(SiPM%2==0 && FEB<_timeShiftFEBsSide0.size()) timeShiftFEB=_timeShiftFEBsSide0[FEB];
        if(SiPM%2==1 && FEB<_timeShiftFEBsSide1.size()) timeShiftFEB=_timeShiftFEBsSide1[FEB];

        firstSiPMChargeTime += timeShiftFEB;  //Ok, since all SiPMCharge times of this SiPM will be shifted by the same timeShiftFEB

        double startTime = floor(firstSiPMChargeTime / _digitizationPrecision) * _digitizationPrecision;  //start time of the waveform 
                                                                                                          //in multiples of the 
                                                                                                          //digitization interval (12.5ns)

        startTime -= samplingPointShift;  //random shift of start time (same shift for all FEBs of this event)

        const std::vector<CrvSiPMCharges::CrvSingleCharge> &timesAndCharges = siPMCharges.GetSiPMCharges(SiPM);
        std::vector<double> times, charges;
        for(size_t i=0; i<timesAndCharges.size(); i++)
        {
          times.push_back(timesAndCharges[i]._time + timeShiftFEB);
          charges.push_back(timesAndCharges[i]._charge); 
        }

        //first create the full waveform
        std::vector<double> fullWaveform;
        _makeCrvWaveforms->MakeWaveform(times, charges, fullWaveform, startTime, _digitizationPrecision);
        _makeCrvWaveforms->AddElectronicNoise(fullWaveform, _noise, _randGaussQ);

        //break the waveform apart into short pieces (_digitizationPoints)
        //and apply the zero suppression, i.e. set all waveform digi points to zero which are below the minimum voltage, 
        //if the neighboring digi points are also below the minimum voltage
        std::vector<CrvDigiMC::CrvSingleWaveform> &singleWaveforms = crvDigiMC.GetSingleWaveforms(SiPM);

        for(size_t i=0; i<fullWaveform.size(); i++)
        {
          if(SingleWaveformStart(fullWaveform, i)) //acts as a zero suppression
          {
            //start new single waveform
            CrvDigiMC::CrvSingleWaveform singleWaveform;
            singleWaveform._startTime=startTime+i*_digitizationPrecision;
            //collect voltages
            for(int singleWaveformIndex=0; singleWaveformIndex<_digitizationPoints; i++, singleWaveformIndex++)
            {
              if(i<fullWaveform.size()) singleWaveform._voltages.push_back(fullWaveform[i]);
              else singleWaveform._voltages.push_back(0);   //so that all single waveforms have always the right number of entries
            }

            //collect StepPointMCs and SimParticles responsible for this single waveform
            std::set<art::Ptr<StepPointMC> > steps;
            std::map<art::Ptr<SimParticle>, int> simparticles;
            for(size_t j=0; j<timesAndCharges.size(); j++)
            {
              if(timesAndCharges[j]._time>=singleWaveform._startTime-50.0 && timesAndCharges[j]._time>=singleWaveform._startTime+50.0)  //FIXME
              {
                steps.insert(timesAndCharges[j]._step);
                if(timesAndCharges[j]._step.isNonnull()) simparticles[timesAndCharges[j]._step->simParticle()]++;
              }
            }
            //loop through the steps to fill the single waveform
            std::set<art::Ptr<StepPointMC> >::iterator stepIter;
            for(stepIter=steps.begin(); stepIter!=steps.end(); stepIter++)
            {
              singleWaveform._steps.push_back(*stepIter);
            }
            //find the most likely SimParticle
            //if no SimParticle was recorded for this single waveform, then it was caused either by noise hits (if the threshold is low enough), or is the tail end of the peak.
            //in that case, _simparticle will be null (set by the default constructor of art::Ptr)
            std::map<art::Ptr<SimParticle>,int >::iterator simparticleIter;
            int simparticleCount=0;
            for(simparticleIter=simparticles.begin(); simparticleIter!=simparticles.end(); simparticleIter++)
            {
              if(simparticleIter->second>simparticleCount)
              {
                simparticleCount=simparticleIter->second;
                singleWaveform._simparticle=simparticleIter->first;
              }
            }

            i--;
            singleWaveforms.push_back(singleWaveform);
            empty=false;
          }
        }
      } //SiPM

      //2 options:
      //(1) -create a crvDigiMC object as a reference to a crvDigiMCCollection map entry at the beginning for all counters
      //    -fill this crvDigiMC object
      //    -if the crvDigiMC object stays empty, erase the map entry in crvDigiMCCollection
      //(2) -create a standalone crvDigiMC object
      //    -fill this crvDigiMC object
      //    -if the crvDigiMC didn't stay empty, create a new map entry in crvDigiMCCollection and fill its content with
      //     the new crvDigiMC  <---- too time consuming, therefore use option (1)

      if(empty) crvDigiMCCollection->erase(barIndex);  
    }

    event.put(std::move(crvDigiMCCollection));
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

    if(i>=1)
    {
      if(fullWaveform[i-1]>_minVoltage) return true;  //the previous point was above the threshold --> continue recording
    }

    if(i>=2)
    {
      if(fullWaveform[i-2]>_minVoltage) return true;  //the point before the previous point was above the threshold --> continue recording
    }

    return false;
  }

} // end namespace mu2e

using mu2e::CrvWaveformsGenerator;
DEFINE_ART_MODULE(CrvWaveformsGenerator)
