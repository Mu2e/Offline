//
// A module to create CRV waveforms from CRV SiPMCharges
//
//
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvWaveforms.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/CrvParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvSiPMCharges.hh"
#include "MCDataProducts/inc/CrvDigiMC.hh"
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
    void beginRun(art::Run &run);

    private:
    std::string _crvSiPMChargesModuleLabel;
    std::string _singlePEWaveformFileName;
 
    boost::shared_ptr<mu2eCrv::MakeCrvWaveforms> _makeCrvWaveforms;

    double                              _digitizationPeriod;
    double                              _FEBtimeSpread;
    double                              _minVoltage;
    double                              _noise;
    double                              _singlePEWaveformMaxTime;

    CLHEP::HepRandomEngine&             _engine;
    CLHEP::RandFlat                     _randFlat;
    CLHEP::RandGaussQ                   _randGaussQ;

    
    std::vector<double> _timeShiftFEBsSide0, _timeShiftFEBsSide1;


    bool SingleWaveformStart(std::vector<double> &fullWaveform, size_t i);
  };

  CrvWaveformsGenerator::CrvWaveformsGenerator(fhicl::ParameterSet const& pset) :
    EDProducer{pset},
    _crvSiPMChargesModuleLabel(pset.get<std::string>("crvSiPMChargesModuleLabel")),
    _singlePEWaveformFileName(pset.get<std::string>("singlePEWaveformFileName")),
    _FEBtimeSpread(pset.get<double>("FEBtimeSpread")),         //2.0 ns (due to cable lengths differences, etc.)
    _minVoltage(pset.get<double>("minVoltage")),               //0.022V (corresponds to 3.5PE)
    _noise(pset.get<double>("noise")),

    _singlePEWaveformMaxTime(pset.get<double>("singlePEWaveformMaxTime")),        //100ns
    _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())},
    _randFlat{_engine},
    _randGaussQ{_engine}
  {
    double singlePEWaveformPrecision(pset.get<double>("singlePEWaveformPrecision"));    //1.0 ns
    double singlePEWaveformStretchFactor(pset.get<double>("singlePEWaveformStretchFactor"));    //1.047
    double singlePEReferenceCharge(pset.get<double>("singlePEReferenceCharge")); //2.652e-13 C (the charge which was used to generate the above 1PE waveform)
    ConfigFileLookupPolicy configFile;
    _singlePEWaveformFileName = configFile(_singlePEWaveformFileName);
    _makeCrvWaveforms = boost::shared_ptr<mu2eCrv::MakeCrvWaveforms>(new mu2eCrv::MakeCrvWaveforms());
    _makeCrvWaveforms->LoadSinglePEWaveform(_singlePEWaveformFileName, singlePEWaveformPrecision, singlePEWaveformStretchFactor, 
                                            _singlePEWaveformMaxTime, singlePEReferenceCharge);
    produces<CrvDigiMCCollection>();
  }

  void CrvWaveformsGenerator::beginRun(art::Run &run)
  {
    mu2e::ConditionsHandle<mu2e::CrvParams> crvPar("ignored");
    _digitizationPeriod  = crvPar->digitizationPeriod;
  }

  void CrvWaveformsGenerator::produce(art::Event& event)
  {
    std::unique_ptr<CrvDigiMCCollection> crvDigiMCCollection(new CrvDigiMCCollection);

    art::Handle<CrvSiPMChargesCollection> crvSiPMChargesCollection;
    event.getByLabel(_crvSiPMChargesModuleLabel,"",crvSiPMChargesCollection);

    double samplingPointShift = _randFlat.fire()*_digitizationPeriod;

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
      int SiPM = iter->GetSiPMNumber();
      CRSScintillatorBarIndex barIndex = iter->GetScintillatorBarIndex();
      unsigned int FEB=barIndex.asUint()/32.0; //assume that the counters are ordered in the correct way,
                                               //i.e. that all counters beloning to the same FEB are grouped together

      const std::vector<CrvSiPMCharges::SingleCharge> &timesAndCharges = iter->GetCharges();
      double firstChargeTime = NAN;
      for(size_t i=0; i<timesAndCharges.size(); i++)
      {
        if(isnan(firstChargeTime) || firstChargeTime>timesAndCharges[i]._time) firstChargeTime=timesAndCharges[i]._time;
      }
      if(isnan(firstChargeTime)) continue;

      double timeShiftFEB=0;
      if(SiPM%2==0 && FEB<_timeShiftFEBsSide0.size()) timeShiftFEB=_timeShiftFEBsSide0[FEB];
      if(SiPM%2==1 && FEB<_timeShiftFEBsSide1.size()) timeShiftFEB=_timeShiftFEBsSide1[FEB];

      firstChargeTime += timeShiftFEB;  //Ok, since all SiPMCharge times of this SiPM will be shifted by the same timeShiftFEB

      double startTime = floor(firstChargeTime / _digitizationPeriod) * _digitizationPeriod;  //start time of the waveform
                                                                                              //in multiples of the
                                                                                              //digitization period (12.55ns)

      startTime -= samplingPointShift;  //random shift of start time (same shift for all FEBs of this event)

      std::vector<double> times, charges;
      for(size_t i=0; i<timesAndCharges.size(); i++)
      {
        times.push_back(timesAndCharges[i]._time + timeShiftFEB);
        charges.push_back(timesAndCharges[i]._charge);
      }

      //first create the full waveform
      std::vector<double> fullWaveform;
      _makeCrvWaveforms->MakeWaveform(times, charges, fullWaveform, startTime, _digitizationPeriod);
      _makeCrvWaveforms->AddElectronicNoise(fullWaveform, _noise, _randGaussQ);

      //break the waveform apart into short pieces (CrvDigiMC::NSamples)
      //and apply the zero suppression, i.e. set all waveform digi points to zero which are below the minimum voltage,
      //if the neighboring digi points are also below the minimum voltage
      for(size_t i=0; i<fullWaveform.size(); i++)
      {
        if(SingleWaveformStart(fullWaveform, i)) //acts as a zero suppression
        {
          //start new single waveform
          double digiStartTime=startTime+i*_digitizationPeriod;

          //collect voltages
          std::array<double,CrvDigiMC::NSamples> voltages;
          for(size_t singleWaveformIndex=0; singleWaveformIndex<CrvDigiMC::NSamples; i++, singleWaveformIndex++)
          {
            if(i<fullWaveform.size()) voltages[singleWaveformIndex]=fullWaveform[i];
            else voltages[singleWaveformIndex]=0.0;  //so that all unused single waveform samples are set to zero
          }

          //collect CrvSteps and SimParticles responsible for this single waveform
          std::set<art::Ptr<CrvStep> > steps;  //use a set to remove dublicate steppoints
          std::map<art::Ptr<SimParticle>, int> simparticles;
          for(size_t j=0; j<timesAndCharges.size(); j++)
          {
            if(timesAndCharges[j]._time>=digiStartTime-_singlePEWaveformMaxTime && 
               timesAndCharges[j]._time<=digiStartTime+CrvDigiMC::NSamples*_digitizationPeriod)
            {
              steps.insert(timesAndCharges[j]._step);
              if(timesAndCharges[j]._step.isNonnull()) simparticles[timesAndCharges[j]._step->simParticle()]++;
            }
          }

          //loop through the steps to fill the single waveform
          std::vector<art::Ptr<CrvStep> > stepVector;
          std::set<art::Ptr<CrvStep> >::iterator stepIter;
          for(stepIter=steps.begin(); stepIter!=steps.end(); stepIter++) stepVector.push_back(*stepIter);

          //find the most likely SimParticle
          //if no SimParticle was recorded for this single waveform, then it was caused either by noise hits (if the threshold is low enough),
          //or is the tail end of the peak. in that case, _simparticle will be null (set by the default constructor of art::Ptr)
          art::Ptr<SimParticle> simParticle;
          std::map<art::Ptr<SimParticle>,int >::iterator simparticleIter;
          int simparticleCount=0;
          for(simparticleIter=simparticles.begin(); simparticleIter!=simparticles.end(); simparticleIter++)
          {
            if(simparticleIter->second>simparticleCount)
            {
              simparticleCount=simparticleIter->second;
              simParticle=simparticleIter->first;
            }
          }

          i--;
          crvDigiMCCollection->emplace_back(voltages, stepVector, simParticle, digiStartTime, barIndex, SiPM);
        }
      }
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
