//
// A module to create CRV waveforms from CRV SiPMCharges
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CRVResponse/inc/MakeCrvWaveforms.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CRVConditions/inc/CRVDigitizationPeriod.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/MCDataProducts/inc/CrvSiPMCharges.hh"
#include "Offline/MCDataProducts/inc/CrvDigiMC.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/SeedService/inc/SeedService.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
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
    std::string _eventWindowMarkerLabel;
    std::string _protonBunchTimeMCLabel;

    boost::shared_ptr<mu2eCrv::MakeCrvWaveforms> _makeCrvWaveforms;

    double                              _digitizationStart, _digitizationEnd;
    double                              _FEBtimeSpread;
    double                              _minVoltage;
    double                              _noise;
    double                              _singlePEWaveformMaxTime;

    CLHEP::HepRandomEngine&             _engine;
    CLHEP::RandFlat                     _randFlat;
    CLHEP::RandGaussQ                   _randGaussQ;


    std::vector<double> _digitizationPointShiftFEBsSide0, _digitizationPointShiftFEBsSide1;
    std::vector<double> _timeShiftFEBsSide0, _timeShiftFEBsSide1;

    struct ChargeCluster
    {
      std::vector<std::pair<double,double> > timesAndCharges;
    };

    void FindChargeClusters(const std::vector<CrvSiPMCharges::SingleCharge> &timesAndCharges,
                            std::vector<ChargeCluster> &chargeClusters,
                            double timeShiftFEB);
    bool SingleWaveformStart(std::vector<double> &fullWaveform, size_t i);
  };

  CrvWaveformsGenerator::CrvWaveformsGenerator(fhicl::ParameterSet const& pset) :
    EDProducer{pset},
    _crvSiPMChargesModuleLabel(pset.get<std::string>("crvSiPMChargesModuleLabel")),
    _singlePEWaveformFileName(pset.get<std::string>("singlePEWaveformFileName")),
    _eventWindowMarkerLabel(pset.get<std::string>("eventWindowMarker","EWMProducer")),
    _protonBunchTimeMCLabel(pset.get<std::string>("protonBunchTimeMC","EWMProducer")),
    _digitizationStart(pset.get<double>("digitizationStart")),       //400ns
    _digitizationEnd(pset.get<double>("digitizationEnd")),           //1750ns
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
  }

  void CrvWaveformsGenerator::produce(art::Event& event)
  {
    art::Handle<EventWindowMarker> eventWindowMarker;
    event.getByLabel(_eventWindowMarkerLabel,"",eventWindowMarker);
    EventWindowMarker::SpillType spillType = eventWindowMarker->spillType();

    art::Handle<ProtonBunchTimeMC> protonBunchTimeMC;
    event.getByLabel(_protonBunchTimeMCLabel, protonBunchTimeMC);
    double TDC0time = -protonBunchTimeMC->pbtime_;

    ProditionsHandle<EventTiming> eventTimingHandle;
    const EventTiming &eventTiming = eventTimingHandle.get(event.id());
    double jitter = TDC0time - eventTiming.timeFromProtonsToDRMarker();

    double digitizationStart=_digitizationStart+jitter;
    double digitizationEnd=_digitizationEnd+jitter;
    if(spillType!=EventWindowMarker::SpillType::onspill)
    {
      double eventWindowLength = eventWindowMarker->eventLength();
      digitizationStart = TDC0time;
      digitizationEnd = digitizationStart + eventWindowLength;
    }

    std::unique_ptr<CrvDigiMCCollection> crvDigiMCCollection(new CrvDigiMCCollection);

    art::Handle<CrvSiPMChargesCollection> crvSiPMChargesCollection;
    event.getByLabel(_crvSiPMChargesModuleLabel,"",crvSiPMChargesCollection);

    GeomHandle<CosmicRayShield> CRS;
    _digitizationPointShiftFEBsSide0.clear();
    _digitizationPointShiftFEBsSide1.clear();
    _timeShiftFEBsSide0.clear();
    _timeShiftFEBsSide1.clear();
    unsigned int nCounters = CRS->getAllCRSScintillatorBars().size();
    unsigned int nFEBs = rint(ceil(nCounters/32.0));
    for(unsigned int i=0; i<nFEBs; ++i)
    {
      //the closest digitization point with respect to a certain time
      //can happen anywhere within the digitization period of 12.55ns.
      //this time is different for each FEB.
      _digitizationPointShiftFEBsSide0.emplace_back(_randFlat.fire()*CRVDigitizationPeriod);
      _digitizationPointShiftFEBsSide1.emplace_back(_randFlat.fire()*CRVDigitizationPeriod);
      //the FEBs will be synchronized to account for cable length differences etc.,
      //but there may still be small time differences between the FEBs.
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

      double timeShiftFEB=0;
      if(SiPM%2==0 && FEB<_timeShiftFEBsSide0.size()) timeShiftFEB=_timeShiftFEBsSide0[FEB];
      if(SiPM%2==1 && FEB<_timeShiftFEBsSide1.size()) timeShiftFEB=_timeShiftFEBsSide1[FEB];

      const std::vector<CrvSiPMCharges::SingleCharge> &timesAndCharges = iter->GetCharges();
      std::vector<ChargeCluster> chargeClusters;
      FindChargeClusters(timesAndCharges, chargeClusters, timeShiftFEB);

      //need to find where this FEB's TDC=0 is located with respect to the global time
      //can be anywhere within the digitization period
      double digitizationPointShiftFEB=0;
      if(SiPM%2==0 && FEB<_digitizationPointShiftFEBsSide0.size()) digitizationPointShiftFEB=_digitizationPointShiftFEBsSide0[FEB];
      if(SiPM%2==1 && FEB<_digitizationPointShiftFEBsSide1.size()) digitizationPointShiftFEB=_digitizationPointShiftFEBsSide1[FEB];
      double TDC0timeAdjusted=TDC0time+digitizationPointShiftFEB;  //that's the time when TDC=0 for this FEB

      for(size_t iCluster=0; iCluster<chargeClusters.size(); ++iCluster)
      {
        //if the number of charges in this cluster cannot achieve the minimum voltage, skip this cluster
        if(chargeClusters[iCluster].timesAndCharges.size()*_makeCrvWaveforms->GetSinglePEMaxVoltage()<_minVoltage) continue;

        //find the TDC time when the first charge occurs (adjusted for this FEB)
        double firstChargeTime=chargeClusters[iCluster].timesAndCharges.front().first;
        firstChargeTime-=1.0*CRVDigitizationPeriod;
        double TDCstartTimeAdjusted=ceil((firstChargeTime-TDC0timeAdjusted)/CRVDigitizationPeriod) * CRVDigitizationPeriod + TDC0timeAdjusted;

        //first create the full waveform
        std::vector<double> fullWaveform;
        _makeCrvWaveforms->MakeWaveform(chargeClusters[iCluster].timesAndCharges,
                                        fullWaveform, TDCstartTimeAdjusted, CRVDigitizationPeriod);
        _makeCrvWaveforms->AddElectronicNoise(fullWaveform, _noise, _randGaussQ);

        //break the waveform apart into short pieces (CrvDigiMC::NSamples) and apply the zero suppression
        //don't digitize outside of digitizationStart and digitizationEnd
        for(size_t i=0; i<fullWaveform.size(); ++i)
        {
          if(SingleWaveformStart(fullWaveform, i)) //acts as a zero suppression
          {
            //start new single waveform
            double digiStartTime=TDCstartTimeAdjusted+i*CRVDigitizationPeriod;
            if(digiStartTime<digitizationStart) continue; //digis cannot start before the digitization interval
            if(digiStartTime>digitizationEnd) continue; //digis cannot start after the digitization interval
//            if(digiStartTime+(CrvDigiMC::NSamples-1)*CRVDigitizationPeriod>digitizationEnd) continue; //digis cannot end after the digitization interval

            //collect voltages
            std::array<double,CrvDigiMC::NSamples> voltages;
            for(size_t singleWaveformIndex=0; singleWaveformIndex<CrvDigiMC::NSamples; ++i, ++singleWaveformIndex)
            {
              if(i<fullWaveform.size() && TDCstartTimeAdjusted+i*CRVDigitizationPeriod<=digitizationEnd) voltages[singleWaveformIndex]=fullWaveform[i];
              else voltages[singleWaveformIndex]=0.0;  //so that all unused single waveform samples are set to zero
            }

            //collect CrvSteps and SimParticles responsible for this single waveform
            std::set<art::Ptr<CrvStep> > steps;  //use a set to remove dublicate steppoints
            std::map<art::Ptr<SimParticle>, int> simparticles;
            for(size_t j=0; j<timesAndCharges.size(); ++j)
            {
              if(timesAndCharges[j]._time>=digiStartTime-_singlePEWaveformMaxTime &&
                 timesAndCharges[j]._time<=digiStartTime+CrvDigiMC::NSamples*CRVDigitizationPeriod)
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
            for(simparticleIter=simparticles.begin(); simparticleIter!=simparticles.end(); ++simparticleIter)
            {
              if(simparticleIter->second>simparticleCount)
              {
                simparticleCount=simparticleIter->second;
                simParticle=simparticleIter->first;
              }
            }

            --i;
            crvDigiMCCollection->emplace_back(voltages, stepVector, simParticle, digiStartTime, TDC0timeAdjusted, barIndex, SiPM);
          }
        }
      }
    }

    event.put(std::move(crvDigiMCCollection));
  } // end produce

  void CrvWaveformsGenerator::FindChargeClusters(const std::vector<CrvSiPMCharges::SingleCharge> &timesAndCharges,
                                                   std::vector<ChargeCluster> &chargeClusters,
                                                   double timeShiftFEB)
  {
    chargeClusters.reserve(timesAndCharges.size());
    for(size_t i=0; i<timesAndCharges.size(); ++i)
    {
      //No check whether times are within digitizationStart-_digitizationMargin and digitizationEnd
      double timeTmp=timesAndCharges[i]._time+timeShiftFEB;  //apply timeShift to account for inaccuraries in the FEB time calibration

      if(chargeClusters.empty())
      {
        chargeClusters.resize(1);
        chargeClusters.back().timesAndCharges.reserve(timesAndCharges.size());
      }
      else
      {
        if(timeTmp-chargeClusters.back().timesAndCharges.back().first>_singlePEWaveformMaxTime+4.0*CRVDigitizationPeriod)
        {
          chargeClusters.resize(chargeClusters.size()+1);
          chargeClusters.back().timesAndCharges.reserve(timesAndCharges.size());
        }
      }
      chargeClusters.back().timesAndCharges.emplace_back(timeTmp,timesAndCharges[i]._charge);
    }
  }

  bool CrvWaveformsGenerator::SingleWaveformStart(std::vector<double> &fullWaveform, size_t i)
  {
    //test up to 2 points before and after i, if the waveform is/was above _minVoltage
    //to decide whether to start or continue to record digis
    for(size_t j=std::max<size_t>(i,2)-2; j<std::min<size_t>(i+3,fullWaveform.size()); ++j)
    {
      if(fullWaveform[j]>_minVoltage) return true;
    }

    return false;
  }

} // end namespace mu2e

using mu2e::CrvWaveformsGenerator;
DEFINE_ART_MODULE(CrvWaveformsGenerator)
