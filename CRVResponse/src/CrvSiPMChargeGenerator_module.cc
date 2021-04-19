//
// A module to create CRV SiPM charges from CRV photons
//
//
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvSiPMCharges.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "DAQConditions/inc/EventTiming.hh"
#include "DataProducts/inc/EventWindowMarker.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvPhotons.hh"
#include "MCDataProducts/inc/CrvSiPMCharges.hh"
#include "MCDataProducts/inc/ProtonBunchTimeMC.hh"
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
  class CrvSiPMChargeGenerator : public art::EDProducer
  {

    public:
    explicit CrvSiPMChargeGenerator(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginRun(art::Run &run);

    private:
    std::string _crvPhotonsModuleLabel;
    double      _deadSiPMProbability;
    int         _nPixelsX;
    int         _nPixelsY;
    double      _overvoltage;
    double      _timeConstant;
    double      _capacitance;
    double      _digitizationStart, _digitizationEnd;
    std::string _eventWindowMarkerLabel;
    std::string _protonBunchTimeMCLabel;

    mu2eCrv::MakeCrvSiPMCharges::ProbabilitiesStruct _probabilities;
    std::vector<std::pair<int,int> >   _inactivePixels;

    boost::shared_ptr<mu2eCrv::MakeCrvSiPMCharges> _makeCrvSiPMCharges;

    CLHEP::HepRandomEngine& _engine;
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;

    std::string             _photonMapFileName;
    ConfigFileLookupPolicy  _resolveFullPath;
  };

  CrvSiPMChargeGenerator::CrvSiPMChargeGenerator(fhicl::ParameterSet const& pset) :
    EDProducer{pset},
    _crvPhotonsModuleLabel(pset.get<std::string>("crvPhotonsModuleLabel")),
    _deadSiPMProbability(pset.get<double>("deadSiPMProbability")),   //0.01
    _nPixelsX(pset.get<int>("nPixelsX")),                            //40
    _nPixelsY(pset.get<int>("nPixelsY")),                            //40
    _overvoltage(pset.get<double>("overvoltage")),                   //3.0V
    _timeConstant(pset.get<double>("timeConstant")),                 //12.0ns
    _capacitance(pset.get<double>("capacitance")),                   //8.84e-14F (per pixel)
    _digitizationStart(pset.get<double>("digitizationStart")),       //400ns
    _digitizationEnd(pset.get<double>("digitizationEnd")),           //1750ns
    _eventWindowMarkerLabel(pset.get<std::string>("eventWindowMarker","EWMProducer")),
    _protonBunchTimeMCLabel(pset.get<std::string>("protonBunchTimeMC","EWMProducer")),
    _inactivePixels(pset.get<std::vector<std::pair<int,int> > >("inactivePixels")),      //{18,18},....,{21,21}
    _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())},
    _randFlat{_engine},
    _randPoissonQ{_engine},
    _photonMapFileName(pset.get<std::string>("photonMapFileName"))
  {
    produces<CrvSiPMChargesCollection>();
    _probabilities._avalancheProbParam1 = pset.get<double>("AvalancheProbParam1");  //0.65
    _probabilities._avalancheProbParam2 = pset.get<double>("AvalancheProbParam2");  //2.7
    _probabilities._trapType0Prob = pset.get<double>("TrapType0Prob");              //0
    _probabilities._trapType1Prob = pset.get<double>("TrapType1Prob");              //0
    _probabilities._trapType0Lifetime = pset.get<double>("TrapType0Lifetime");      //5.0ns
    _probabilities._trapType1Lifetime = pset.get<double>("TrapType1Lifetime");      //50.0ns
    _probabilities._thermalRate = pset.get<double>("ThermalRate");                  //3.0e-4 ns^-1   300MHz for entire SiPM
    _probabilities._crossTalkProb = pset.get<double>("CrossTalkProb");              //0.05

    std::string fullPhotonMapFileName(_resolveFullPath(_photonMapFileName));
    _makeCrvSiPMCharges = boost::shared_ptr<mu2eCrv::MakeCrvSiPMCharges>(new mu2eCrv::MakeCrvSiPMCharges(_randFlat, _randPoissonQ, fullPhotonMapFileName));
  }

  void CrvSiPMChargeGenerator::beginRun(art::Run &run)
  {
    _makeCrvSiPMCharges->SetSiPMConstants(_nPixelsX, _nPixelsY, _overvoltage, _timeConstant, _capacitance, _probabilities, _inactivePixels);
  }

  void CrvSiPMChargeGenerator::produce(art::Event& event)
  {
    std::unique_ptr<CrvSiPMChargesCollection> crvSiPMChargesCollection(new CrvSiPMChargesCollection);

    art::Handle<CrvPhotonsCollection> crvPhotonsCollection;
    event.getByLabel(_crvPhotonsModuleLabel,"",crvPhotonsCollection);

    art::Handle<EventWindowMarker> eventWindowMarker;
    event.getByLabel(_eventWindowMarkerLabel,"",eventWindowMarker);
    EventWindowMarker::SpillType spillType = eventWindowMarker->spillType();

    art::Handle<ProtonBunchTimeMC> protonBunchTimeMC;
    event.getByLabel(_protonBunchTimeMCLabel, protonBunchTimeMC);
    double TDC0time = -protonBunchTimeMC->pbtime_;

    ProditionsHandle<EventTiming> eventTimingHandle;
    const EventTiming &eventTiming = eventTimingHandle.get(event.id());
    double jitter = TDC0time - eventTiming.timeFromProtonsToDRMarker();

    double startTime=_digitizationStart+jitter;
    double endTime=_digitizationEnd+jitter;
    if(spillType!=EventWindowMarker::SpillType::onspill)
    {
      double eventWindowLength = eventWindowMarker->eventLength();
      startTime = TDC0time;
      endTime = startTime + eventWindowLength;
    }

    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator iter;
    for(iter=counters.begin(); iter!=counters.end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = (*iter)->index();

      for(int SiPM=0; SiPM<4; SiPM++)
      {
        if(!(*iter)->getBarDetail().hasCMB(SiPM%2)) continue;  //no SiPM charges at non-existing SiPMs
                                                               //SiPM%2 returns the side of the CRV counter
                                                               //0 ... negative side
                                                               //1 ... positive side

        if(_randFlat.fire() < _deadSiPMProbability) continue;  //assume that this random SiPM is dead

        //time wrapping happened in the photon generator
        std::vector<std::pair<double,size_t> > photonTimesNew;   //pair of photon time and index in the original photon vector
        CrvPhotonsCollection::const_iterator crvPhotons;
        for(crvPhotons=crvPhotonsCollection->begin(); crvPhotons!=crvPhotonsCollection->end(); crvPhotons++)
        {
          if(crvPhotons->GetScintillatorBarIndex()==barIndex && crvPhotons->GetSiPMNumber()==SiPM)
          {
            const std::vector<CrvPhotons::SinglePhoton> &photonTimes = crvPhotons->GetPhotons();
            for(size_t iphoton=0; iphoton<photonTimes.size(); iphoton++)
            {
              double time = photonTimes[iphoton]._time;
              photonTimesNew.emplace_back(time,iphoton);
            }
            break;
          }
        }

        std::vector<mu2eCrv::SiPMresponse> SiPMresponseVector;
        _makeCrvSiPMCharges->Simulate(photonTimesNew, SiPMresponseVector, startTime, endTime);

        if(SiPMresponseVector.size()>0)
        {
          crvSiPMChargesCollection->emplace_back(barIndex,SiPM);
          std::vector<CrvSiPMCharges::SingleCharge> &charges = crvSiPMChargesCollection->back().GetCharges();

          std::vector<mu2eCrv::SiPMresponse>::const_iterator responseIter;
          for(responseIter=SiPMresponseVector.begin(); responseIter!=SiPMresponseVector.end(); responseIter++)
          {
            double time=responseIter->_time;
            double charge=responseIter->_charge;
            double chargeInPEs=responseIter->_chargeInPEs;
            int photonIndex=responseIter->_photonIndex;
            bool darkNoise=responseIter->_darkNoise;
            if(!darkNoise)
            {
              const std::vector<CrvPhotons::SinglePhoton> &photonTimes = crvPhotons->GetPhotons();
              charges.emplace_back(time, charge, chargeInPEs, photonTimes[photonIndex]._step);
            }
            else charges.emplace_back(time, charge, chargeInPEs);
          }
        }//non-empty SiPM charges
      }//SiPM
    }//barIndex

    event.put(std::move(crvSiPMChargesCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvSiPMChargeGenerator;
DEFINE_ART_MODULE(CrvSiPMChargeGenerator)
