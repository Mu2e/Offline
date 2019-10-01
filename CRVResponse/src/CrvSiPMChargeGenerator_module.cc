//
// A module to create CRV SiPM charges from CRV photons
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
//
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvSiPMCharges.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvPhotonsCollection.hh"
#include "MCDataProducts/inc/CrvSiPMChargesCollection.hh"
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
    double      _blindTime;             //time window during which the SiPM is blind
    double      _microBunchPeriod;

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
    _blindTime(pset.get<double>("blindTime")),                       //500ns
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
    mu2e::ConditionsHandle<mu2e::AcceleratorParams> accPar("ignored");
    _microBunchPeriod = accPar->deBuncherPeriod;
    _makeCrvSiPMCharges->SetSiPMConstants(_nPixelsX, _nPixelsY, _overvoltage, _blindTime, _microBunchPeriod,
                                            _timeConstant, _capacitance, _probabilities, _inactivePixels);
  }

  void CrvSiPMChargeGenerator::produce(art::Event& event)
  {
    std::unique_ptr<CrvSiPMChargesCollection> crvSiPMChargesCollection(new CrvSiPMChargesCollection);

    art::Handle<CrvPhotonsCollection> crvPhotonsCollection;
    event.getByLabel(_crvPhotonsModuleLabel,"",crvPhotonsCollection);

    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator iter;
    for(iter=counters.begin(); iter!=counters.end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = (*iter)->index();
      CrvPhotonsCollection::const_iterator crvPhotons=crvPhotonsCollection->find(barIndex);

      CrvSiPMCharges &crvSiPMCharges = (*crvSiPMChargesCollection)[barIndex];

      for(int SiPM=0; SiPM<4; SiPM++)
      {

        if(!(*iter)->getBarDetail().hasCMB(SiPM%2)) continue;  //no SiPM charges at non-existing SiPMs
                                                               //SiPM%2 returns the side of the CRV counter
                                                               //0 ... negative side
                                                               //1 ... positive side

        if(_randFlat.fire() < _deadSiPMProbability) continue;  //assume that this random SiPM is dead

        std::vector<std::pair<double,size_t> > photonTimesAdjusted;   //pair of photon time and index in the original photon vector
        if(crvPhotons!=crvPhotonsCollection->end())  //if there are no photons at this SiPM, then we still need to continue to simulate dark noise
        {
          const std::vector<CrvPhotons::SinglePhoton> &photonTimes = crvPhotons->second.GetPhotons(SiPM);
          for(size_t iphoton=0; iphoton<photonTimes.size(); iphoton++)
          {
            double time = photonTimes[iphoton]._time;
            time = fmod(time,_microBunchPeriod);
            if(time>_blindTime) photonTimesAdjusted.push_back(std::pair<double,size_t>(time,iphoton)); //wrapped time
            //no ghost hits, since the SiPMs are off during the blind time (which is longer than the "ghost time")
          }
        }

        std::vector<mu2eCrv::SiPMresponse> SiPMresponseVector;
        _makeCrvSiPMCharges->Simulate(photonTimesAdjusted, SiPMresponseVector);

        std::vector<CrvSiPMCharges::CrvSingleCharge> &chargesOneSiPM = crvSiPMCharges.GetSiPMCharges(SiPM);

        std::vector<mu2eCrv::SiPMresponse>::const_iterator responseIter;
        for(responseIter=SiPMresponseVector.begin(); responseIter!=SiPMresponseVector.end(); responseIter++)
        {
          //time in SiPMresponseVector is between blindTime and microBunchPeriod
          //no additional time wrapping and check for blind time is required
          double time=responseIter->_time;
          double charge=responseIter->_charge;
          double chargeInPEs=responseIter->_chargeInPEs;
          int photonIndex=responseIter->_photonIndex;
          bool darkNoise=responseIter->_darkNoise;
          if(!darkNoise)
          {
            const std::vector<CrvPhotons::SinglePhoton> &photonTimes = crvPhotons->second.GetPhotons(SiPM);
            chargesOneSiPM.emplace_back(time, charge, chargeInPEs,photonTimes[photonIndex]._step);
          }
          else chargesOneSiPM.emplace_back(time, charge, chargeInPEs);
//std::cout<<"SiPM charge   bar index: "<<barIndex<<"   SiPM: "<<SiPM<<"   time: "<<time<<std::endl;
        }

      }

      //2 options:
      //(1) -create a crvSiPMCharges object as a reference to a crvSiPMChargesCollection map entry at the beginning for all counters
      //    -fill this this crvSiPMCharges object
      //    -if the crvSiPMCharges object stays empty, erase the map entry in crvSiPMChargesCollection
      //(2) -create a standalone crvSiPMCharges object
      //    -fill this this crvSiPMCharges object
      //    -if the crvSiPMCharges didn't stay empty, create a new map entry in crvSiPMChargesCollection and fill its content with
      //     the new crvSiPMCharges  <---- too time consuming, therefore use option (1)

      if(crvSiPMCharges.IsEmpty()) crvSiPMChargesCollection->erase(barIndex);
    }

    event.put(std::move(crvSiPMChargesCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvSiPMChargeGenerator;
DEFINE_ART_MODULE(CrvSiPMChargeGenerator)
