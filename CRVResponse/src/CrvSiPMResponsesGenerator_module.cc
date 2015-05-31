//
// A module to create CRV SiPM responses from CRV photon arrivals
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvSiPMResponses.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvPhotonArrivalsCollection.hh"
#include "MCDataProducts/inc/CrvSiPMResponsesCollection.hh"
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
  class CrvSiPMResponsesGenerator : public art::EDProducer 
  {

    public:
    explicit CrvSiPMResponsesGenerator(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    std::string _crvPhotonArrivalsModuleLabel;
    int         _numberPixels;
    double      _bias;
    double      _scaleFactor;
    double      _minCharge;
    double      _blindTime;             //time during which the SiPM is blind
    double      _microBunchPeriod;
    MakeCrvSiPMResponses::ProbabilitiesStruct _probabilities;

    boost::shared_ptr<MakeCrvSiPMResponses> _makeCrvSiPMResponses;

    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;
  };

  CrvSiPMResponsesGenerator::CrvSiPMResponsesGenerator(fhicl::ParameterSet const& pset) :
    _crvPhotonArrivalsModuleLabel(pset.get<std::string>("crvPhotonArrivalsModuleLabel")),
    _numberPixels(pset.get<int>("numberPixels")),   //1600
    _bias(pset.get<double>("bias")),                //2.5V
    _scaleFactor(pset.get<double>("scaleFactor")),  //0.08 (based on a time step of 1.0ns)
    _minCharge(pset.get<double>("minCharge")),      //3.0PE
    _blindTime(pset.get<double>("blindTime")),      //500ns
    _randFlat(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _randPoissonQ(art::ServiceHandle<art::RandomNumberGenerator>()->getEngine())
  {
    produces<CrvSiPMResponsesCollection>();
    _probabilities._constGeigerProbCoef = pset.get<double>("GeigerProbCoef");  //2.0
    _probabilities._constGeigerProbVoltScale = pset.get<double>("GeigerProbVoltScale");  //3.0
    _probabilities._constTrapType0Prob = pset.get<double>("TrapType0Prob");  //trap_prob*trap_type0_prob=0.2*0.7=0.14
    _probabilities._constTrapType1Prob = pset.get<double>("TrapType1Prob");  //trap_prob*trap_type1_prob=0.2*0.3=0.06
    _probabilities._constTrapType0Lifetime = pset.get<double>("TrapType0Lifetime");   //5.0ns
    _probabilities._constTrapType1Lifetime = pset.get<double>("TrapType1Lifetime");  //50.0ns
    _probabilities._constThermalProb = pset.get<double>("ThermalProb"); //6.25e-7ns^1     1MHz at SiPM --> 1MHz/#pixel=625Hz at Pixel --> 625 s^-1 = 6.25e-7 ns^-1   //exp(-E_th/T)=1.6e-6
    _probabilities._constPhotonProduction = pset.get<double>("PhotonProduction");  //0.1  //0.4
  }

  void CrvSiPMResponsesGenerator::beginJob()
  {
  }

  void CrvSiPMResponsesGenerator::beginRun(art::Run &run)
  {
    mu2e::ConditionsHandle<mu2e::AcceleratorParams> accPar("ignored");
    _microBunchPeriod = accPar->deBuncherPeriod;
    _makeCrvSiPMResponses = boost::shared_ptr<MakeCrvSiPMResponses>(new MakeCrvSiPMResponses(_randFlat, _randPoissonQ));
    _makeCrvSiPMResponses->SetSiPMConstants(_numberPixels, _bias, _blindTime, _microBunchPeriod, 
                                            _scaleFactor, _probabilities);
  }

  void CrvSiPMResponsesGenerator::endJob()
  {
  }

  void CrvSiPMResponsesGenerator::produce(art::Event& event) 
  {
    std::unique_ptr<CrvSiPMResponsesCollection> crvSiPMResponsesCollection(new CrvSiPMResponsesCollection);

    art::Handle<CrvPhotonArrivalsCollection> crvPhotonArrivalsCollection;
    event.getByLabel(_crvPhotonArrivalsModuleLabel,"",crvPhotonArrivalsCollection);

    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator iter; 
    for(iter=counters.begin(); iter!=counters.end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = (*iter)->index();
      CrvPhotonArrivalsCollection::const_iterator crvPhotons=crvPhotonArrivalsCollection->find(barIndex); 

      CrvSiPMResponses crvSiPMResponses;
      bool minChargeReached = false;

      for(int SiPM=0; SiPM<4; SiPM++) 
      {
        std::vector<CrvSiPMResponses::CrvSingleSiPMResponse> &responsesOneSiPM = crvSiPMResponses.GetSiPMResponses(SiPM);

        std::vector<double> photonArrivalTimesAdjusted;
        if(crvPhotons!=crvPhotonArrivalsCollection->end())
        {
          const std::vector<double> &photonArrivalTimes = crvPhotons->second.GetPhotonArrivalTimes(SiPM);
          std::vector<double>::const_iterator timeIter;
          for(timeIter=photonArrivalTimes.begin(); timeIter!=photonArrivalTimes.end(); timeIter++)
          {
            double time = *timeIter;
	    time = fmod(time,_microBunchPeriod);  //"ghost hits" are not used / "splitting of hits" is acceptable, 
                                                  //since we are blind outside the window
            if(time>_blindTime) photonArrivalTimesAdjusted.push_back(time);
          }
        }

        std::vector<SiPMresponse> SiPMresponseVector;
        _makeCrvSiPMResponses->Simulate(photonArrivalTimesAdjusted, SiPMresponseVector);

        double totalCharge=0;
        std::vector<SiPMresponse>::const_iterator responseIter;
        for(responseIter=SiPMresponseVector.begin(); responseIter!=SiPMresponseVector.end(); responseIter++)
        {
          responsesOneSiPM.emplace_back(responseIter->_time, responseIter->_charge);
          totalCharge+=responseIter->_charge;
        }
        if(totalCharge>=_minCharge) minChargeReached=true;
      }

      if(minChargeReached)
      {
        (*crvSiPMResponsesCollection)[barIndex] = crvSiPMResponses;
      }
    }

    event.put(std::move(crvSiPMResponsesCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvSiPMResponsesGenerator;
DEFINE_ART_MODULE(CrvSiPMResponsesGenerator)
