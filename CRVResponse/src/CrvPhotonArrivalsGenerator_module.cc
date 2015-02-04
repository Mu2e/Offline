//
// A module to create CRV photons arriving at the SiPMs (using StepPointMCs)
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvPhotonArrivals.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "MCDataProducts/inc/CrvPhotonArrivalsCollection.hh"

#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VModularPhysicsList.hh"
#include "G4ParticleTable.hh"

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

#include <TRandom3.h>
#include <TMath.h>

namespace mu2e 
{
  class CrvPhotonArrivalsGenerator : public art::EDProducer 
  {

    public:
    explicit CrvPhotonArrivalsGenerator(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void endJob();

    private:
    std::string _g4ModuleLabel;
    std::string _lookupTableFileName;

    boost::shared_ptr<MakeCrvPhotonArrivals> _makeCrvPhotonArrivals;

    double      _scintillationYield;
    double      _scintillatorDecayTimeFast;
    double      _scintillatorDecayTimeSlow;
    double      _fiberDecayTime;

    G4ParticleTable *_particleTable;
  };

  CrvPhotonArrivalsGenerator::CrvPhotonArrivalsGenerator(fhicl::ParameterSet const& pset) :
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
    _lookupTableFileName(pset.get<std::string>("lookupTableFileName")),
    _scintillationYield(pset.get<double>("scintillationYield",820.0)),    //820.0 photons per MeV
    _scintillatorDecayTimeFast(pset.get<double>("scintillatorDecayTimeFast",3.0)),  //3.0 ns
    _scintillatorDecayTimeSlow(pset.get<double>("scintillatorDecayTimeSlow",10.0)), //10.0 ns
    _fiberDecayTime(pset.get<double>("fiberDecayTime",7.4))     //7.4 ns
  {
//    ConfigFileLookupPolicy configFile;
//    _lookupTableFileName = configFile(_lookupTableFileName);
    _makeCrvPhotonArrivals = boost::shared_ptr<MakeCrvPhotonArrivals>(new MakeCrvPhotonArrivals());
    _makeCrvPhotonArrivals->LoadLookupTable(_lookupTableFileName.c_str());
    _makeCrvPhotonArrivals->SetScintillationYield(_scintillationYield);
    _makeCrvPhotonArrivals->SetScintillatorDecayTimeFast(_scintillatorDecayTimeFast);
    _makeCrvPhotonArrivals->SetScintillatorDecayTimeSlow(_scintillatorDecayTimeSlow);
    _makeCrvPhotonArrivals->SetFiberDecayTime(_fiberDecayTime);
    produces<CrvPhotonArrivalsCollection>();

    std::string physName = "QGSP_BERT_EMV";
    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = factory.GetReferencePhysList(physName);
    if(!phys) throw std::logic_error("can't find physics list "+physName);
    for (int i=0; ; i++) 
    {
       G4VPhysicsConstructor* elem = const_cast<G4VPhysicsConstructor*>(phys->GetPhysics(i));
       if (elem == NULL) break;
       elem->ConstructParticle();;
    }
    _particleTable = G4ParticleTable::GetParticleTable();
    _particleTable->SetReadiness();
  }

  void CrvPhotonArrivalsGenerator::beginJob()
  {
  }

  void CrvPhotonArrivalsGenerator::endJob()
  {
  }

  void CrvPhotonArrivalsGenerator::produce(art::Event& event) 
  {
    std::unique_ptr<CrvPhotonArrivalsCollection> crvPhotonArrivalsCollection(new CrvPhotonArrivalsCollection);

    GeomHandle<CosmicRayShield> CRS;
    StepInstanceName CRVInstance(StepInstanceName::CRV);
    art::Handle<StepPointMCCollection> CRVSteps;
    event.getByLabel(_g4ModuleLabel,CRVInstance.name(),CRVSteps);

    for(StepPointMCCollection::const_iterator iter=CRVSteps->begin(); iter!=CRVSteps->end(); iter++)
    {
      StepPointMC const& step(*iter);

      const CLHEP::Hep3Vector &p1 = step.position();
      CLHEP::Hep3Vector p2 = p1 + step.momentum().unit()*step.stepLength();
      double energyDepositedTotal= step.totalEDep();
      double energyDepositedNonIonizing = step.nonIonizingEDep();

      int PDGcode = step.simParticle()->pdgId();
      if(_particleTable->FindParticle(PDGcode)==NULL)
      {
        std::cerr<<"Error in CrvPhotonArrivalsGenerator: Found a PDG code which is not in the GEANT particle table: ";
        std::cerr<<PDGcode<<std::endl;
        std::cerr<<"FIXME: Skipping this StepPoint."<<std::endl;
        continue;
      }
      double mass = _particleTable->FindParticle(PDGcode)->GetPDGMass();  //MeV/c^2
      double charge = _particleTable->FindParticle(PDGcode)->GetPDGCharge(); 

      double momentum1 = step.momentum().mag(); //MeV/c
      double energy1 = sqrt(momentum1*momentum1 + mass*mass); //MeV
//FIXME: does not take the energy of daughter particles into account
      double energy2 = energy1 - energyDepositedTotal; //MeV  

      double gamma1 = energy1 / mass;
      double gamma2 = energy2 / mass;
      double beta1 = sqrt(1.0-1.0/(gamma1*gamma1));
      double beta2 = sqrt(1.0-1.0/(gamma2*gamma2));
      double beta = (beta1+beta2)/2.0;
      double velocity = beta*CLHEP::c_light;
      double t1 = step.time();
      double t2 = t1 + step.stepLength()/velocity;

//if there is a following step point, it will give a more realistic energy and time
      StepPointMCCollection::const_iterator iterNextStep = iter;
      iterNextStep++;
      if(iterNextStep!=CRVSteps->end())
      {
        StepPointMC const& nextStep(*iterNextStep);
        if(nextStep.barIndex()==step.barIndex() && nextStep.simParticle()->id()==step.simParticle()->id())
        {
          p2 = nextStep.position();
          double momentum2 = nextStep.momentum().mag(); //MeV/c
          energy2 = sqrt(momentum2*momentum2 + mass*mass); //MeV
          gamma2 = energy2 / mass;
          beta2 = sqrt(1.0-1.0/(gamma2*gamma2));
          beta = (beta1+beta2)/2.0;
          velocity = beta*CLHEP::c_light;
          t2 = nextStep.time();
        }
      }

      const CRSScintillatorBar &CRSbar = CRS->getBar(step.barIndex());
      CLHEP::Hep3Vector p1Local = CRSbar.toLocal(p1);
      CLHEP::Hep3Vector p2Local = CRSbar.toLocal(p2);

      _makeCrvPhotonArrivals->SetActualHalfLength(CRSbar.getHalfLength());
      _makeCrvPhotonArrivals->MakePhotons(p1Local, p2Local, t1, t2,  
                              PDGcode, beta, charge,
                              energyDepositedTotal,
                              energyDepositedNonIonizing);

      bool needToStore = false;
      if(iterNextStep==CRVSteps->end()) needToStore=true;
      else
      {
        if(iterNextStep->barIndex()!=step.barIndex()) needToStore=true;
      }

      if(needToStore)
      {
        CrvPhotonArrivals &crvPhotons = (*crvPhotonArrivalsCollection)[step.barIndex()];
        for(int SiPM=0; SiPM<4; SiPM++)
        {
          const std::vector<double> &times=_makeCrvPhotonArrivals->GetArrivalTimes(SiPM);
          crvPhotons.GetPhotonArrivalTimes(SiPM).insert(crvPhotons.GetPhotonArrivalTimes(SiPM).begin(),times.begin(),times.end());
        }
        _makeCrvPhotonArrivals->Reset();
      }
    }

    event.put(std::move(crvPhotonArrivalsCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvPhotonArrivalsGenerator;
DEFINE_ART_MODULE(CrvPhotonArrivalsGenerator)
