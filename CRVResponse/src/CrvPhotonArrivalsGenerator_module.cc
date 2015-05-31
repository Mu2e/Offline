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

#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CrvPhotonArrivalsCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
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
  class CrvPhotonArrivalsGenerator : public art::EDProducer 
  {

    public:
    explicit CrvPhotonArrivalsGenerator(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void endJob();
    void beginRun(art::Run& r);

    private:
    std::vector<std::string> _g4ModuleLabels;
    std::vector<std::string> _processNames;

    std::vector<std::string>                                    _lookupTableFileNames;
    std::vector<double>                                         _lookupTableCounterLengths;
    std::map<double, boost::shared_ptr<MakeCrvPhotonArrivals> > _makeCrvPhotonArrivals;

    double      _scintillationYield;
    double      _scintillatorDecayTimeFast;
    double      _scintillatorDecayTimeSlow;
    double      _fiberDecayTime;

    double      _startTime;             //StepPoint times before this time will be ignored to reduce computation times
                                        //(in particular by ignoring hits during the beam flash).
                                        //This time should be at least 100ns before the end of the SiPM's blind time
                                        //to account for the travel time of the photons inside the CRV bar.
                                        //Default is 0.
    SimParticleTimeOffset _timeOffsets;

    CLHEP::RandFlat       _randFlat;
  };

  CrvPhotonArrivalsGenerator::CrvPhotonArrivalsGenerator(fhicl::ParameterSet const& pset) :
    _g4ModuleLabels(pset.get<std::vector<std::string> >("g4ModuleLabels")),
    _processNames(pset.get<std::vector<std::string> >("processNames")),
    _lookupTableFileNames(pset.get<std::vector<std::string> >("lookupTableFileNames")),
    _lookupTableCounterLengths(pset.get<std::vector<double> >("lookupTableCounterLengths")),
    _scintillationYield(pset.get<double>("scintillationYield")),    //850.0 photons per MeV
    _scintillatorDecayTimeFast(pset.get<double>("scintillatorDecayTimeFast")), //3.0 ns
    _scintillatorDecayTimeSlow(pset.get<double>("scintillatorDecayTimeSlow")), //10.0 ns
    _fiberDecayTime(pset.get<double>("fiberDecayTime")),     //7.4 ns
    _startTime(pset.get<double>("startTime")),               //0.0 ns
    _timeOffsets(pset.get<fhicl::ParameterSet>("timeOffsets", fhicl::ParameterSet())),
    _randFlat(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
  {
    if(_g4ModuleLabels.size()!=_processNames.size()) throw std::logic_error("ERROR: mismatch between specified selectors (g4ModuleLabels/processNames)");

    if(_lookupTableFileNames.size()!=_lookupTableCounterLengths.size()) throw std::logic_error("ERROR: mismatch between specified lookup tables (lookupTableFileNames/lookupTableCounterLengths)");
    for(unsigned int i=0; i<_lookupTableFileNames.size(); i++)
    {
      double counterLength = _lookupTableCounterLengths[i];
      _makeCrvPhotonArrivals.emplace(counterLength, boost::shared_ptr<MakeCrvPhotonArrivals>(new MakeCrvPhotonArrivals(_randFlat)));
      std::map<double, boost::shared_ptr<MakeCrvPhotonArrivals> >::iterator iterCPA=_makeCrvPhotonArrivals.find(counterLength);
      iterCPA->second->LoadLookupTable(_lookupTableFileNames[i].c_str());
      iterCPA->second->SetScintillationYield(_scintillationYield);
      iterCPA->second->SetScintillatorDecayTimeFast(_scintillatorDecayTimeFast);
      iterCPA->second->SetScintillatorDecayTimeSlow(_scintillatorDecayTimeSlow);
      iterCPA->second->SetFiberDecayTime(_fiberDecayTime);
    }
    produces<CrvPhotonArrivalsCollection>();
  }

  void CrvPhotonArrivalsGenerator::beginJob()
  {
  }

  void CrvPhotonArrivalsGenerator::endJob()
  {
  }

  void CrvPhotonArrivalsGenerator::beginRun(art::Run& rr)
  {
    GeomHandle<CosmicRayShield> CRS;
    std::vector<CRSScintillatorShield> const &shields = CRS->getCRSScintillatorShields();
    std::vector<CRSScintillatorShield>::const_iterator ishield;
    for(ishield=shields.begin(); ishield!=shields.end(); ++ishield) 
    {
      if(ishield->getCRSScintillatorBarDetail().getMaterialName()!="G4_POLYSTYRENE")
      {
        throw std::logic_error("scintillator material is not the expected G4_POLYSTYRENE which is used in the look-up tables");
      }
    }
  }

  void CrvPhotonArrivalsGenerator::produce(art::Event& event) 
  {
    _timeOffsets.updateMap(event);

    std::unique_ptr<CrvPhotonArrivalsCollection> crvPhotonArrivalsCollection(new CrvPhotonArrivalsCollection);

    GeomHandle<CosmicRayShield> CRS;
    std::string productInstanceName("CRV");

    std::vector<art::Handle<StepPointMCCollection> > CRVStepsVector;
    std::unique_ptr<art::Selector> selector;
    for(size_t j=0; j<_g4ModuleLabels.size(); j++)
    {
      if(_g4ModuleLabels[j]!="" && _g4ModuleLabels[j]!="*")
        selector = std::unique_ptr<art::Selector>(new art::Selector(art::ProductInstanceNameSelector("CRV") &&
                                                                    art::ModuleLabelSelector(_g4ModuleLabels[j]) && 
                                                                    art::ProcessNameSelector(_processNames[j])));
      else
        selector = std::unique_ptr<art::Selector>(new art::Selector(art::ProductInstanceNameSelector("CRV") &&
                                                                    art::ProcessNameSelector(_processNames[j])));
      //the ProcessNameSelector allows "*" and ""

      event.getMany(*selector, CRVStepsVector);
      for(size_t i=0; i<CRVStepsVector.size(); i++)
      {
        const art::Handle<StepPointMCCollection> &CRVSteps = CRVStepsVector[i];

        for(StepPointMCCollection::const_iterator iter=CRVSteps->begin(); iter!=CRVSteps->end(); iter++)
        {
          StepPointMC const& step(*iter);

          if(step.time()<_startTime) continue;   //Ignore this StepPoint to reduce computation time.

          const CLHEP::Hep3Vector &p1 = step.position();
          CLHEP::Hep3Vector p2 = p1 + step.momentum().unit()*step.stepLength();
          double energyDepositedTotal= step.totalEDep();
          double energyDepositedNonIonizing = step.nonIonizingEDep();

          GlobalConstantsHandle<ParticleDataTable> particleDataTable;
          int PDGcode = step.simParticle()->pdgId();
          ParticleDataTable::maybe_ref particle = particleDataTable->particle(PDGcode);
          if(!particle) 
          {
            std::cerr<<"Error in CrvPhotonArrivalsGenerator: Found a PDG code which is not in the GEANT particle table: ";
            std::cerr<<PDGcode<<std::endl;
            continue;
          }
          double mass = particle.ref().mass();  //MeV/c^2
          double charge = particle.ref().charge(); //in units of elementary charges 

          double momentum1 = step.momentum().mag(); //MeV/c
          double energy1 = sqrt(momentum1*momentum1 + mass*mass); //MeV
//FIXME: does not take the energy of daughter particles into account
          double energy2 = energy1 - energyDepositedTotal; //MeV  
          if(energy2<mass) energy2=mass;

          double gamma1 = energy1 / mass;
          double gamma2 = energy2 / mass;
          double beta1 = sqrt(1.0-1.0/(gamma1*gamma1));
          double beta2 = sqrt(1.0-1.0/(gamma2*gamma2));
          double beta = (beta1+beta2)/2.0;
          double velocity = beta*CLHEP::c_light;
          double t1 = _timeOffsets.timeWithOffsetsApplied(step);
          double t2 = t1 + step.stepLength()/velocity;

//if there is a following step point, it will give a more realistic energy and time
          StepPointMCCollection::const_iterator iterNextStep = iter+1;
          if(iterNextStep!=CRVSteps->end())
          {
            StepPointMC const& nextStep(*iterNextStep);
            if(nextStep.barIndex()==step.barIndex() && nextStep.simParticle()->id()==step.simParticle()->id())
            {
	      t2 = _timeOffsets.timeWithOffsetsApplied(nextStep);
              velocity = (p2-p1).mag()/(t2-t1);
              beta = velocity/CLHEP::c_light;
            }
          }

          const CRSScintillatorBar &CRSbar = CRS->getBar(step.barIndex());
          const CLHEP::Hep3Vector &p1Local = CRSbar.toLocal(p1);
          const CLHEP::Hep3Vector &p2Local = CRSbar.toLocal(p2);

          double counterLength = CRSbar.getHalfLength()*2.0;
          std::map<double, boost::shared_ptr<MakeCrvPhotonArrivals> >::iterator iterCPA=_makeCrvPhotonArrivals.find(counterLength);
          if(iterCPA!=_makeCrvPhotonArrivals.end())
          {
            iterCPA->second->MakePhotons(p1Local, p2Local, t1, t2,  
                                        PDGcode, beta, charge,
                                        energyDepositedTotal,
                                        energyDepositedNonIonizing);

            CrvPhotonArrivals &crvPhotons = (*crvPhotonArrivalsCollection)[step.barIndex()];
            for(int SiPM=0; SiPM<4; SiPM++)
            {
              const std::vector<double> &times=iterCPA->second->GetArrivalTimes(SiPM);
              crvPhotons.GetPhotonArrivalTimes(SiPM).insert(crvPhotons.GetPhotonArrivalTimes(SiPM).end(),times.begin(),times.end());
            }
          }
          else std::cout<<"ERROR: Did not find the matching lookuptable file for the CRV counter with the length "<<counterLength<<"."<<std::endl;

        } //loop over StepPointMCs in the StepPointMC collection
      } //loop over all StepPointMC collections
    } //loop over all module labels / process names from the fcl file

    event.put(std::move(crvPhotonArrivalsCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvPhotonArrivalsGenerator;
DEFINE_ART_MODULE(CrvPhotonArrivalsGenerator)
