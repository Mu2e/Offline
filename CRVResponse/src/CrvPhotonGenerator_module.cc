//
// A module to create CRV photons arriving at the SiPMs (using StepPointMCs)
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
//
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvPhotons.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CrvPhotonsCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
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

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>


namespace mu2e
{
  class CrvPhotonGenerator : public art::EDProducer
  {

    public:
    explicit CrvPhotonGenerator(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginRun(art::Run& r);

    private:
    std::vector<std::string> _g4ModuleLabels;
    std::vector<std::string> _processNames;

    ConfigFileLookupPolicy                                     _resolveFullPath;
    std::vector<std::string>                                   _lookupTableFileNames;
    std::vector<int>                                           _lookupTableReflectors;
    std::vector<std::string>                                   _lookupTableCRVSectors;
    std::vector<boost::shared_ptr<mu2eCrv::MakeCrvPhotons> >   _makeCrvPhotons;

    double      _scintillationYield;
    double      _scintillationYieldVariation;
    double      _scintillationYieldVariationCutoff;

    double      _startTime;             //StepPoint times before this time will be ignored to reduce computation times
                                        //(in particular by ignoring hits during the beam flash).
                                        //This time should be at least 100ns before the end of the SiPM's blind time
                                        //to account for the travel time of the photons inside the CRV bar.
                                        //Default is 0.

    std::string _visibleEnergyAdjustmentFileName;

    SimParticleTimeOffset _timeOffsets;

    CLHEP::HepRandomEngine& _engine;
    CLHEP::RandFlat       _randFlat;
    CLHEP::RandGaussQ     _randGaussQ;
    CLHEP::RandPoissonQ   _randPoissonQ;

    std::map<CRSScintillatorBarIndex,double>  _scintillationYieldAdjustments;
  };

  CrvPhotonGenerator::CrvPhotonGenerator(fhicl::ParameterSet const& pset) :
    EDProducer{pset},
    _g4ModuleLabels(pset.get<std::vector<std::string> >("g4ModuleLabels")),
    _processNames(pset.get<std::vector<std::string> >("processNames")),
    _lookupTableFileNames(pset.get<std::vector<std::string> >("lookupTableFileNames")),
    _lookupTableReflectors(pset.get<std::vector<int> >("reflectors")),
    _lookupTableCRVSectors(pset.get<std::vector<std::string> >("CRVSectors")),
    _scintillationYield(pset.get<double>("scintillationYield")),    //5000.0 photons per MeV
    _scintillationYieldVariation(pset.get<double>("scintillationYieldVariation")),    //20.0%
    _scintillationYieldVariationCutoff(pset.get<double>("scintillationYieldVariationCutoff")),    //20.0%
    _startTime(pset.get<double>("startTime")),               //0.0 ns
    _visibleEnergyAdjustmentFileName(pset.get<std::string>("visibleEnergyAdjustmentFileName")),
    _timeOffsets(pset.get<fhicl::ParameterSet>("timeOffsets", fhicl::ParameterSet())),
    _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())},
    _randFlat(_engine),
    _randGaussQ(_engine),
    _randPoissonQ(_engine)
  {
    if(_g4ModuleLabels.size()!=_processNames.size()) throw std::logic_error("ERROR: mismatch between specified selectors (g4ModuleLabels/processNames)");

    if(_lookupTableFileNames.size()!=_lookupTableCRVSectors.size()) throw std::logic_error("ERROR: mismatch between specified lookup tables (lookupTableFileNames/CRVSectors)");
    if(_lookupTableReflectors.size()!=_lookupTableCRVSectors.size()) throw std::logic_error("ERROR: mismatch between specified lookup tables (reflectors/CRVSectors)");

    ConfigFileLookupPolicy configFile;
    _visibleEnergyAdjustmentFileName = configFile(_visibleEnergyAdjustmentFileName);

    for(size_t i=0; i<_lookupTableFileNames.size(); i++)
    {
      bool tableLoaded=false;
      for(size_t j=0; j<i; j++)
      {
        if(_lookupTableFileNames[i]==_lookupTableFileNames[j])
        {
           tableLoaded=true;
           _makeCrvPhotons.emplace_back(_makeCrvPhotons[j]);
           std::cout<<"CRV sector "<<i<<" ("<<_lookupTableCRVSectors[i]<<") uses "<<_makeCrvPhotons.back()->GetFileName()<<std::endl;
           break;
        }
      }
      if(tableLoaded) continue;

      _makeCrvPhotons.emplace_back(boost::shared_ptr<mu2eCrv::MakeCrvPhotons>(new mu2eCrv::MakeCrvPhotons(_randFlat, _randGaussQ, _randPoissonQ)));
      boost::shared_ptr<mu2eCrv::MakeCrvPhotons> &photonMaker=_makeCrvPhotons.back();
      photonMaker->LoadLookupTable(_resolveFullPath(_lookupTableFileNames[i]));
      photonMaker->SetScintillationYield(_scintillationYield);
      photonMaker->LoadVisibleEnergyAdjustmentTable(_visibleEnergyAdjustmentFileName);
      std::cout<<"CRV sector "<<i<<" ("<<_lookupTableCRVSectors[i]<<") uses "<<_makeCrvPhotons.back()->GetFileName()<<std::endl;
    }

    produces<CrvPhotonsCollection>();
  }

  void CrvPhotonGenerator::beginRun(art::Run& rr)
  {
    GeomHandle<CosmicRayShield> CRS;
    std::vector<CRSScintillatorShield> const &shields = CRS->getCRSScintillatorShields();
    if(shields.size()!=_lookupTableCRVSectors.size()) throw std::logic_error("ERROR: mismatch between the geometry and the specified lookup table CRVSectors");

    for(size_t i=0; i<shields.size(); i++)
    {
      if(shields[i].getCRSScintillatorBarDetail().getMaterialName()!="G4_POLYSTYRENE")
        throw std::logic_error("ERROR: scintillator material is not the expected G4_POLYSTYRENE which is used in the look-up tables");
      if(shields[i].getName().substr(4)!=_lookupTableCRVSectors[i]) throw std::logic_error("ERROR: mismatch between the geometry and the specified lookup table CRVSectors");
                            //substr(4) removes the "CRV_" part of the sector name
    }
  }

  void CrvPhotonGenerator::produce(art::Event& event)
  {
    _timeOffsets.updateMap(event);

    _scintillationYieldAdjustments.clear();

    std::unique_ptr<CrvPhotonsCollection> crvPhotonsCollection(new CrvPhotonsCollection);

    GeomHandle<CosmicRayShield> CRS;

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
        for(size_t istep=0; istep<CRVSteps->size(); istep++)
        {
          StepPointMC const& step(CRVSteps->at(istep));

          double t1 = _timeOffsets.timeWithOffsetsApplied(step);
          if(t1<_startTime) continue;   //Ignore this StepPoint to reduce computation time.

          const CLHEP::Hep3Vector &p1 = step.position();
          CLHEP::Hep3Vector p2 = p1 + step.momentum().unit()*step.stepLength();    //this stepLength does not necessarily take us to the right p2 due to scattering
          double energyDepositedTotal= step.totalEDep();
          double energyDepositedNonIonizing = step.nonIonizingEDep();

          GlobalConstantsHandle<ParticleDataTable> particleDataTable;
          int PDGcode = step.simParticle()->pdgId();
          ParticleDataTable::maybe_ref particle = particleDataTable->particle(PDGcode);
          if(!particle)
          {
            std::cerr<<"Error in CrvPhotonGenerator: Found a PDG code which is not in the GEANT particle table: ";
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
          double t2 = t1 + step.stepLength()/velocity;

          const CRSScintillatorBar &CRSbar = CRS->getBar(step.barIndex());
          const CLHEP::Hep3Vector &p1Local = CRSbar.toLocal(p1);
          const CLHEP::Hep3Vector &p2Local = CRSbar.toLocal(p2);

          if(_scintillationYieldAdjustments.find(step.barIndex())==_scintillationYieldAdjustments.end())
          {
            double adjustment=0;
            do
            {
              adjustment=_randGaussQ.fire(0, _scintillationYield*_scintillationYieldVariation);
            } while(adjustment<-_scintillationYield*_scintillationYieldVariationCutoff);
            _scintillationYieldAdjustments[step.barIndex()] = adjustment;
          }
          double scintillationYieldAdjustment = _scintillationYieldAdjustments[step.barIndex()];

          const CRSScintillatorBarId &barId = CRSbar.id();
          int CRVSectorNumber=barId.getShieldNumber();
          boost::shared_ptr<mu2eCrv::MakeCrvPhotons> &photonMaker=_makeCrvPhotons.at(CRVSectorNumber);
          photonMaker->MakePhotons(p1Local, p2Local, t1, t2,
                                        PDGcode, beta, charge,
                                        energyDepositedTotal,
                                        energyDepositedNonIonizing,
                                        step.stepLength(),
                                        scintillationYieldAdjustment,
                                        _lookupTableReflectors[CRVSectorNumber]);

          CrvPhotons &crvPhotons = (*crvPhotonsCollection)[step.barIndex()];
          for(int SiPM=0; SiPM<4; SiPM++)
          {
            const std::vector<double> &times=photonMaker->GetArrivalTimes(SiPM);
            art::Ptr<StepPointMC> spmcptr(CRVSteps,istep);
            for(size_t itime=0; itime<times.size(); itime++)
            {
              CrvPhotons::SinglePhoton photon;
              photon._time = times[itime];
              photon._step = spmcptr;
              crvPhotons.GetPhotons(SiPM).push_back(photon);
            }
          }

        } //loop over StepPointMCs in the StepPointMC collection
      } //loop over all StepPointMC collections
    } //loop over all module labels / process names from the fcl file

/* photnns into the event */

    event.put(std::move(crvPhotonsCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvPhotonGenerator;
DEFINE_ART_MODULE(CrvPhotonGenerator)
