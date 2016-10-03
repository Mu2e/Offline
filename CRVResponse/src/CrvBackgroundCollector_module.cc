//
// A module to collect background StepPoints and put the deposited energies into an NTuple
// (Ignored the Cerenkov part)
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>

#include <TTree.h>


namespace mu2e 
{
  class CrvBackgroundCollector : public art::EDAnalyzer 
  {

    public:
    explicit CrvBackgroundCollector(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& e);
    void beginRun(const art::Run& r);
    void endJob();

    private:
    std::vector<std::string> _g4ModuleLabels;
    std::vector<std::string> _processNames;

    std::vector<std::string> _CRVSectorNames;
    double                   _startTime;
    int                      _counters;

    TTree                   *_tree;
    ULong64_t                _eventNumber, _volumeId;
    double                   _t1, _x1, _y1, _z1;
    double                   _t2, _x2, _y2, _z2;
    double                   _beta, _charge, _totalE, _nonIonizingE;
    int                      _PDGcode;
  };

  CrvBackgroundCollector::CrvBackgroundCollector(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _g4ModuleLabels(pset.get<std::vector<std::string> >("g4ModuleLabels")),
    _processNames(pset.get<std::vector<std::string> >("processNames")),
    _CRVSectorNames(pset.get<std::vector<std::string> >("CRVSectorNames")),
    _startTime(pset.get<double>("startTime"))
  {
    if(_g4ModuleLabels.size()!=_processNames.size()) throw std::logic_error("ERROR: mismatch between specified selectors (g4ModuleLabels/processNames)");

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("background");
    _tree = tfdir.make<TTree>("background","background");
    _tree->Branch("eventNumber", &_eventNumber);
    _tree->Branch("volumeId", &_volumeId);
    _tree->Branch("t1", &_t1);
    _tree->Branch("x1", &_x1);
    _tree->Branch("y1", &_y1);
    _tree->Branch("z1", &_z1);
    _tree->Branch("t2", &_t2);
    _tree->Branch("x2", &_x2);
    _tree->Branch("y2", &_y2);
    _tree->Branch("z2", &_z2);
    _tree->Branch("beta", &_beta);
    _tree->Branch("charge", &_charge);
    _tree->Branch("totalE", &_totalE);
    _tree->Branch("nonIonizingE", &_nonIonizingE);
    _tree->Branch("PDGcode", &_PDGcode);
  }

  void CrvBackgroundCollector::beginRun(const art::Run& r)
  {
    GeomHandle<CosmicRayShield> CRS;

    _counters=0;

    const std::vector<CRSScintillatorShield> &sectors = CRS->getCRSScintillatorShields();
    std::vector<CRSScintillatorShield>::const_iterator iterSectors;
    for(iterSectors=sectors.begin(); iterSectors!=sectors.end(); iterSectors++)
    {
      const CRSScintillatorShield &sector = *iterSectors;
      if(std::find(_CRVSectorNames.begin(), _CRVSectorNames.end(), sector.getName())==_CRVSectorNames.end()) continue;

      const std::vector<CRSScintillatorModule> &modules = sector.getCRSScintillatorModules();
      std::vector<CRSScintillatorModule>::const_iterator iterModules;
      for(iterModules=modules.begin(); iterModules!=modules.end(); iterModules++)
      {
        const CRSScintillatorModule &module = *iterModules;
        const std::vector<CRSScintillatorLayer> &layers = module.getLayers();
        std::vector<CRSScintillatorLayer>::const_iterator iterLayers;
        for(iterLayers=layers.begin(); iterLayers!=layers.end(); iterLayers++)
        {
          const CRSScintillatorLayer &layer = *iterLayers;
          int n = layer.nBars();
          _counters+=n;
        }
      }
    }
  }

  void CrvBackgroundCollector::endJob()
  {
    std::cout<<"#counters: "<<_counters<<std::endl;
  }

  void CrvBackgroundCollector::analyze(const art::Event& event) 
  {
    GeomHandle<CosmicRayShield> CRS;

    mu2e::ConditionsHandle<mu2e::AcceleratorParams> accPar("ignored");
    double microBunchPeriod = accPar->deBuncherPeriod;

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

          _t1 = step.time();
          _t1 = fmod(_t1,microBunchPeriod);
          if(_t1<_startTime) continue;

          const CRSScintillatorBar &CRSbar = CRS->getBar(step.barIndex());

          int sectorNumber = CRSbar.id().getShieldNumber();
          const std::string sectorName = CRS->getCRSScintillatorShield(sectorNumber).getName();

          if(std::find(_CRVSectorNames.begin(), _CRVSectorNames.end(), sectorName)==_CRVSectorNames.end()) continue;

          const CLHEP::Hep3Vector &p1 = step.position();
          CLHEP::Hep3Vector p2 = p1 + step.momentum().unit()*step.stepLength();
          _totalE = step.totalEDep();
          _nonIonizingE = step.nonIonizingEDep();

          GlobalConstantsHandle<ParticleDataTable> particleDataTable;
          _PDGcode = step.simParticle()->pdgId();
          ParticleDataTable::maybe_ref particle = particleDataTable->particle(_PDGcode);
          if(!particle) 
          {
            std::cerr<<"Error in CrvPhotonArrivalsGenerator: Found a PDG code which is not in the GEANT particle table: ";
            std::cerr<<_PDGcode<<std::endl;
            continue;
          }
          double mass = particle.ref().mass();  //MeV/c^2
          _charge = particle.ref().charge(); //in units of elementary charges 

          double momentum1 = step.momentum().mag(); //MeV/c
          double energy1 = sqrt(momentum1*momentum1 + mass*mass); //MeV
//FIXME: does not take the energy of daughter particles into account
          double energy2 = energy1 - _totalE; //MeV  
          if(energy2<mass) energy2=mass;

          double gamma1 = energy1 / mass;
          double gamma2 = energy2 / mass;
          double beta1 = sqrt(1.0-1.0/(gamma1*gamma1));
          double beta2 = sqrt(1.0-1.0/(gamma2*gamma2));
          _beta = (beta1+beta2)/2.0;
          double velocity = _beta*CLHEP::c_light;
          _t2 = _t1 + step.stepLength()/velocity;

          const CLHEP::Hep3Vector &p1Local = CRSbar.toLocal(p1);
          const CLHEP::Hep3Vector &p2Local = CRSbar.toLocal(p2);

          _x1 = p1Local.x();
          _y1 = p1Local.y();
          _z1 = p1Local.z();
          _x2 = p2Local.x();
          _y2 = p2Local.y();
          _z2 = p2Local.z();

          _volumeId = step.volumeId();
          _eventNumber = event.event();

          _tree->Fill();
        }
      }
    }
  }

} // end namespace mu2e

using mu2e::CrvBackgroundCollector;
DEFINE_ART_MODULE(CrvBackgroundCollector)
