// Original author  K.L.Genser
// A constructor implementing Mu2e specific customizations

// C++ includes
#include <set>
#include <string>
#include <vector>

// Mu2e includes
#include "Offline/Mu2eG4/inc/Mu2eG4CustomizationPhysicsConstructor.hh"
#include "Offline/Mu2eG4/inc/customizeChargedPionDecay.hh"
#include "Offline/Mu2eG4/inc/toggleProcesses.hh"

// Framework includes
#include "cetlib_except/exception.h"

// G4 includes
#include "Geant4/G4BuilderType.hh"
#include "Geant4/G4Threading.hh"
#include "Geant4/G4EmParameters.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4ProductionCuts.hh"
#include "Geant4/G4Region.hh"
#include "Geant4/G4RegionStore.hh"
#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4VEmProcess.hh"
#include "Geant4/G4VEnergyLossProcess.hh"

// factory
#include "Geant4/G4PhysicsConstructorFactory.hh"
//

namespace {

  // Activate Geant4 EM secondary biasing for each requested process in every listed region.
  // G4EmParameters (and so /process/em/setSecBiasing) activates only the first region
  // given for a process, so the process objects are configured directly here.
  void activateSecondaryBiasing(const mu2e::Mu2eG4Config::Physics& phys,
                                const mu2e::Mu2eG4Config::Debug& debug) {

    std::vector<mu2e::Mu2eG4Config::SecondaryBiasing> requests;
    if (!phys.secondaryBiasing(requests)) { return; }

    // Geant4 maps material-cuts couples to a biased region by the production-cuts pointer,
    // and ActivateSecondaryBiasing only prints a warning for an unknown region.  Check both
    // here: every region must exist and own production cuts that no other region uses,
    // otherwise biasing would silently be skipped or leak into other regions.
    G4RegionStore* regionStore = G4RegionStore::GetInstance();
    std::set<std::string> biasedRegions;
    for (const auto& req : requests) {
      for (const auto& rname : req.regions()) { biasedRegions.insert(rname); }
    }
    for (const auto& rname : biasedRegions) {
      const G4Region* region = regionStore->GetRegion(rname, false);
      if (region == nullptr) {
        throw cet::exception("BADINPUT")
          << "physics.secondaryBiasing: unknown G4Region " << rname
          << " (regions are made by physics.minRangeRegionCuts)\n";
      }
      const G4ProductionCuts* cuts = region->GetProductionCuts();
      if (cuts == nullptr) {
        throw cet::exception("BADINPUT")
          << "physics.secondaryBiasing: G4Region " << rname
          << " has no production cuts of its own, so Geant4 would bias every volume that"
          << " shares the default cuts; give it an entry in physics.minRangeRegionCuts\n";
      }
      for (const G4Region* other : *regionStore) {
        if (other != region && other->GetProductionCuts() == cuts
            && biasedRegions.count(other->GetName()) == 0) {
          throw cet::exception("BADINPUT")
            << "physics.secondaryBiasing: G4Region " << rname
            << " shares its production cuts with unbiased G4Region " << other->GetName() << "\n";
        }
      }
    }

    const bool print = G4Threading::IsMasterThread() || debug.diagLevel() > 0;
    for (const auto& req : requests) {
      const std::string pname = req.process();
      const double factor = req.factor();
      const double energyLimit = req.energyLimit()*CLHEP::MeV;
      if (!(factor > 0.)) {
        throw cet::exception("BADINPUT")
          << "physics.secondaryBiasing: factor for " << pname << " must be > 0, got " << factor << "\n";
      }

      // a process object may be shared by several particles
      std::set<const G4VProcess*> done;
      G4ParticleTable::G4PTblDicIterator* iter = G4ParticleTable::GetParticleTable()->GetIterator();
      iter->reset();
      while ( (*iter)() ) {
        G4ParticleDefinition* particle = iter->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        if (pmanager == nullptr) { continue; }
        G4ProcessVector* plist = pmanager->GetProcessList();
        for (G4int j = 0; j < pmanager->GetProcessListLength(); ++j) {
          G4VProcess* proc = (*plist)[j];
          if (proc->GetProcessName() != pname || done.count(proc) > 0) { continue; }
          done.insert(proc);
          auto* eloss = dynamic_cast<G4VEnergyLossProcess*>(proc);
          auto* emproc = dynamic_cast<G4VEmProcess*>(proc);
          if (eloss == nullptr && emproc == nullptr) {
            throw cet::exception("BADINPUT")
              << "physics.secondaryBiasing: process " << pname << " of "
              << particle->GetParticleName()
              << " is neither a G4VEnergyLossProcess nor a G4VEmProcess\n";
          }
          for (const auto& rname : req.regions()) {
            if (eloss != nullptr) { eloss->ActivateSecondaryBiasing(rname, factor, energyLimit); }
            else { emproc->ActivateSecondaryBiasing(rname, factor, energyLimit); }
          }
          if (print) {
            G4cout << "Mu2eG4 secondaryBiasing: " << pname << " of " << particle->GetParticleName()
                   << " factor " << factor << " energyLimit " << energyLimit/CLHEP::MeV
                   << " MeV in G4Regions:";
            for (const auto& rname : req.regions()) { G4cout << " " << rname; }
            G4cout << G4endl;
          }
        }
      }
      if (done.empty()) {
        throw cet::exception("BADINPUT")
          << "physics.secondaryBiasing: no particle has a process named " << pname << "\n";
      }
    }
  }

} // namespace

namespace mu2e {

  G4_DECLARE_PHYSCONSTR_FACTORY(Mu2eG4CustomizationPhysicsConstructor);

  Mu2eG4CustomizationPhysicsConstructor::
  Mu2eG4CustomizationPhysicsConstructor(const Mu2eG4Config::Physics* phys
                                        , const Mu2eG4Config::Debug* debug
                                        , const Mu2eG4ResourceLimits* lim
                                        )
    : G4VPhysicsConstructor("Mu2eCustomizations",bUnknown)
    , phys_(phys)
    , debug_(debug)
    , mu2elimits_(lim)
      // bUnknown disables duplication check for physics builder type
  {}

  // default constructor required by G4PhysicsConstructorFactory; fixme: check, also pset_?
  Mu2eG4CustomizationPhysicsConstructor::
  Mu2eG4CustomizationPhysicsConstructor()
    : G4VPhysicsConstructor("Mu2eCustomizations",bUnknown)
    , phys_(nullptr)
    , debug_(nullptr)
    , mu2elimits_(nullptr)
  {}

  void Mu2eG4CustomizationPhysicsConstructor::ConstructParticle()
  {
    // Empty on purpose, for now
  }

  void Mu2eG4CustomizationPhysicsConstructor::ConstructProcess()
  {

    if (debug_->diagLevel()>0) {
      G4cout << "Mu2eG4CustomizationPhysicsConstructor::"
             << __func__ << " Called" << G4endl;
    }

    if(G4Threading::IsMasterThread()) {

      // G4 does not include pi+ -> e+ nu + cc. Fix that in one of several ways.
      // decay tables are common
      customizeChargedPionDecay(*phys_, *debug_);

    }

    // Switch off the decay of some particles
    switchDecayOff(*phys_, *debug_);

    // add conditional and special mu2e processes
    addUserProcesses(*phys_, *debug_, *mu2elimits_);

    // swap Bertini Cascade with Precompound model in G4MuonMinusCapture if requested
    switchCaptureDModel(*phys_, *debug_);

    // EM secondary splitting in several G4Regions if requested
    activateSecondaryBiasing(*phys_, *debug_);

    // use more accurate boundary crossing algorithm for muons and hadrons
    // place this late in the setup sequence to avoid a subsequent change
    if (phys_->setMuHadLateralDisplacement()) {
      G4EmParameters* params = G4EmParameters::Instance();
      params->SetMuHadLateralDisplacement(true);
    }

    { bool disableEnergyLossFluctuations = false; // they are enabled by deafult
      if (phys_->disableEnergyLossFluctuations(disableEnergyLossFluctuations)) {
        if (disableEnergyLossFluctuations) {
          if (debug_->diagLevel()>0) {
            G4cout << __func__
                   << " Turning off energy loss fluctuations "
                   << G4endl;
          }
          G4EmParameters* params = G4EmParameters::Instance();
          params->SetLossFluctuations(!disableEnergyLossFluctuations);
        }
      }
    }

#if G4VERSION>4110
    { unsigned int energyLossFluctuationModel = 2; // 2 is geant4 11.1 default
      if (phys_->setEnergyLossFluctuationModel(energyLossFluctuationModel)) {
        if (debug_->diagLevel()>0) {
          G4cout << __func__
                 << " Setting energy loss fluctuations model to: "
                 << energyLossFluctuationModel
                 << G4endl;
        }
        G4EmParameters* params = G4EmParameters::Instance();
        // a potential maintance point, 3 options as of Geant4 11.1
        G4EmFluctuationType fluctType = fDummyFluctuation; // actually a model
        if (energyLossFluctuationModel == 0) { fluctType = fDummyFluctuation; }
        else if (energyLossFluctuationModel == 1) { fluctType = fUniversalFluctuation; }
        else if (energyLossFluctuationModel == 2) { fluctType = fUrbanFluctuation; }
        else {
          throw cet::exception("BADINPUT")
            << "Mu2eG4CustomizationPhysicsConstructor::ConstructProcess() Error: unknown energyLossFluctuationModel: "
            << energyLossFluctuationModel
            << "\n";
        }
        params->SetFluctuationType(fluctType);
      }
    }
#endif

  }

}
