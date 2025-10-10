//
// Decide which physics list to use.
//
//
// Original author Rob Kutschke
//
// Notes:
// 1) Extract the name of the requested physics list from the
//    config file, instantiate the requested physics list and
//    return a bare pointer to it.
//
// 2) The caller receives the pointer and immediately passes it
//    to G4, which takes ownership of the physics list object.
//    The G4 interface requires a bare pointer.
//
// 3) There are same special names:
//     Minimal - the original Mu2e minimal physics list
//     *_MU2E* Mu2e custom lists
//
// 4) All other names are presumed to be valid names for physics lists that
//    can be created by the PhysListFactory.  At this writing ( April 2010),
//    the PhysListFactory has a default if it does not recognize name as
//    one of its known lists.  The default is QGSP_BERT 3.3.
//
// C++ includes
#include <string>

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Offline/Mu2eG4/inc/physicsListDecider.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4DecayMuonsWithSpinPhysicsConstructor.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4MinimalModularPhysicsList.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4MinDEDXModularPhysicsList.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4StepLimiterPhysicsConstructor.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4CustomizationPhysicsConstructor.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4BiasedRDPhysics.hh"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// G4 includes
#include "Geant4/G4PhysListFactory.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4RadioactiveDecayPhysics.hh"
#include "Geant4/G4ThermalNeutrons.hh"
#include "Geant4/G4ErrorPhysicsList.hh"
#include "Geant4/G4EmStandardPhysics_option4.hh"

#if G4VERSION>4103
#include "Geant4/G4EmParameters.hh"
#endif

#if G4VERSION>4106
#include "Geant4/G4HadronicParameters.hh"
#endif

using namespace std;

namespace mu2e{

  G4VUserPhysicsList* physicsListDecider(const Mu2eG4Config::Physics& phys
                                         , const Mu2eG4Config::Debug& debug
                                         , const Mu2eG4ResourceLimits& lim) {

    G4VModularPhysicsList* tmpPL(nullptr);
    const string name = phys.physicsListName();

    debug.diagLevel()>-1 && G4cout << __func__ << " invoked with " << name << G4endl;

    // special cases
    if ( name  == "Minimal" ) {
      return new Mu2eG4MinimalModularPhysicsList();
    }

    else if ( name  == "MinDEDX" ) {
      return new Mu2eG4MinDEDXModularPhysicsList(); // limited EM Processes
    }

    else if ( name  == "ErrorPhysicsList" ) {
      // rather special case of G4VUserPhysicsList for Track Error
      // Propagation, with special Energy Loss implementation
      // (see User's Guide: For Application Developers)
      return new G4ErrorPhysicsList();
    }

    // General case
    else {

    // the Hadronic params have to be set before defining the list

#if G4VERSION>4112
    { bool BertiniAs11_2 = true; // restores 11.2 behavior in Geant4 11.3.p02
      phys.setBertiniAs11_2(BertiniAs11_2);
      if (BertiniAs11_2) {
        if (debug.diagLevel()>0) {
          G4cout << __func__
                 << " Setting Bertini model behavior to one of 11.2 "
                 << G4endl;
        }
        G4HadronicParameters::Instance()->SetBertiniAs11_2(BertiniAs11_2);
      }
    }
#endif

      G4PhysListFactory physListFactory;
      physListFactory.SetVerbose(debug.diagLevel());
      tmpPL = physListFactory.GetReferencePhysList(name);

    }

    if ( tmpPL==nullptr ) {
      throw cet::exception("G4CONTROL")
        << "Unable to load physics list named: "
        << name
        << "\n";
    }

    // The modular physics list takes ownership of the Mu2eG4StepLimiterPhysicsConstructor.
    tmpPL->RegisterPhysics( new Mu2eG4StepLimiterPhysicsConstructor() );

    // Mu2e Customizations
    tmpPL->RegisterPhysics( new Mu2eG4CustomizationPhysicsConstructor(&phys, &debug, &lim));

    if ( phys.turnOffRadioactiveDecay() && phys.turnOnRadioactiveDecay() ) {
      mf::LogError("Config") << "Inconsistent config";
      G4cout << "Error: turnOnRadioactiveDecay & turnOffRadioactiveDecay on" << G4endl;
      throw cet::exception("BADINPUT")<<" decide on turnOn/OffRadioactiveDecay\n";
    }

    if ( phys.turnOnRadioactiveDecay() && phys.radiationVRmode() ) {
      mf::LogError("Config") << "Inconsistent config";
      G4cout << "Error: turnOnRadioactiveDecay & radiationVRmode on" << G4endl;
      throw cet::exception("BADINPUT")<<" decide on one option\n";
    }

    if ( phys.turnOffRadioactiveDecay() && phys.radiationVRmode() ) {
      mf::LogError("Config") << "Inconsistent config";
      G4cout << "Error: turnOffRadioactiveDecay & radiationVRmode on" << G4endl;
      throw cet::exception("BADINPUT")<<" decide on one option\n";
    }

    if (phys.turnOffRadioactiveDecay()) {
      tmpPL->RemovePhysics("G4RadioactiveDecay");
    }

    if (phys.turnOnRadioactiveDecay()) {
      // turn it off first to avoid warning if this process is already included in the
      // physics list, does nothing if the process is absent
      tmpPL->RemovePhysics("G4RadioactiveDecay");
      tmpPL->RegisterPhysics(new G4RadioactiveDecayPhysics(debug.diagLevel()));
    }

    if (phys.radiationVRmode()){
      tmpPL->RemovePhysics("G4RadioactiveDecay");
      tmpPL->RegisterPhysics(new Mu2eG4BiasedRDPhysics(&phys, debug.diagLevel()));
    }

    if (phys.turnOnThermalNeutronPhysics()) {
      tmpPL->RegisterPhysics(new G4ThermalNeutrons(debug.diagLevel()));
    }

#if G4VERSION>4104

    // Changing MSC model transition energy if requested
    { double mscModelTransitionEnergy(std::numeric_limits<double>::max());  // initializing
      // to something distinct before fetching the requested value if present and only using it then
      if (phys.mscModelTransitionEnergy(mscModelTransitionEnergy)) {
        if (debug.diagLevel()>0) {
          G4cout << __func__
                 << " Changing MscEnergyLimit to "
                 << mscModelTransitionEnergy << " MeV" << G4endl;
        }
        G4EmParameters* emParams = G4EmParameters::Instance();
        emParams->SetMscEnergyLimit(mscModelTransitionEnergy*CLHEP::MeV);
      }
    }

    if ( phys.useEmOption4InTracker() && (name.find("_EMZ") == std::string::npos) ) {
      // assign EmStandard_opt4 to the tracker
      if (debug.diagLevel()>0) {
        G4cout << __func__
               << " Assigning EmStandardPhysics_option4 to the TrackerMother" << G4endl;
      }
      G4EmParameters* emParams = G4EmParameters::Instance();
      emParams->AddPhysics("TrackerMother", "G4EmStandard_opt4");
    }

#endif

    // Muon Spin and Radiative decays plus pion muons with spin
    if ( phys.decayMuonsWithSpin() ) {

      // requires spin tracking
      if ( phys.stepper() != "G4ClassicalRK4WSpin" &&
           phys.stepper() != "G4DormandPrince745WSpin" &&
           phys.stepper() != "G4TDormandPrince45WSpin" ) {
        mf::LogError("Config") << "Inconsistent config";
        G4cout << "Error: Mu2eG4DecayMuonsWithSpinPhysicsConstructor requires enabling spin tracking" << G4endl;
        throw cet::exception("BADINPUT")<<" Mu2eG4DecayMuonsWithSpinPhysicsConstructor requires enabling spin tracking\n";
      }

      tmpPL->RegisterPhysics( new Mu2eG4DecayMuonsWithSpinPhysicsConstructor(debug.diagLevel()));
    }

    G4double productionCut = phys.minRangeCut();
    G4double protonProductionCut = phys.protonProductionCut();
    mf::LogInfo("GEOM_MINRANGECUT")
      << "Setting production cut to " << productionCut
      << ", protonProductionCut to " << protonProductionCut << " mm";
    if (debug.diagLevel() > 0) {
      G4cout << __func__ << " Setting gamma, e- and e+ production cut"
             << " to " << productionCut << " mm and for proton to "
             << protonProductionCut << " mm" << G4endl;
    }
    //setCutCmd equivalent:
    tmpPL->SetDefaultCutValue(productionCut);
    tmpPL->SetCutValue(protonProductionCut, "proton");
    // regional cuts (if any) are set during the geometry construction

    if (debug.diagLevel() > 0) tmpPL->DumpCutValuesTable();

    return dynamic_cast<G4VUserPhysicsList*>(tmpPL);

  }

} // end namespace mu2e
