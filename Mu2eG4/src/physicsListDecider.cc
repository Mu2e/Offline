//
// Decide which physics list to use.
//
// $Id: physicsListDecider.cc,v 1.18 2014/03/04 04:26:36 genser Exp $
// $Author: genser $
// $Date: 2014/03/04 04:26:36 $
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
// 3) There are two special names:
//     Minimal - the original Mu2e minimal physics list
//     N02     - the physics list copied from the G4 novice example N02.
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
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/DecayMuonsWithSpin.hh"
#include "Mu2eG4/inc/MinimalPhysicsList.hh"
#include "Mu2eG4/inc/StepLimiterPhysConstructor.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "fhiclcpp/ParameterSet.h"

//tmp arrangement
#include "Mu2eG4/inc/QGSP_BERT_HP_MU2E00.hh"
#include "Mu2eG4/inc/QGSP_BERT_MU2E00.hh"
#include "Mu2eG4/inc/Shielding_MU2E00.hh"
#include "Mu2eG4/inc/Shielding_MU2E01.hh"
#include "Mu2eG4/inc/Shielding_MU2E02.hh"
#include "Mu2eG4/inc/FTFP_BERT_PBAR_MU2E02.hh"

// G4 includes
#include "G4PhysListFactory.hh"
#include "G4VUserPhysicsList.hh"
#if G4VERSION<4099
#include "QGSP.hh"
#endif

using namespace std;

namespace mu2e{
  namespace {
    std::string getPhysicsListName(const SimpleConfig& config) {
      return config.getString("g4.physicsListName");
    }

    std::string getPhysicsListName(const fhicl::ParameterSet& pset) {
      return pset.get<std::string>("physics.physicsListName");
    }

    bool turnOffRadioactiveDecay(const SimpleConfig& config) {
      return config.getBool("g4.turnOffRadioactiveDecay",false);
    }

    bool turnOffRadioactiveDecay(const fhicl::ParameterSet& pset) {
      return pset.get<bool>("physics.turnOffRadioactiveDecay",false);
    }

    int getDiagLevel(const SimpleConfig& config) {
      return config.getInt("g4.diagLevel");
    }

    int getDiagLevel(const fhicl::ParameterSet& pset) {
      return pset.get<int>("debug.diagLevel");
    }

    std::string getStepperName(const SimpleConfig& config) {
       return config.getString("g4.stepper");
    }

    std::string getStepperName(const fhicl::ParameterSet& pset) {
      return pset.get<std::string>("physics.stepper");
    }

  }


  template<class Config>
  G4VUserPhysicsList* physicsListDecider (const Config& config){

    G4VUserPhysicsList* physicsList(0);

    const string name = getPhysicsListName(config);

    // Two special cases
    if ( name  == "Minimal" ) {
      physicsList = dynamic_cast<G4VUserPhysicsList*>(new MinimalPhysicsList );
    }

#if G4VERSION<4099
    else if ( name == "QGSP" ){
      G4VModularPhysicsList* tmp = new QGSP();
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );
      physicsList = tmp;
    }
#endif

    else if ( name == "QGSP_BERT_MU2E00" ){
      G4VModularPhysicsList* tmp = new QGSP_BERT_MU2E00();
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );
      mf::LogWarning("PHYS") << "This Mu2e Physics List has not been certified";
      G4cout << "Warning: This Mu2e Physics List has not been certified" << G4endl;
      physicsList = tmp;
    }

    else if ( name == "QGSP_BERT_HP_MU2E00" ){
      G4VModularPhysicsList* tmp = new TQGSP_BERT_HP_MU2E00<G4VModularPhysicsList,Config>(config);
      mf::LogWarning("PHYS") << "This Mu2e Physics List has not been certified";
      G4cout << "Warning: This Mu2e Physics List has not been certified" << G4endl;
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );
      physicsList = tmp;
    }

    else if ( name == "Shielding_MU2E00" ){
      G4VModularPhysicsList* tmp = new Shielding_MU2E00();
#if G4VERSION>4099
      mf::LogWarning("PHYS") << "This Mu2e Physics List has not been certified for use with Geant4 v10+.";
      G4cout << "Warning: This Mu2e Physics List has not been certified for use with Geant4 v10+." << G4endl;
#endif
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );
      physicsList = tmp;
    }

    else if ( name == "Shielding_MU2E01" ){
      G4VModularPhysicsList* tmp = new TShielding_MU2E01<G4VModularPhysicsList,Config>(config);
#if G4VERSION>4099
      mf::LogWarning("PHYS") << "This Mu2e Physics List has not been certified for use with Geant4 v10+.";
      cout << "Warning: This Mu2e Physics List has not been certified for use with Geant4 v10+." << endl;
#endif
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );
      physicsList = tmp;
    }

    else if ( name == "Shielding_MU2E02" ){
      G4VModularPhysicsList* tmp = new TShielding_MU2E02<G4VModularPhysicsList,Config>(config);
#if G4VERSION>4099
      mf::LogWarning("PHYS") << "This Mu2e Physics List has not been certified for use with Geant4 v10+.";
      cout << "Warning: This Mu2e Physics List has not been certified for use with Geant4 v10+." << endl;
#endif
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );
      physicsList = tmp;
    }

    else if ( name == "FTFP_BERT_PBAR_MU2E02" ){
      G4VModularPhysicsList* tmp = new TFTFP_BERT_PBAR_MU2E02<G4VModularPhysicsList>;
#if G4VERSION>4099
      mf::LogWarning("PHYS") << "This Mu2e Physics List has not been certified for use with Geant4 v10+.";
      cout << "Warning: This Mu2e Physics List has not been certified for use with Geant4 v10+." << endl;
#endif
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );
      physicsList = tmp;
    }

    // General case
    else {
      G4PhysListFactory physListFactory;
      G4VModularPhysicsList* tmp = physListFactory.GetReferencePhysList(name);

      // The modular physics list takes ownership of the StepLimiterPhysConstructor.
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );

      physicsList = tmp;
    }

    if ( !physicsList ){
      throw cet::exception("G4CONTROL")
        << "Unable to load physics list named: "
        << name
        << "\n";
    }

    if (turnOffRadioactiveDecay(config)) {
      (dynamic_cast<G4VModularPhysicsList*>(physicsList))->RemovePhysics("G4RadioactiveDecay");
    }

    // Muon Spin and Radiative decays plus pion muons with spin
    if ( getDecayMuonsWithSpin(config) ) {

      // requires spin tracking: G4ClassicalRK4WSpin
      if ( getStepperName(config) !=  "G4ClassicalRK4WSpin") {
        mf::LogError("Config") << "Inconsistent config";
        G4cout << "Error: DecayMuonsWithSpin requires G4ClassicalRK4WSpin stepper" << G4endl;
        throw cet::exception("BADINPUT")<<" DecayMuonsWithSpin requires G4ClassicalRK4WSpin stepper\n";
      }

      (dynamic_cast<G4VModularPhysicsList*>(physicsList))->
        RegisterPhysics( new DecayMuonsWithSpin(getDiagLevel(config)));
    }

    return physicsList;

  }

  template G4VUserPhysicsList* physicsListDecider(const SimpleConfig& config);
  template G4VUserPhysicsList* physicsListDecider(const fhicl::ParameterSet& pset);

} // end namespace mu2e
