//
// Set the G4 range/production cuts including the ones for regions if any
//
// Original author  K.L.Genser

// C++ includes
#include <iostream>
#include <vector>

// Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// G4 includes
#include "G4VModularPhysicsList.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

// Mu2e includes
#include "Mu2eG4/inc/setMu2eG4ProductionCuts.hh"

namespace mu2e{

  // to be called while constructing a physics list
  void setMu2eG4ProductionCuts(const fhicl::ParameterSet& pset, G4VModularPhysicsList* mPL){

    mPL->SetCutsWithDefault();
    // the above calls  parent class SetCutsWithDefault, which calls 
    // SetDefaultCutValue(defaultCutValue); *G4VUserPhysicsList::SetCuts()*;

    int verboseLevel = mPL->GetVerboseLevel();

    double productionCut = pset.get<double>("physics.minRangeCut");
    double protonProductionCut = pset.get<double>("physics.protonProductionCut");
    mf::LogInfo("GEOM_MINRANGECUT")
      << "Setting production cut to " << productionCut 
      << ", protonProductionCut to " << protonProductionCut << " mm";
    if (pset.get<int>("debug.diagLevel") > 0) {
      G4cout << __func__ << " Setting gamma, e- and e+ production cut"
             << " to " << productionCut << " mm and for proton to "
             << protonProductionCut << " mm" << G4endl;
    }
    //setCutCmd equivalent:
    mPL->SetDefaultCutValue(productionCut);
    mPL->SetCutValue(protonProductionCut, "proton");

    // special cuts per region
    const fhicl::ParameterSet& regionProductionCutsPSet{
      pset.get<fhicl::ParameterSet>("physics.minRangeRegionCuts",fhicl::ParameterSet())};

    if (!regionProductionCutsPSet.is_empty()) {

      const std::vector<std::string> regionNames{regionProductionCutsPSet.get_names()};

      for(const auto& regionName : regionNames) {
        G4Region* region = G4RegionStore::GetInstance()->GetRegion(regionName);
        if (region!=nullptr) {
          G4ProductionCuts* cuts = new G4ProductionCuts();
          double productionCut = regionProductionCutsPSet.get<double>(regionName);
          cuts->SetProductionCut(productionCut); // same cut for gamma, e- and e+, proton/ions
          double protonProductionCut = pset.get<double>("physics.protonProductionCut");
          cuts->SetProductionCut(protonProductionCut,"proton");
          if ( verboseLevel > 0 ) {
            std::cout << __func__ << " Setting gamma, e- and e+ production cut for "
                      << regionName << " to " << productionCut << " mm and for proton to "
                      << protonProductionCut << " mm" << std::endl;
            std::cout << __func__ << " Resulting cuts for gamma, e-, e+, proton: ";
            for (auto const& rcut : cuts->GetProductionCuts() ) {
              std::cout << " " << rcut;
            }
            std::cout << std::endl;
          }
          region->SetProductionCuts(cuts);
        } else {
          if ( verboseLevel > -1 ) {
            std::cout << __func__ << " Did not find requested region: "
                      << regionName << std::endl;
          }
        }
      }
    }
  }
}  // end namespace mu2e
