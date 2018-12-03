//
// Set the G4 minimum range cut as specified in the geometry file.
//
// $Id: setMinimumRangeCut.cc,v 1.2 2012/07/15 22:06:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:17 $
//

// C++ includes
#include <iostream>

// Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// G4 includes
#include "G4VUserPhysicsList.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

// Mu2e includes
#include "Mu2eG4/inc/setMinimumRangeCut.hh"

namespace mu2e{

  void setMinimumRangeCut(const fhicl::ParameterSet& pset, G4VUserPhysicsList* mPL){
    double minRangeCut = pset.get<double>("physics.minRangeCut");
    mf::LogInfo("GEOM_MINRANGECUT")
      << "Setting minRange cut to " << minRangeCut << " mm";
    //setCutCmd equivalent:
    mPL->SetDefaultCutValue(minRangeCut);
    mPL->SetCuts(); // does not work for proton production cut when
                    // called during pre_init or after initialization;
                    // the selective SetCutValue does work

    // special cuts per region
    const fhicl::ParameterSet& minRangeRegionCutsPSet{
      pset.get<fhicl::ParameterSet>("physics.minRangeRegionCuts",fhicl::ParameterSet())};

    if (!minRangeRegionCutsPSet.is_empty()) {

      int verbosityLevel = pset.get<int>("debug.diagLevel",0);
      const std::vector<std::string> regionNames{minRangeRegionCutsPSet.get_names()};

      for(const auto& regionName : regionNames) {
        G4Region* region = G4RegionStore::GetInstance()->GetRegion(regionName);
        if (region!=nullptr) {
          G4ProductionCuts* cuts = new G4ProductionCuts();
          double rangeCut = minRangeRegionCutsPSet.get<double>(regionName);
          cuts->SetProductionCut(rangeCut); // same cut for gamma, e- and e+, proton/ions
          double protonProductionCut = pset.get<double>("physics.protonProductionCut");
          cuts->SetProductionCut(protonProductionCut,"proton");
          if ( verbosityLevel > 0 ) {
            std::cout << __func__ << " Setting gamma, e- and e+ production cut for "
                      << regionName << " to " << rangeCut << " mm and for proton to "
                      << protonProductionCut << " mm" << std::endl;
            std::cout << __func__ << " Resulting cuts for gamma, e-, e+, proton: ";
            for (auto const& rcut : cuts->GetProductionCuts() ) {
              std::cout << " " << rcut;
            }
            std::cout << std::endl;
          }
          region->SetProductionCuts(cuts);
        } else {
          if ( verbosityLevel > -1 ) {
            std::cout << __func__ << " Did not find requested region: "
                      << regionName << std::endl;
          }
        }
      }
    }
  }

}  // end namespace mu2e
