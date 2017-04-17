//
// Set the G4 minimum range/production cut as specified in the config file
//
// Author: KLG
//

#include <iomanip>
#include <iostream>
#include <sstream>

#include "Mu2eG4/inc/setParticleCut.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "G4VUserPhysicsList.hh"

namespace mu2e{


  namespace {

    double getProductionCut(const SimpleConfig& config, std::string const& pName){
      std::ostringstream os;
      os << "physics." << pName << "ProductionCut";
      return config.getDouble(os.str());
    }

    double getProductionCut(const fhicl::ParameterSet& pset, std::string const& pName){
      std::ostringstream os;
      os << "physics." << pName << "ProductionCut";
      return pset.get<double>(os.str());
    }

  }

  void setParticleCut(double theCut,  std::string const& pName, G4VUserPhysicsList* mPL ) {
    mf::LogInfo("GEOM")
      << "Setting "<< pName << " cut to " << theCut << " mm\n";
    // we have e+,e-,gamma range cuts and proton production cut;
    // one can write now the functions for the specific range cuts if needed,
    // but they are covered by the common setMinimumRangeCut for now

    // before doing it one would need to solve the problem with
    // fhicl atoms like e+, so one would need to do a translation from
    // electron to e- etc...; creating another hierarchy level may also be in order

    // the proton production cut which is special
    mPL->SetCutValue(theCut, pName); 

    // proton cut is special, e.g. for list using HP; and needs to be
    // done post initialization and using the above call; it is a
    // production cut also applied to nuclei; see "Setting the cuts"
    // section of the Geant4 User's Guide for Application Developers

  }

  void setParticleCut( SimpleConfig const& config, std::string const& pName, G4VUserPhysicsList* mPL ){
    double pCut = getProductionCut(config,pName);
    setParticleCut(pCut, pName, mPL);
  }

  void setParticleCut(fhicl::ParameterSet const& pset, std::string const& pName, G4VUserPhysicsList* mPL){
    double pCut = getProductionCut(pset,pName);
    setParticleCut(pCut, pName, mPL);
  }


}  // end namespace mu2e
