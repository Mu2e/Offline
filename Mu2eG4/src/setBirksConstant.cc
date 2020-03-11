//
// Set the G4 BirksConstant as specified in the fhicl file.
//
// $Id: setBirksConstant.cc,v 1.2 2015/11/04 22:06:17 genser Exp $
// $Author: genser $
// $Date: 2015/11/04 22:06:17 $
//

#include "Mu2eG4/inc/setBirksConstant.hh"
#include "fhiclcpp/ParameterSet.h"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

// Geant4 includes
#include "G4Material.hh"

// C++ includes

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <utility>

#include "fhiclcpp/ParameterSet.h"

namespace mu2e{

  void setBirksConstant(const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug) {

    fhicl::ParameterSet birksConstsPSet;
    if(phys.BirksConsts.get_if_present(birksConstsPSet)) {

      const std::vector<std::string> matNames{birksConstsPSet.get_names()};

      int verbosityLevel = debug.diagLevel();

      // in principle we could do it without this map and set values direcly from pset
      std::map<std::string,double> birksConstsMap;

      for(const auto& mat: matNames) {
        birksConstsMap[mat] = birksConstsPSet.get<double>(mat);

        if ( verbosityLevel > 0) {
          mf::LogInfo("PHYS")
            << "setting Birks constant for " <<  mat
            << " to " << birksConstsMap[mat] << " mm/MeV";
        }
        G4Material *gmat = findMaterialOrThrow( mat );
        gmat->GetIonisation()->SetBirksConstant(birksConstsMap[mat]*CLHEP::mm/CLHEP::MeV);
      }
    }
  }

}  // end namespace mu2e
