//
// Free function to construct a placeholder for the stopping target.
// Useful for some low detail graphics.
//
//
// Original author Rob Kutschke
//

// C++ includes
#include <iostream>
#include <string>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/constructDummyStoppingTarget.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestTubs.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4ThreeVector.hh"


using namespace std;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructDummyStoppingTarget( VolumeInfo   const& mother,
                                           SimpleConfig const& config ){

    // A helper class for parsing the config file.
    MaterialFinder materialFinder(config);

    // Parse the config file.
    double rIn           = config.getDouble("dummyStoppingTarget.rIn",           0.);
    double rOut          = config.getDouble("dummyStoppingTarget.rOut",        100.);
    double halfLength    = config.getDouble("dummyStoppingTarget.halfLength",  400.);
    double z0            = config.getDouble("dummyStoppingTarget.z0" ,        5900.);
    G4Material* material = materialFinder.get("dummyStoppingTarget.materialName","WAGVacuum");

    bool doSurfaceCheck = config.getBool("g4.doSurfaceCheck",false);

    // Parameters of a G4Tubs.
    TubsParams params(rIn, rOut, halfLength );

    // Position of the tracker within its mother volume
    G4ThreeVector offset(G4ThreeVector(0.,0.,z0) - mother.centerInMu2e());

    VolumeInfo info = nestTubs( "StoppingTargetMother",
                                params,
                                material,
                                0,
                                offset,
                                mother,
                                0,
                                config.getBool("stoppingTarget.visible",true),
                                G4Color::Yellow(),
                                config.getBool("stoppingTarget.solid",true),
                                config.getBool("g4.forceAuxEdgeVisible",false),
                                true,
                                doSurfaceCheck
                                );

    return info;
  }

}  // end namespace mu2e

