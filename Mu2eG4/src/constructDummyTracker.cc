//
// Free function to construct a placeholder for the tracker.
// Useful for some low detail graphics.
//
// $Id: constructDummyTracker.cc,v 1.6 2012/07/15 22:06:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:17 $
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
#include "Mu2eG4/inc/constructDummyTracker.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestTubs.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"


using namespace std;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructDummyTracker( G4LogicalVolume* mother,
                                    double zOff,
                                    SimpleConfig const& config ){

    // A helper class for parsing the config file.
    MaterialFinder materialFinder(config);

    // Parse the config file.
    double rIn           = config.getDouble("dummytracker.rIn",            0.);
    double rOut          = config.getDouble("dummytracker.rOut",         800.);
    double halfLength    = config.getDouble("dummytracker.halfLength",  1300.);
    double z0            = config.getDouble("dummytracker.z0" ,        10200.);
    G4Material* material = materialFinder.get("dummytracker.materialName","WAGVacuum");

    bool doSurfaceCheck = config.getBool("g4.doSurfaceCheck",false);

    // Parameters of a G4Tubs.
    TubsParams params(rIn, rOut, halfLength);

    // Position of the tracker within its mother volume.
    G4ThreeVector trackerOffset(0.,0.,z0-zOff);

    VolumeInfo info = nestTubs( "TrackerMother",
                                params,
                                material,
                                0,
                                trackerOffset,
                                mother,
                                0,
                                true,
                                G4Color::Red(),
                                false,
                                true,
                                true,
                                doSurfaceCheck
                                );
    return info;
  }

}  // end namespace mu2e

