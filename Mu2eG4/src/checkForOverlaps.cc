//
// Free function to do Geant4 overlap check
//
// $Id: checkForOverlaps.cc,v 1.1 2013/08/30 15:52:33 genser Exp $
// $Author: genser $
// $Date: 2013/08/30 15:52:33 $
//
// Original author KLG
//

// C++ includes
#include <iostream>
#include <iomanip>

// Mu2e includes
#include "Mu2eG4/inc/checkForOverlaps.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// G4 includes
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

namespace mu2e{

  bool checkForOverlaps( G4VPhysicalVolume* const pv,
                         SimpleConfig const& config, bool verbose){

    static G4double const minSurfaceCheckPoints =      config.getInt("g4.minSurfaceCheckPoints",     100);
    static G4double const maxSurfaceCheckPoints =      config.getInt("g4.maxSurfaceCheckPoints",10000000);
    static G4double const nSurfaceCheckPointsPercmsq = config.getInt("g4.nSurfaceCheckPointsPercmsq",  1);
    // make sure it is > 0
    static G4int const iSurfaceCheckPointsPercmsq  = std::max(1.,nSurfaceCheckPointsPercmsq);
    static G4double const mmsqPerPoint = 100./iSurfaceCheckPointsPercmsq;

    G4double vsa = pv->GetLogicalVolume()->GetSolid()->GetSurfaceArea();

    G4int nSurfaceCheckPoints = std::max(std::min(vsa/mmsqPerPoint, maxSurfaceCheckPoints),
                                         minSurfaceCheckPoints);

    if (verbose) {
      G4cout << __func__ << " SurfaceArea is: " << std::setw(12) << G4long(vsa) 
             << ", surfaceCheckPoints " << std::setw(12) << nSurfaceCheckPoints 
             << " for " << pv->GetName() << std::endl;
    }

    return pv->CheckOverlaps( nSurfaceCheckPoints, 0.0, true );

  }

}  // end namespace mu2e

