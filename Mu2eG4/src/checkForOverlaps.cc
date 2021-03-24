//
// Free function to do Geant4 overlap check
//
//
// Original author KLG
//

// C++ includes
#include <iostream>
#include <iomanip>

// Framework includes
#include "cetlib/cpu_timer.h"

// Mu2e includes
#include "Mu2eG4/inc/checkForOverlaps.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// G4 includes
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4VSolid.hh"

namespace mu2e{

  bool checkForOverlaps( G4VPhysicalVolume* const pv,
                         SimpleConfig const& config, bool verbose){

    // verbose = true; // tmp override

    static G4double const minSurfaceCheckPoints =      config.getInt("g4.minSurfaceCheckPoints",     100);
    static G4double const maxSurfaceCheckPoints =      config.getInt("g4.maxSurfaceCheckPoints",10000000);
    static G4double const nSurfaceCheckPointsPercmsq = config.getInt("g4.nSurfaceCheckPointsPercmsq",  1);
    // make sure it is > 0
    static G4int const iSurfaceCheckPointsPercmsq  = std::max(1.,nSurfaceCheckPointsPercmsq);
    static G4double const mmsqPerPoint = 100./iSurfaceCheckPointsPercmsq;

    G4double vsa = pv->GetLogicalVolume()->GetSolid()->GetSurfaceArea();

    G4int nSurfaceCheckPoints = std::max(std::min(vsa/mmsqPerPoint, maxSurfaceCheckPoints),
                                         minSurfaceCheckPoints);

    cet::cpu_timer ocTimer;

    if (verbose) {
      auto cflags = G4cout.flags();
      G4cout << __func__ << " SurfaceArea is: " << std::setw(12) << G4long(vsa) 
             << ", surfaceCheckPoints " << std::setw(12) << nSurfaceCheckPoints 
             << " for " << pv->GetName() << std::endl;
      G4cout.flags(cflags);
      ocTimer.start();
    }

    bool ocRC = pv->CheckOverlaps( nSurfaceCheckPoints, 0.0, true );

    if (verbose) {
      ocTimer.stop();
      double realElapsed = ocTimer.elapsed_real_time();
      double cpuElapsed  = ocTimer.elapsed_cpu_time();
      auto cflags = G4cout.flags();
      auto cprecision = G4cout.precision();
      auto cwidth = G4cout.width();
      G4cout << __func__
	     << " RealElapsed: " << std::fixed << std::setw(12) << std::setprecision(6) << realElapsed
	     << " CpuElapsed:  " << cpuElapsed
	     << " CpuElapsed/surfaceCheckPoint: " << std::setprecision(9) << cpuElapsed/nSurfaceCheckPoints
	     << std::endl;
      G4cout.flags(cflags);
      G4cout.precision(cprecision);
      G4cout.width(cwidth);
    }

    return ocRC;

  }

}  // end namespace mu2e

