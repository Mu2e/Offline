//
//  checkConfigRelics.cc
//  
//  Checks the configuration of a geomtry file.  Throws an exception
//  if an "old style" geomtry file that is no longer supported is used.
//
//  Author: Lisa Goodenough
//  Date: 2017/4/19
//
//

// C++ includes
#include <vector>
#include <string>

// framework includes
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/checkConfigRelics.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {
    
    void checkConfigRelics(const SimpleConfig& config) {
        
        static const std::vector<std::string> keys = {
            
            // G4_module
            "g4.printPhysicsProcessSummary",
            "g4.pointTrajectoryMinSteps",
            
            // old post G4InitializeTasks() call tree
            "g4.PiENuPolicy",
            "g4.PiENuPolicyVerbosity",
            "g4.minRangeCut",
            "g4.noDecay",
            "g4.doMuMinusConversionAtRest",
            "g4.useNewMuMinusAtomicCapture",
            
            // physicsListDecider() call tree
            "g4.physicsListName",
            "g4.useNewMuMinusAtomicCapture",
            "g4.decayMuonsWithSpin",
            
            // Mu2eWorld
            "world.verbosityLevel",
            "tracker.ActiveWr_Wl_SD",
            "writeGDML",
            "GDMLFileName",
            "g4.stepper",
            "bfield.maxStep",
            
            // ConstructMaterials
            "mu2e.standardDetector",
            "g4.printElements",
            "g4.printMaterials",
            
            // old StackingAction
            "g4.doCosmicKiller",
            "g4.cosmicKillLevel",
            "g4.cosmicVerbose",
            "g4.cosmicPcut",
            "g4.yaboveDirtYmin",
            "g4.stackPrimaryOnly",
            "g4.killLowEKine",
            "g4.killPitchToLowToStore",
            "g4.minPitch",
            "g4.stackingActionDropPDG",
            "g4.stackingActionKeepPDG",
            "g4.eKineMin",
            "g4.killLowEKinePDG",
            "g4.eKineMinPDG",
            
            // old TrackingAction
            "g4.particlesSizeLimit",
            "g4.mcTrajectoryMomentumCut",
            "g4.saveTrajectoryMomentumCut",
            "g4.mcTrajectoryMinSteps",
            "g4.printTrackTiming",
            "g4.trackingActionEventList",
            
            // old SteppingAction
            "g4.steppingActionStepsSizeLimit",
            "g4.killLowEKine",
            "g4SteppingAction.killInTheseVolumes",
            "g4SteppingAction.killerVerbose",
            "g4SteppingAction.killInHallAir",
            "g4.eKineMin",
            "g4.killLowEKinePDG",
            "g4.eKineMinPDG",
            "g4.steppingActionEventList",
            "g4.steppingActionTrackList",
            "g4.steppingActionMaxSteps",
            "g4.steppingActionMaxGlobalTime",
            "g4.steppingActionTimeVD",
            "g4.mcTrajectoryVolumes",
            "g4.mcTrajectoryVolumePtDistances",
            "g4.mcTrajectoryDefaultMinPointDistance"
            
        };
        
        std::string present;
        for(const auto k: keys) {
            if(config.hasName(k)) {
                present += k+" ";
            }
        }
        if(!present.empty()) {
            throw cet::exception("CONFIG") << "Please use fcl to configure Mu2eG4_module. "
            << "Detected obsolete SimpleConfig parameters: " << present;
        }
    }
    
} // end namespace mu2e

