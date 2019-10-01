//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eStudyWorld.cc,v 1.8 2013/04/05 20:22:45 genser Exp $
// $Author: genser $
// $Date: 2013/04/05 20:22:45 $
//
// Original author K. Genser based on Mu2eWorld
//
//  Heirarchy is:
//  0      World (vaccum?)
//  1      volume created in the construct... function
//

// C++ includes
#include <iostream>
#include <vector>
#include <cmath>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "art/Utilities/make_tool.h"

// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
// #include "Mu2eG4/inc/constructStudyEnv_v001.hh"
// #include "Mu2eG4/inc/constructStudyEnv_v002.hh"
// #include "Mu2eG4/inc/constructStudyEnv_v003.hh"
// #include "Mu2eG4/inc/constructStudyEnv_v004.hh"
#include "Mu2eG4/inc/Mu2eStudyWorld.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "BFieldGeom/inc/BFieldManager.hh"

// G4 includes
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4ClassicalRK4.hh"
#include "G4ImplicitEuler.hh"
#include "G4ExplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4GDMLParser.hh"

#include "Mu2eG4/inc/Mu2eGlobalField.hh"
#include "Mu2eG4/inc/FieldMgr.hh"


using namespace std;

namespace mu2e {

  Mu2eStudyWorld::Mu2eStudyWorld()
  {}

  Mu2eStudyWorld::~Mu2eStudyWorld(){
    // Do not destruct the solids, logical volumes or physical volumes.
    // G4 looks after that itself.
  }

  Mu2eStudyWorld::Mu2eStudyWorld(const fhicl::ParameterSet& pset,
				       SensitiveDetectorHelper *sdHelper/*no ownership passing*/)
    : sdHelper_(sdHelper)
    , pset_(pset)
    , writeGDML_(pset.get<bool>("debug.writeGDML"))
    , gdmlFileName_(pset.get<std::string>("debug.GDMLFileName"))
    , g4stepperName_(pset.get<std::string>("physics.stepper"))
    , bfieldMaxStep_(pset.get<double>("physics.bfieldMaxStep"))//unused
  {
    _verbosityLevel = pset.get<int>("debug.worldVerbosityLevel");
  }

  // This is the callback called by G4
  G4VPhysicalVolume * Mu2eStudyWorld::construct(){

    // Construct all of the world

    _verbosityLevel =  max(_verbosityLevel,_config.getInt("world.verbosityLevel", 0));

    // we will only use very few elements of the geometry service;
    // mainly its config (SimpleConfig) part and the origin which is
    // used by VolumeInfo and related fuctions/classes

    MaterialFinder materialFinder(_config);

    // Dimensions and material of the world.
    G4Material* worldMaterial = materialFinder.get("world.materialName");

    const bool worldBoxVisible = _config.getBool("world.boxVisible");
    const bool worldBoxSolid   = _config.getBool("world.boxSolid");

    const bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck");
    const bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible");
    const bool placePV             = true;

    G4double worldHalfLength       = _config.getDouble("world.halfLength");
    G4double outerLayerThickness   = _config.getDouble("world.outerLayerThickness");

    vector<double> worldBoundaries(3,worldHalfLength+outerLayerThickness);

    // the canonical World volume
    VolumeInfo worldVInfo(nestBox("World", 
                                  worldBoundaries,
                                  worldMaterial, 
                                  0, 
                                  G4ThreeVector(),
                                  0, // no parent
                                  0, // we assign this volume copy number 0
                                  worldBoxVisible, 
                                  G4Colour::Cyan(), 
                                  worldBoxSolid,
                                  forceAuxEdgeVisible, 
                                  placePV, 
                                  false)); // do not surface check this one

    // Now box almost filling up the world to force a step close to
    // the world boundary

    vector<double> boxInWorldBoundaries(3,worldHalfLength);

    VolumeInfo boxInTheWorldVInfo(nestBox("BoxInTheWorld", 
                                          boxInWorldBoundaries,
                                          worldMaterial, 
                                          0, // no rotation
                                          G4ThreeVector(),
                                          worldVInfo,
                                          1, // we assign this volume copy number 1
                                          worldBoxVisible, 
                                          G4Colour::Cyan(), 
                                          worldBoxSolid,
                                          forceAuxEdgeVisible, 
                                          placePV, 
                                          doSurfaceCheck));

    std::string simulatedDetector = _geom.pset().get<std::string>("simulatedDetector.tool_type");

    if (simulatedDetector != "Mu2e") {
    
      constructEnv_ = art::make_tool<InitEnvToolBase>(_geom.pset().get<fhicl::ParameterSet>("simulatedDetector"));

      if (constructEnv_) constructEnv_->construct(boxInTheWorldVInfo,_config);
      else {
	throw cet::exception("CONFIG") << __func__ << ": unknown study environment: " << simulatedDetector << "\n";
      }
    }
    if ( _verbosityLevel > 0) {
      cout << __func__ << " world half dimensions     : " 
           << worldBoundaries[0] << ", "
           << worldBoundaries[1] << ", "
           << worldBoundaries[2] << ", "
           << endl;
    }

    // if      ( seVer == 1 ) constructStudyEnv_v001(boxInTheWorldVInfo, _config);
    // else if ( seVer == 2 ) constructStudyEnv_v002(boxInTheWorldVInfo, _config);
    // else if ( seVer == 3 ) constructStudyEnv_v003(boxInTheWorldVInfo, _config);
    // else if ( seVer == 4 ) constructStudyEnv_v004(boxInTheWorldVInfo, _config);

    //    sdHelper_->instantiateLVSDs(_config); // needs work in the study case

    // Create magnetic fields and managers only after all volumes have been defined.
    //    constructBFieldAndManagers();
    //    constructStepLimiters();

    // Write out geometry into a gdml file.
    if (writeGDML_) {
      G4GDMLParser parser;
      parser.Write(gdmlFileName_, worldVInfo.logical);
    }

    return worldVInfo.physical;
  }
    

    
    void Mu2eStudyWorld::constructSDandField(){
        
        std::cout << "We are in Mu2eStudyWorld::constructSDandField()" << std::endl;
        
        //constructWorldSD();
    }


  // Adding a step limiter is a two step process.
  // 1) In the physics list constructor add a G4StepLimiter to the list of discrete
  //    physics processes attached to each particle species of interest.
  //
  // 2) In this code, create a G4UserLimits object and attach it to the logical
  //    volumes of interest.
  // The net result is specifying a step limiter for pairs of (logical volume, particle species).
  //
  void Mu2eStudyWorld::constructStepLimiters() {

    // Maximum step length, in mm.

    // AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
    // G4UserLimits* stepLimit = reg.add( G4UserLimits(bfieldMaxStep_) );

    // double maxStep = _config.getDouble("bfield.maxStep", 20.);

    // We may make separate G4UserLimits objects per logical volume but we may choose not to.
    // _stepLimits.push_back( G4UserLimits(maxStep) );
    // G4UserLimits* stepLimit = &(_stepLimits.back());

    // AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
    // G4UserLimits* stepLimit = reg.add( G4UserLimits(maxStep) );
    // tracker->SetUserLimits( stepLimit );

  } // end Mu2eStudyWorld::constructStepLimiters(){

} // end namespace mu2e
