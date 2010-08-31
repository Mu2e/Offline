//
// Maintain up to date geometry information and serve it to
// other services and to the modules.
//
// $Id: GeometryService.cc,v 1.7 2010/08/31 00:24:51 logash Exp $
// $Author: logash $ 
// $Date: 2010/08/31 00:24:51 $
//
// Original author Rob Kutschke
//

// C++ include files
#include <iostream>
#include <typeinfo>

// Framework include files
#include "DataFormats/Provenance/interface/ModuleDescription.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/Provenance/interface/RunID.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


// Mu2e include files
#include "GeometryService/inc/GeometryService.hh"
#include "TargetGeom/inc/Target.hh"
#include "TargetGeom/inc/TargetMaker.hh"
#include "CTrackerGeom/inc/CTracker.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "LTrackerGeom/inc/LTrackerMaker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TTrackerGeom/inc/TTrackerMaker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "ITrackerGeom/inc/ITrackerMaker.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/CalorimeterMaker.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "BFieldGeom/inc/BFieldManagerMaker.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/BeamlineMaker.hh"
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"
#include "VirtualDetectorGeom/inc/VirtualDetectorMaker.hh"

using namespace std;
using namespace mu2e::calorimeter;

namespace mu2e {

  GeometryService::GeometryService(edm::ParameterSet const& iPS, 
                                   edm::ActivityRegistry&iRegistry) :
    _inputfile(iPS.getUntrackedParameter<std::string>("inputfile","geom000.txt")),
    _detectors(),
    _run_count()
  {
    iRegistry.watchPreBeginRun(this, &GeometryService::preBeginRun);
  }

  GeometryService::~GeometryService(){
  }

  void 
  GeometryService::preBeginRun(edm::RunID const& iID, edm::Timestamp const& iTime) {

    if(++_run_count > 1) {
      edm::LogWarning("GEOM") << "This test version does not change geometry on run boundaries.";
      return;
    }

    edm::LogInfo  log("GEOM");
    log << "Geometry input file is: " << _inputfile << "\n";

    _config = auto_ptr<SimpleConfig>(new SimpleConfig(_inputfile));

    if ( _config->getBool("printConfig",false) ){
      log << *_config;
    }

    // Throw if the configuration 
    checkConfig();

    // Make a detector for every component present in the configuration.
    BeamlineMaker beamlinem( *_config );
    addDetector( beamlinem.getBeamlinePtr() );

    if(_config->getBool("hasTarget",false)){
      TargetMaker targm( *_config );
      addDetector( targm.getTargetPtr() );
    }

    if(_config->getBool("hasCTracker",false)){
      addDetector( std::auto_ptr<CTracker>(new CTracker(*_config)) );
    }

    if(_config->getBool("hasLTracker",false)){
      LTrackerMaker ltm( *_config );
      addDetector( ltm.getLTrackerPtr() );
    } else if (_config->getBool("hasITracker",false)){
      ITrackerMaker itm( *_config );
      addDetector( itm.getITrackerPtr() );
    } else if (_config->getBool("hasTTracker",false)){
      TTrackerMaker ttm( *_config );
      addDetector( ttm.getTTrackerPtr() );
    }

    if(_config->getBool("hasCalorimeter",false)){
      CalorimeterMaker calorm( *_config );
      addDetector( calorm.getCalorimeterPtr() );
    }
    
    if(_config->getBool("hasBFieldManager",false)){
      BFieldManagerMaker bfmgr( *_config);
      addDetector( bfmgr.getBFieldManager() );
    }

    VirtualDetectorMaker vdm( *_config );
    addDetector( vdm.getVirtualDetectorPtr() );

  }

  // Check that the configuration is self consistent.
  void GeometryService::checkConfig(){
    int ntrackers(0);
    string allTrackers;
    if ( _config->getBool("hasLTracker",false) ) {
      allTrackers += " LTracker";
      ++ntrackers;
    }
    if ( _config->getBool("hasTTracker",false) ) {
      allTrackers += " TTracker";
      ++ntrackers;
    }
    if ( _config->getBool("hasCTracker",false) ) {
      allTrackers += " CTracker";
      ++ntrackers;
    }
    if ( _config->getBool("hasITracker",false) ) {
      allTrackers += " ITracker";
      ++ntrackers;
    }

    if ( ntrackers > 1 ){
      throw cms::Exception("GEOM")
        << "This configuration has more than one tracker: "
        << allTrackers
        << "\n";
    }

  }

} // end namespace mu2e
