//
// Construct and return an Beamline.
//
//
// $Id: BeamlineMaker.cc,v 1.1 2010/08/31 00:24:51 logash Exp $
// $Author: logash $ 
// $Date: 2010/08/31 00:24:51 $
//
// Original author Peter Shanahan
//

#include <iostream>
#include <iomanip>
#include <cmath>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "BeamlineGeom/inc/BeamlineMaker.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/TransportSolenoid.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

#ifndef __CINT__ 

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  BeamlineMaker::BeamlineMaker(SimpleConfig const& c)
  {
    // Do the real work.
    BuildBeamline(c);
    BuildTS(c);
  }
  
  BeamlineMaker::~BeamlineMaker (){}

  void BeamlineMaker::BuildBeamline(SimpleConfig const& c){

    // Build the global beamline Geometry.

    _beamline = auto_ptr<Beamline>(new Beamline());

    _beamline->_solenoidOffset = c.getDouble("mu2e.solenoidOffset");

  }//::BuildBeamline

  void BeamlineMaker::BuildTS(SimpleConfig const& c){
    
    // Build the TransportSolenoid Geometry.

    // Read or calculate data

    double solenoidOffset = _beamline->_solenoidOffset;

    double rTorus = c.getDouble("toyTS.rTorus");
    double rVac   = c.getDouble("toyTS.rVac");
    double rCryo  = c.getDouble("toyTS.rCryo");

    double ts1HalfLength = c.getDouble("toyTS1.halfLength");
    double ts3HalfLength = solenoidOffset - rTorus;
    double ts5HalfLength = c.getDouble("toyTS5.halfLength");

    double ts1zOffset    = (-rTorus-ts1HalfLength);
    double ts5zOffset    = ( rTorus+ts5HalfLength);

    double coll1HalfLength = c.getDouble("coll1.halfLength");
    double coll3HalfLength = c.getDouble("coll3.halfLength");
    double coll5HalfLength = c.getDouble("coll5.halfLength");
    double coll3Hole       = c.getDouble("coll3.hole");

    // Save results to the geometry objects

    _beamline->_ts._rTorus = rTorus;
    _beamline->_ts._rVac   = rVac;
    _beamline->_ts._rCryo  = rCryo;

    _beamline->_ts._ts1.set(ts1HalfLength,
			    CLHEP::Hep3Vector(solenoidOffset,0.0,ts1zOffset),
			    0 );
    _beamline->_ts._ts3.set(ts3HalfLength,
			    CLHEP::Hep3Vector(0.0,0.0,0.0),
			    new CLHEP::HepRotation(CLHEP::HepRotationY(90.0*CLHEP::degree)) );
    _beamline->_ts._ts5.set(ts5HalfLength,
			    CLHEP::Hep3Vector(-solenoidOffset,0.0,ts5zOffset),
			    0 );

    _beamline->_ts._coll1.set( coll1HalfLength,CLHEP::Hep3Vector(0.0,0.0,-ts1HalfLength+coll1HalfLength));
    _beamline->_ts._coll31.set(coll3HalfLength,CLHEP::Hep3Vector(0.0,0.0,-coll3HalfLength-coll3Hole/2));
    _beamline->_ts._coll32.set(coll3HalfLength,CLHEP::Hep3Vector(0.0,0.0, coll3HalfLength+coll3Hole/2));
    _beamline->_ts._coll5.set( coll5HalfLength,CLHEP::Hep3Vector(0.0,0.0, ts5HalfLength-coll5HalfLength));

  }

} // namespace mu2e

#endif
