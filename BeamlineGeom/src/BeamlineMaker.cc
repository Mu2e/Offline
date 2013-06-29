//
// Construct and return an Beamline.
//
//
// $Id: BeamlineMaker.cc,v 1.13 2013/06/29 00:04:29 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/06/29 00:04:29 $
//
// Original author Peter Shanahan
//

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/BeamlineMaker.hh"
#include "BeamlineGeom/inc/TransportSolenoid.hh"
#include "BeamlineGeom/inc/Collimator_TS1.hh"
#include "BeamlineGeom/inc/Collimator_TS3.hh"
#include "BeamlineGeom/inc/Collimator_TS5.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

#ifndef __CINT__

namespace mu2e {

  using CLHEP::Hep3Vector;
  using CLHEP::HepRotation;
  using CLHEP::HepRotationX;
  using CLHEP::HepRotationY;
  using CLHEP::degree;

  std::unique_ptr<Beamline> BeamlineMaker::make(const SimpleConfig& c) {
    std::unique_ptr<Beamline> res ( new Beamline() );
    BuildBeamline     (c, res.get() );
    BuildTSCryostat   (c, res.get() );
    BuildTSCoils      (c, res.get() );
    BuildTSCollimators(c, &res->_ts );
    BuildPbarWindow   (c, &res->_ts );
    return res;
  }

  void BeamlineMaker::BuildBeamline(const SimpleConfig& c, Beamline* bl) {
    bl->_solenoidOffset = c.getDouble("mu2e.solenoidOffset");
  }

  void BeamlineMaker::BuildTSCryostat(const SimpleConfig& c, Beamline* bl ){

    TransportSolenoid & ts ( bl->_ts );
    ts._rTorus = c.getDouble("ts.rTorus",0.);
    ts._rVac   = c.getDouble("ts.rVac",0.); 
    ts._material = c.getString("ts.materialName");
    ts._insideMaterial = c.getString("ts.insideMaterialName");

    double ts1HalfLength = c.getDouble("ts.ts1.halfLength");
    double ts3HalfLength = bl->solenoidOffset() - ts.torusRadius();
    double ts5HalfLength = c.getDouble("ts.ts5.halfLength");

    double ts1zOffset    = (-ts.torusRadius()-ts1HalfLength);
    double ts5zOffset    = ( ts.torusRadius()+ts5HalfLength);

    std::vector<TSSection*> tmp_inOut;
    // Straight sections
    ts._ts1in.set( ts.innerRadius(), 
                   c.getDouble("ts.ts1in.rOut",0.),
                   ts1HalfLength,
                   Hep3Vector( bl->solenoidOffset(),0.0,ts1zOffset),
                   nullptr );

    ts._ts1out.set( c.getDouble("ts.ts1out.rIn",0.),
                    c.getDouble("ts.ts1out.rOut",0.),
                    ts1HalfLength,
                    Hep3Vector( bl->solenoidOffset(),0.0,ts1zOffset),
                    nullptr );
    
    tmp_inOut.push_back( dynamic_cast<TSSection*>( &ts._ts1in ) );
    tmp_inOut.push_back( dynamic_cast<TSSection*>( &ts._ts1out ) );
    ts._cryoMap.insert( make_pair( TransportSolenoid::TS1, std::move( tmp_inOut ) ) );

    ts._ts3in.set( ts.innerRadius(), 
                   c.getDouble("ts.ts3in.rOut",0.),
                   ts3HalfLength,
                   Hep3Vector(),
                   new HepRotation(HepRotationY(90.0*degree)) );

    ts._ts3out.set( c.getDouble("ts.ts3out.rIn",0.),
                    c.getDouble("ts.ts3out.rOut",0.),
                    ts3HalfLength,
                    Hep3Vector(),
                    new HepRotation(HepRotationY(90.0*degree)) );

    tmp_inOut.push_back( dynamic_cast<TSSection*>( &ts._ts3in ) );
    tmp_inOut.push_back( dynamic_cast<TSSection*>( &ts._ts3out ) );
    ts._cryoMap.insert( make_pair( TransportSolenoid::TS3, std::move( tmp_inOut ) ) );

    ts._ts5in.set( ts.innerRadius(), 
                   c.getDouble("ts.ts5in.rOut",0.),
                   ts5HalfLength,
                   Hep3Vector(-bl->solenoidOffset(),0.0,ts5zOffset),
                   nullptr );

    ts._ts5out.set( c.getDouble("ts.ts5out.rIn",0.),
                    c.getDouble("ts.ts5out.rOut",0.),
                    ts5HalfLength,
                    Hep3Vector(-bl->solenoidOffset(),0.0,ts5zOffset),
                    nullptr );

    tmp_inOut.push_back( dynamic_cast<TSSection*>( &ts._ts5in ) );
    tmp_inOut.push_back( dynamic_cast<TSSection*>( &ts._ts5out ) );
    ts._cryoMap.insert( make_pair( TransportSolenoid::TS5, std::move( tmp_inOut ) ) );

    // Torus sections
    ts._ts2in.set( ts.torusRadius(),
                   ts.innerRadius(),
                   c.getDouble("ts.ts2in.rOut",0.),
                   1.5*CLHEP::pi, CLHEP::halfpi,
                   Hep3Vector( ts.getTS3_in().getHalfLength(), 0.,-ts.torusRadius() ),
                   new HepRotation(HepRotationX(90.0*degree)) );

    ts._ts2out.set( ts.torusRadius(),
                    c.getDouble("ts.ts2out.rIn",0.),
                    c.getDouble("ts.ts2out.rOut",0.),
                    1.5*CLHEP::pi, CLHEP::halfpi,
                    Hep3Vector( ts.getTS3_in().getHalfLength(), 0.,-ts.torusRadius() ),
                    new HepRotation(HepRotationX(90.0*degree)) );

    tmp_inOut.push_back( dynamic_cast<TSSection*>( &ts._ts2in ) );
    tmp_inOut.push_back( dynamic_cast<TSSection*>( &ts._ts2out ) );
    ts._cryoMap.insert( make_pair( TransportSolenoid::TS2, std::move( tmp_inOut ) ) );

    ts._ts4in.set( ts.torusRadius(),
                   ts.innerRadius(),
                   c.getDouble("ts.ts4in.rOut",0.),
                   CLHEP::halfpi, CLHEP::halfpi,
                   Hep3Vector( -ts.getTS3_in().getHalfLength(), 0.,ts.torusRadius() ),
                   new HepRotation(HepRotationX(90.0*degree)) );
    
    ts._ts4out.set( ts.torusRadius(),
                    c.getDouble("ts.ts4out.rIn",0.),
                    c.getDouble("ts.ts4out.rOut",0.),
                    CLHEP::halfpi, CLHEP::halfpi,
                    Hep3Vector( -ts.getTS3_in().getHalfLength(), 0.,ts.torusRadius() ),
                    new HepRotation(HepRotationX(90.0*degree)) );

    tmp_inOut.push_back( dynamic_cast<TSSection*>( &ts._ts4in ) );
    tmp_inOut.push_back( dynamic_cast<TSSection*>( &ts._ts4out ) );
    ts._cryoMap.insert( make_pair( TransportSolenoid::TS4, std::move( tmp_inOut ) ) );

    // - end wall parameters
    ts._rIn_endWallU1 = c.getDouble("ts.tsUendWall1.rIn",0.); 
    ts._rIn_endWallU2 = c.getDouble("ts.tsUendWall2.rIn",0.); 
    ts._rIn_endWallD  = c.getDouble("ts.tsUendWall.rIn",0.);  

    ts._rOut_endWallU1 = c.getDouble("ts.tsUendWall1.rOut",0.); 
    ts._rOut_endWallU2 = c.getDouble("ts.tsUendWall2.rOut",0.); 
    ts._rOut_endWallD  = c.getDouble("ts.tsUendWall.rOut",0.);  

    ts._halfLength_endWallU1 = c.getDouble("ts.tsDendWall1.halfLength",0.); 
    ts._halfLength_endWallU2 = c.getDouble("ts.tsDendWall2.halfLength",0.); 
    ts._halfLength_endWallD  = c.getDouble("ts.tsDendWall.halfLength",0.);  

  }

  void BeamlineMaker::BuildTSCoils (const SimpleConfig& c, Beamline* bl ) {

    TransportSolenoid * ts ( &bl->_ts );

    ts->_coilMaterial = c.getString("ts.coils.material");

    // Loop over TS regions
    for ( int iTS  = TransportSolenoid::TS1 ; iTS <= TransportSolenoid::TS5 ; iTS++ )
      { 

        std::vector<double> tmp_rIn, tmp_rOut, tmp_zLength, tmp_yRotAngle;

        std::ostringstream prefix; 
        prefix << "ts" << iTS+1;
        c.getVectorDouble( prefix.str()+".coils.rIn"      , tmp_rIn       );
        c.getVectorDouble( prefix.str()+".coils.rOut"     , tmp_rOut      );
        c.getVectorDouble( prefix.str()+".coils.zLength"  , tmp_zLength   );
        c.getVectorDouble( prefix.str()+".coils.yRotAngle", tmp_yRotAngle );

        TransportSolenoid::checkSizeOfVectors( ts->getNCoils(iTS),
                                               { tmp_rIn, tmp_rOut, tmp_zLength, tmp_yRotAngle } );

        TransportSolenoid::vector_unique_ptr<Coil> tmp_coilvector;

        const double totCoilLength = accumulate( tmp_zLength.begin(), tmp_zLength.end(), 0.);

        // Loop over coils per TS region
        for ( unsigned i(0) ; i < ts->getNCoils(iTS) ; i++ ) {
          std::unique_ptr<Coil> coil ( new Coil( tmp_rIn.at(i), 
                                                 tmp_rOut.at(i),
                                                 tmp_zLength.at(i)*0.5 ) );
          coil->setRotation( new HepRotation(HepRotationY(-tmp_yRotAngle.at(i)*degree) ) );
          coil->setOrigin  ( ts->getNCoils(iTS), 
                             ts->getTSCryo(iTS,TransportSolenoid::IN),
                             i, coil->getRotation(), totCoilLength );

          tmp_coilvector.push_back( std::move( coil ) );
        }

        ts->_coilMap.insert( std::make_pair( iTS, std::move( tmp_coilvector ) ) );

      }

  }

  void BeamlineMaker::BuildTSCollimators (const SimpleConfig& c, TransportSolenoid* ts) {

    // Set collimators
    double coll1HalfLength = c.getDouble("coll1.halfLength");
    double coll3HalfLength = c.getDouble("coll3.halfLength");
    double coll5HalfLength = c.getDouble("coll5.halfLength");
    double coll3Hole       = c.getDouble("coll3.hole");

    CollimatorTS1 & coll1  ( ts->_coll1 );
    CollimatorTS3 & coll31 ( ts->_coll31 );
    CollimatorTS3 & coll32 ( ts->_coll32 );
    CollimatorTS5 & coll5  ( ts->_coll5 );

    coll1 .set(coll1HalfLength,Hep3Vector(0.,0.,0.));
    coll31.set(coll3HalfLength,Hep3Vector(0.,0.,-coll3HalfLength-coll3Hole/2));
    coll32.set(coll3HalfLength,Hep3Vector(0.,0., coll3HalfLength+coll3Hole/2));
    coll5 .set(coll5HalfLength,Hep3Vector(0.,0.,0.));

    // TS1
    coll1._rIn1      = c.getDouble("coll1.innerRadius1",0.);
    coll1._rIn2      = c.getDouble("coll1.innerRadius2",0.);
    coll1._rIn3      = c.getDouble("coll1.innerRadius3",0.);
    coll1._material1 = c.getString("coll1.material1Name");
    coll1._material2 = c.getString("coll1.material2Name");

    // TS3
    coll32._hole             = coll31._hole              = c.getDouble("coll3.hole",0.);
    coll32._holeRadius       = coll31._holeRadius        = c.getDouble("coll3.holeRadius",0.);
    coll32._holeHalfHeight   = coll31._holeHalfHeight    = c.getDouble("coll3.holeHalfHeight",0.);
    coll32._holeDisplacement = coll31._holeDisplacement  = c.getDouble("coll3.holeDisplacement",0.);
    coll32._rotationAngle    = coll31._rotationAngle     = c.getDouble("coll3.rotationAngle",0.);
    coll32._material         = coll31._material          = c.getString("coll3.materialName");

    // TS5
    coll5._rIn         = c.getDouble("coll5.innerRadius");
    coll5._rMid1       = c.getDouble("coll5.midRadius1");
    coll5._rMid2       = c.getDouble("coll5.midRadius2");
    coll5._rOut        = c.getDouble("coll5.outerRadius");
    coll5._halfLengthU = c.getDouble("coll5.halfLengthU");
    coll5._halfLengthD = c.getDouble("coll5.halfLengthD");

    coll5._material    = c.getString("coll5.materialName");
    coll5._absMaterial = c.getString("coll5.absorberMaterialName");

  }

  void BeamlineMaker::BuildPbarWindow(const SimpleConfig& c, TransportSolenoid* ts){

    PbarWindow & pbarWindow ( ts->_pbarWindow );

    pbarWindow._shape    = c.getString("pbar.Type","disk");
    pbarWindow._material = c.getString("pbar.materialName");
    pbarWindow.set( c.getDouble("pbar.halfLength"),
                    Hep3Vector() );
    pbarWindow._rOut     = ts->innerRadius();

    // Parameters for wedge
    pbarWindow._y0       = c.getDouble("pbarwedge.y0",0.);
    pbarWindow._y1       = c.getDouble("pbarwedge.y1",0.);
    pbarWindow._dz0      = c.getDouble("pbarwedge.dz0",0.);
    pbarWindow._dz1      = c.getDouble("pbarwedge.dz1",0.);

  }

} // namespace mu2e

#endif
