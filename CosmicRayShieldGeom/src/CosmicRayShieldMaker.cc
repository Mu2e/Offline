//
// Construct and return CosmicRayShield
//
// $Id: CosmicRayShieldMaker.cc,v 1.1 2011/01/25 16:43:52 genser Exp $
// $Author: genser $ 
// $Date: 2011/01/25 16:43:52 $
//
// Original author KLG based on Rob Kutschke's construct... functions
//

// c++ includes
#include <iostream>
#include <iomanip>
#include <cmath>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "CosmicRayShieldGeom/inc/CosmicRayShieldMaker.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BeamlineGeom/inc/Beamline.hh"

#include "TrackerGeom/inc/TubsParams.hh"

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  CosmicRayShieldMaker::CosmicRayShieldMaker(SimpleConfig const& _config)
  {

    int static const diagLevel = 0;

    _crs = auto_ptr<CosmicRayShield>(new CosmicRayShield());

    if( ! _config.getBool("hasCosmicRayShield",false) ) return;

    // first make the steel (fluxreturn) (based on the old constructSteel)

    // Extract information from the config file
    double CRSSteelShieldHalfThick = _config.getDouble("fluxcrv.HallSteelHalfThick");
    double CRSSteelShieldHalfLenXY = _config.getDouble("fluxcrv.HallSteelHalfLengthXY");
    double CRSSteelShieldHalfLenZ  = _config.getDouble("fluxcrv.HallSteelHalfLengthZ");

    // Compute dimensions of 5 sides in Mu2e coordinates
    double CRSSteelShieldTopHalfX   = CRSSteelShieldHalfLenXY + CRSSteelShieldHalfThick;
    double CRSSteelShieldTopHalfY   = CRSSteelShieldHalfThick;
    double CRSSteelShieldTopHalfZ   = CRSSteelShieldHalfLenZ;
    double CRSSteelShieldSideHalfX  = CRSSteelShieldHalfThick;
    double CRSSteelShieldSideHalfY  = CRSSteelShieldHalfLenXY - CRSSteelShieldHalfThick;
    double CRSSteelShieldSideHalfZ  = CRSSteelShieldHalfLenZ; 
    double CRSSteelShieldFrontHalfX = CRSSteelShieldHalfLenXY - CRSSteelShieldHalfThick;
    double CRSSteelShieldFrontHalfY = CRSSteelShieldHalfLenXY - CRSSteelShieldHalfThick;
    double CRSSteelShieldFrontHalfZ = CRSSteelShieldHalfThick;

    double CRSSteelShieldTopDims[3] ={
      CRSSteelShieldTopHalfX,
      CRSSteelShieldTopHalfY,
      CRSSteelShieldTopHalfZ
    };

    double CRSSteelShieldSideDims[3] ={
      CRSSteelShieldSideHalfX,
      CRSSteelShieldSideHalfY,
      CRSSteelShieldSideHalfZ
    };

    double CRSSteelShieldFrontDims[3] ={
      CRSSteelShieldFrontHalfX,
      CRSSteelShieldFrontHalfY,           
      CRSSteelShieldFrontHalfZ                        
    };

    if ( diagLevel > 0) {

      CLHEP::Hep3Vector CRSSteelShieldTopDimsV(
                                          CRSSteelShieldTopHalfX,
                                          CRSSteelShieldTopHalfY,
                                          CRSSteelShieldTopHalfZ
                                          );

      CLHEP::Hep3Vector CRSSteelShieldSideDimsV(
                                           CRSSteelShieldSideHalfX,
                                           CRSSteelShieldSideHalfY,
                                           CRSSteelShieldSideHalfZ
                                           );

      CLHEP::Hep3Vector CRSSteelShieldFrontDimsV(
                                            CRSSteelShieldFrontHalfX,
                                            CRSSteelShieldFrontHalfY,           
                                            CRSSteelShieldFrontHalfZ                        
                                            );

      cout << __func__ << " CRSSteelShieldTopDimsV   : " << CRSSteelShieldTopDimsV   << endl;
      cout << __func__ << " CRSSteelShieldSideDimsV  : " << CRSSteelShieldSideDimsV  << endl;
      cout << __func__ << " CRSSteelShieldFrontDimsV : " << CRSSteelShieldFrontDimsV << endl;

    }

    // Get positions of each side. Assuming view from target foils 
    double dsCoilZ0          = _config.getDouble("toyDS.z0");

    GeomHandle<Beamline> beamg;
    double solenoidOffset = beamg->solenoidOffset();

    CLHEP::Hep3Vector detSolCoilPosition(+solenoidOffset, 0., -dsCoilZ0);

    if ( diagLevel > 0) {
      cout << __func__ << " detSolCoilPosition  : " << detSolCoilPosition  << endl;
    }

    std::vector<double> CRSSteelShieldOffsetSTDV;
    _config.getVectorDouble("fluxcrv.HallSteelOffset", CRSSteelShieldOffsetSTDV, 3);

    CLHEP::Hep3Vector CRSSteelShieldOffset(CRSSteelShieldOffsetSTDV[0],
                                      CRSSteelShieldOffsetSTDV[1],
                                      CRSSteelShieldOffsetSTDV[2]);

    if ( diagLevel > 0) {
      cout << __func__ << " CRSSteelShieldOffset     : " << CRSSteelShieldOffset  << endl;
    }

    // reconstructing the position of the hall air in the world etc...

    vector<double> worldHLen;
    _config.getVectorDouble("world.halfLengths", worldHLen, 3);

    // Get parameters related to the overall dimensions of the hall and to
    // the earthen overburden.

    double floorThick           =  _config.getDouble("hall.floorThick");
    double ceilingThick         =  _config.getDouble("hall.ceilingThick");
    double overburdenDepth      =  _config.getDouble("dirt.overburdenDepth");
    vector<double> hallInHLen;
    _config.getVectorDouble("hall.insideHalfLengths",hallInHLen,3);

    // Top of the floor in G4 world coordinates.
    double yFloor   = -worldHLen[1] + floorThick;

    // Bottom of the ceiling in G4 world coordinates.
    double yCeilingInSide = yFloor + 2.*hallInHLen[1];
    
    // Top of the ceiling in G4 world coordinates.
    double yCeilingOutside  = yCeilingInSide + ceilingThick;

    // Surface of the earth in G4 world coordinates.
    double ySurface  = yCeilingOutside + overburdenDepth;
    
    // Half length and y origin of the dirt box.
    double yLDirt = ( ySurface + worldHLen[1] )/2.;
    double y0Dirt = -worldHLen[1] + yLDirt;
    
    // Center of the dirt box, in the G4 world system.
    CLHEP::Hep3Vector dirtOffset(0.,y0Dirt,0.);

    if ( diagLevel > 0) {
      cout << __func__ << " dirtOffset          : " << dirtOffset    << endl;
    }

    // Position of the center of the hall in the world volume.
    vector<double> hallPosition;
    _config.getVectorDouble("hall.offset", hallPosition,3);

    if ( diagLevel > 0) {
      cout << __func__ << " hallPosition/Offset : " << CLHEP::Hep3Vector(hallPosition[0],
                                                                         hallPosition[1],
                                                                         hallPosition[2])
           << endl;
    }

    double hallY0 = yFloor + hallInHLen[1] + hallPosition[1];

    // Center of the concrete volume in the coordinate system of the dirt.
    // It is a "local" offset

    CLHEP::Hep3Vector wallOffset = 
      CLHEP::Hep3Vector(hallPosition[0], hallY0, hallPosition[2]) - dirtOffset; //

    if ( diagLevel > 0) {
      cout << __func__ << " wallOffset          : " << wallOffset    << endl;
    }

    if ( diagLevel > 0) {
      cout << __func__ << " wallOffsetInWorld   : " << wallOffset + dirtOffset  << endl;
    }

    // Origin of the hall air volume in the system of the hall concrete volume.
    CLHEP::Hep3Vector hallAirOffset = CLHEP::Hep3Vector( 0., (floorThick-ceilingThick)/2., 0.);

    if ( diagLevel > 0) {
      cout << __func__ << " hallAirOffset          : " << hallAirOffset << endl;
    }

    if ( diagLevel > 0) {
      cout << __func__ << " hallAirOffsetInWorld   : " << hallAirOffset + wallOffset + dirtOffset << endl;
    }

    // Position of the origin of the mu2e coordinate system
    CLHEP::Hep3Vector mu2eOrigin = CLHEP::Hep3Vector(
                                                     _config.getDouble("world.mu2eOrigin.xoffset"),
                                                     _config.getDouble("world.mu2eOrigin.height") + yFloor,
                                                     _config.getDouble("world.mu2eOrigin.zoffset")
                                                     );
    if ( diagLevel > 0) {
      cout << __func__ << " mu2eOrigin          : " << mu2eOrigin    << endl;
    }

    // Imagine a box that exactly contains the flux return steel.
    // This is the center of that box, in the coordinate system of the mother volume (the hall air).

    CLHEP::Hep3Vector CRSOffset = - (hallAirOffset + wallOffset + dirtOffset - mu2eOrigin + detSolCoilPosition + CRSSteelShieldOffset);

    if ( diagLevel > 0) {
      cout << __func__ << " CRSOffset           : " << CRSOffset    << endl;
    }

    CLHEP::Hep3Vector CRSOffsetInWorld(CRSOffset + hallAirOffset + wallOffset + dirtOffset);

    if ( diagLevel > 0) {
      cout << __func__ << " CRSOffsetInWorld    : " << CRSOffsetInWorld << endl;
    }

    CLHEP::Hep3Vector TopShieldOffset   (0.,      CRSSteelShieldSideHalfY + CRSSteelShieldHalfThick,  0.);
    CLHEP::Hep3Vector BottomShieldOffset(0.,    -(CRSSteelShieldSideHalfY + CRSSteelShieldHalfThick), 0.);
    CLHEP::Hep3Vector LeftShieldOffset  (         CRSSteelShieldSideHalfY + CRSSteelShieldHalfThick,  0., 0.);
    CLHEP::Hep3Vector RightShieldOffset (       -(CRSSteelShieldSideHalfY + CRSSteelShieldHalfThick), 0., 0.);
    CLHEP::Hep3Vector BackShieldOffset  (0., 0.,  CRSSteelShieldSideHalfZ - CRSSteelShieldHalfThick);
    CLHEP::Hep3Vector FrontShieldOffset (0., 0.,-(CRSSteelShieldSideHalfZ - CRSSteelShieldHalfThick));

    if ( diagLevel > 0) {

      cout << __func__ << " TopShieldOffset     : " << TopShieldOffset    << endl;
      cout << __func__ << " BottomShieldOffset  : " << BottomShieldOffset << endl;
      cout << __func__ << " LeftShieldOffset    : " << LeftShieldOffset   << endl;
      cout << __func__ << " RightShieldOffset   : " << RightShieldOffset  << endl;
      cout << __func__ << " BackShieldOffset    : " << BackShieldOffset   << endl;
      cout << __func__ << " FrontShieldOffset   : " << FrontShieldOffset  << endl;

    }

    // finaly create the steel shield objects

     _crs->addSteelShield("CRSSteelTopShield",
                          TopShieldOffset + CRSOffset, // this is local  offset in the hall Air
                          0,                           // HepRotation*
                          TopShieldOffset + CRSOffsetInWorld, // this is global offset in World
                          CRSSteelShieldTopDims);

     _crs->addSteelShield("CRSSteelBottomShield",
                          BottomShieldOffset + CRSOffset,
                          0,
                          BottomShieldOffset + CRSOffsetInWorld,
                          CRSSteelShieldTopDims);

     _crs->addSteelShield("CRSSteelLeftShield",
                          LeftShieldOffset + CRSOffset,
                          0,
                          LeftShieldOffset + CRSOffsetInWorld,
                          CRSSteelShieldSideDims);

     _crs->addSteelShield("CRSSteelRightShield",
                          RightShieldOffset + CRSOffset,
                          0,
                          RightShieldOffset + CRSOffsetInWorld,
                          CRSSteelShieldSideDims);

     _crs->addSteelShield("CRSSteelBackShield",
                          BackShieldOffset + CRSOffset,
                          0,
                          BackShieldOffset + CRSOffsetInWorld,
                          CRSSteelShieldFrontDims);

     _crs->addSteelShield("CRSSteelFrontShield",
                          FrontShieldOffset + CRSOffset,
                          0,
                          FrontShieldOffset + CRSOffsetInWorld,
                          CRSSteelShieldFrontDims,
                          _config.getDouble("fluxcrv.HallSteelHoleRadius")
                          );

  }
  
  CosmicRayShieldMaker::~CosmicRayShieldMaker (){}

} // namespace mu2e

