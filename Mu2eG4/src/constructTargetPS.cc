//
// Free function to create  Production Solenoid and Production Target.
//
// $Id: constructTargetPS.cc,v 1.2 2014/09/19 19:15:13 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/19 19:15:13 $
//
// Original author KLG based on Mu2eWorld constructPS
//
// Notes:
// Construct the PS. Parent volume is the air inside of the hall.

// C++ includes
#include <iostream>
#include <sstream>
#include <cmath>

// Mu2e includes.
#include "BeamlineGeom/inc/Beamline.hh"
#include "ProductionSolenoidGeom/inc/ProductionSolenoid.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/Polycone.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/constructTargetPS.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/finishNesting.hh"
//#include "Mu2eG4/inc/SensitiveDetectorName.hh"

#include "ProductionSolenoidGeom/inc/PSVacuum.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Polycone.hh"
#include "G4SDManager.hh"
#include "G4Trd.hh"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  void constructTargetPS(VolumeInfo const & parent, SimpleConfig const & _config) {

    int verbosityLevel                  = _config.getInt("PS.verbosityLevel");

    //G4Material* psVacVesselMaterial = findMaterialOrThrow(psVacVesselInnerParams.materialName());

    verbosityLevel >0 &&
      cout << __func__ << " verbosityLevel                   : " << verbosityLevel  << endl;

    //bool psVacuumSensitive = _config.getBool("PS.Vacuum.Sensitive", false);

    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();


    // Build the production target.
    GeomHandle<ProductionTarget> tgt;

    double const _clamp_supWheel_rOut         = _config.getDouble("clamp_supWheel_rOut");
    double const _clamp_supWheel_halfLength   = _config.getDouble("clamp_supWheel_halfLength");
    double const _supWheel_trgtPS_halfLength  = _config.getDouble("supWheel_trgtPS_halfLength");

    double envHalfLength = 2.0*_clamp_supWheel_halfLength+_supWheel_trgtPS_halfLength;
    if (tgt->envHalfLength()>envHalfLength) { envHalfLength = tgt->envHalfLength(); }

    TubsParams prodTargetMotherParams( 0., _clamp_supWheel_rOut, envHalfLength);
    VolumeInfo prodTargetMotherInfo   = nestTubs( "ProductionTargetMother",
                                                  prodTargetMotherParams,
                                                  parent.logical->GetMaterial(),
                                                  0,
                                                  tgt->position() - parent.centerInMu2e(),
                                                  parent,
                                                  0,
                                                  G4Colour::Blue(),
                                                  "PS"
                                                  );

    double const _supWheel_trgtPS_rIn         = _config.getDouble("supWheel_trgtPS_rIn");
    double const _supWheel_trgtPS_rOut        = _config.getDouble("supWheel_trgtPS_rOut");
    G4Material* _suppWheelMaterial = findMaterialOrThrow(_config.getString("supWheel_trgtPS_material"));

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "ProductionTargetSupportWheel", "supWheel_trgtPS" );

    G4ThreeVector _loclCenter(0.0,0.0,0.0);

    TubsParams suppWheelParams( _supWheel_trgtPS_rIn, _supWheel_trgtPS_rOut, _supWheel_trgtPS_halfLength);
    VolumeInfo suppWheelInfo   = nestTubs( "ProductionTargetSupportWheel",
                                           suppWheelParams,
                                           _suppWheelMaterial,
                                           0,
                                           _loclCenter,
                                           prodTargetMotherInfo,
                                           0,
                                           G4Colour::Gray()
                                           );

    double const _clamp_supWheel_rIn         = _config.getDouble("clamp_supWheel_rIn");
    G4Material* _clampSpWheelMaterial = findMaterialOrThrow(_config.getString("clamp_supWheel_material"));
    geomOptions->loadEntry( _config, "ClampSupportWheel_R", "clamp_supWheel" );

    G4ThreeVector _clampPosR(0.0,0.0,_supWheel_trgtPS_halfLength+_clamp_supWheel_halfLength);

    TubsParams clampSpWheelParams( _clamp_supWheel_rIn, _clamp_supWheel_rOut, _clamp_supWheel_halfLength);
    VolumeInfo clampSpWheelInfoR   = nestTubs( "ClampSupportWheel_R",
                                               clampSpWheelParams,
                                               _clampSpWheelMaterial,
                                               0,
                                               _clampPosR,
                                               prodTargetMotherInfo,
                                               0,
                                               G4Colour::Gray()
                                               );

    //    G4ThreeVector _clampPosL(0.0,0.0,-(_supWheel_trgtPS_halfLength+_clamp_supWheel_halfLength));
    //
    //    geomOptions->loadEntry( _config, "ClampSupportWheel_L", "clamp_supWheel" );
    //    VolumeInfo clampSpWheelInfoL   = nestTubs( "ClampSupportWheel_L",
    //                                            clampSpWheelParams,
    //                                            _clampSpWheelMaterial,
    //                                            0,
    //                                            _clampPosL,
    //                                            prodTargetMotherInfo,
    //                                            0,
    //                                            G4Colour::Gray()
    //                                           );

    TubsParams prodTargetParams( 0., tgt->rOut(), tgt->halfLength());

    G4Material* prodTargetMaterial = findMaterialOrThrow(_config.getString("targetPS_materialName"));

    geomOptions->loadEntry( _config, "ProductionTarget", "targetPS" );

    VolumeInfo prodTargetInfo   = nestTubs( "ProductionTarget",
                                            prodTargetParams,
                                            prodTargetMaterial,
                                            &tgt->productionTargetRotation(),
                                            _loclCenter,
                                            prodTargetMotherInfo,
                                            0,
                                            G4Colour::Magenta()
                                            );

    // Now add fins for version 2
    if ( tgt->version() > 1 ) {
      // Length of fins must be calculated:

      double finHalfLength = (tgt->halfLength() * 2.0 - tgt->hubDistUS() 
			      - tgt->hubDistDS() - tgt->hubLenUS() - tgt->hubLenDS())/2.0;  // on the inner side, 
                                         // adjacent to the main target body.

      // Use the steeper of the two hub angles to angle the ends of the fin
      double theAngle = tgt->hubAngleUS();
      if ( tgt->hubAngleDS() > theAngle ) theAngle = tgt->hubAngleDS();

      double finHalfLengthOut = finHalfLength + tgt->finHeight()*std::cos(theAngle*CLHEP::degree);  

      G4Trd * myTrd = new G4Trd("FinTrapezoid",
				finHalfLength, finHalfLengthOut,
				tgt->finThickness()/2.0, tgt->finThickness()/2.0,
				tgt->finHeight()/2.0);

      // std::vector<double> finDims = {tgt->finThickness()/2.0,tgt->finHeight()/2.0,finHalfLength};
      double finZoff = (tgt->hubDistUS() - tgt->hubDistDS() + tgt->hubLenUS() - tgt->hubLenDS())/2.0; // z-offset for fin
      double rToFin = tgt->rOut()+tgt->finHeight()/2.0+0.1;

      double xMove = finZoff * sin(tgt->productionTargetRotation().theta());
      CLHEP::Hep3Vector finOffset1(xMove,-rToFin,finZoff);
      CLHEP::Hep3Vector finOffset2(rToFin*cos(M_PI/6.0)+xMove,rToFin*sin(M_PI/6.0),finZoff-1.0);
      CLHEP::Hep3Vector finOffset3(-rToFin*cos(M_PI/6.0)+xMove-0.1,rToFin*sin(M_PI/6.0)+0.02,finZoff+0.2);
      CLHEP::HepRotation* rotFinBase = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
      rotFinBase->rotateX(90.0*CLHEP::degree);
      rotFinBase->rotateZ(90.0*CLHEP::degree);

      CLHEP::HepRotation* rotFin1 = new CLHEP::HepRotation((*rotFinBase)*tgt->productionTargetRotation());
      CLHEP::HepRotation* rotFin2 = new CLHEP::HepRotation((*rotFinBase)*tgt->productionTargetRotation());
      CLHEP::HepRotation* rotFin3 = new CLHEP::HepRotation((*rotFinBase)*tgt->productionTargetRotation());
      rotFin1->rotateX(M_PI);
      rotFin2->rotateX(M_PI/3.0);
      rotFin3->rotateX(5.0*M_PI/3.0);
      VolumeInfo fin1Vol("ProductionTargetFin1",
			 _loclCenter,
			 prodTargetMotherInfo.centerInWorld);

      VolumeInfo fin2Vol("ProductionTargetFin2",
			 _loclCenter,
			 prodTargetMotherInfo.centerInWorld);

      VolumeInfo fin3Vol("ProductionTargetFin3",
			 _loclCenter,
			 prodTargetMotherInfo.centerInWorld);

      fin1Vol.solid = myTrd;
      fin2Vol.solid = myTrd;
      fin3Vol.solid = myTrd;

      finishNesting( fin1Vol,
		     prodTargetMaterial,
		     rotFin1,
		     finOffset1,
		     prodTargetMotherInfo.logical,
		     0,
		     G4Colour::Magenta(),
		     "PS"
		     );

      finishNesting( fin2Vol,
		     prodTargetMaterial,
		     rotFin2,
		     finOffset2,
		     prodTargetMotherInfo.logical,
		     0,
		     G4Colour::Magenta(),
		     "PS"
		     );

      finishNesting( fin3Vol,
		     prodTargetMaterial,
		     rotFin3,
		     finOffset3,
		     prodTargetMotherInfo.logical,
		     0,
		     G4Colour::Magenta(),
		     "PS"
		     );
		    
    }



    // Using the old terms "right" and "left" to mean "downstream" (DS)
    // and "upstream" (US), respectively.

    Polycone const & pHubRgtParams = *tgt->getHubsRgtPtr();
    VolumeInfo prodTargetHubRgtInfo  = nestPolycone("ProductionTargetHubRgt",
                                                    pHubRgtParams.getPolyconsParams(),
                                                    findMaterialOrThrow(pHubRgtParams.materialName()),
                                                    &tgt->productionTargetRotation(),
                                                    pHubRgtParams.originInMu2e(),
                                                    prodTargetMotherInfo,
                                                    0,
                                                    G4Colour::Magenta(),
                                                    "ProductionTarget"
                                                    );

    Polycone const & pHubLftParams = *tgt->getHubsLftPtr();
    VolumeInfo prodTargetHubLftInfo  = nestPolycone("ProductionTargetHubLft",
                                                    pHubLftParams.getPolyconsParams(),
                                                    findMaterialOrThrow(pHubLftParams.materialName()),
                                                    &tgt->productionTargetRotation(),
                                                    pHubLftParams.originInMu2e(),
                                                    prodTargetMotherInfo,
                                                    0,
                                                    G4Colour::Magenta(),
                                                    "ProductionTarget"
                                                    );

    CLHEP::Hep3Vector zax(0,0,1);
    double spokeRad = 0.5*_config.getDouble("targetPS_Spoke_diameter");
    G4Material* spokeMaterial = findMaterialOrThrow(_config.getString("targetPS_Spoke_materialName"));
    double hub_overhang_angle = _config.getDouble("targetPS_Hub_overhang_angle")*CLHEP::deg;
    double normalOverhangRgt = CLHEP::halfpi+hub_overhang_angle;
    double normalOverhangLft = CLHEP::halfpi-hub_overhang_angle;
    CLHEP::HepRotation invTrgtRot = tgt->productionTargetRotation().inverse();

    std::map<double,CLHEP::Hep3Vector> const & anchoringPntsRgt = tgt->anchoringPntsRgt();

    int iSpk(0);

    // Create variable to avoid multiple look-up
    const bool prodTargetVisible   = geomOptions->isVisible( "ProductionTarget" );
    const bool prodTargetSolid     = geomOptions->isSolid  ( "ProductionTarget" );
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible( "ProductionTarget" );
    const bool placePV             = geomOptions->placePV( "ProductionTarget" );
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck( "ProductionTarget" );

    for (std::map<double,CLHEP::Hep3Vector>::const_iterator iPnt=anchoringPntsRgt.begin(); iPnt!=anchoringPntsRgt.end(); ++iPnt) {
      double tmpAngle = iPnt->first * CLHEP::deg;
      CLHEP::Hep3Vector normalax(cos(tmpAngle),sin(tmpAngle),0.0);
      CLHEP::Hep3Vector tmpAnchPnt = _supWheel_trgtPS_rIn*normalax;
      CLHEP::Hep3Vector tmpDirVec = tmpAnchPnt-iPnt->second;
      CLHEP::Hep3Vector tmpMidPnt(0.5*(tmpAnchPnt.x()+iPnt->second.x()),0.5*(tmpAnchPnt.y()+iPnt->second.y()),0.5*(tmpAnchPnt.z()+iPnt->second.z()));
      double tmpSpokeLength = tmpDirVec.mag();
      tmpDirVec = tmpDirVec.unit();
      CLHEP::Hep3Vector tmpRotAxis = tmpDirVec.cross(zax);
      double rotAngle = tmpDirVec.angle(zax);
      double deepAngleWheel = tmpDirVec.angle(normalax);

      CLHEP::Hep3Vector normalaxHub(sin(normalOverhangRgt),0,cos(normalOverhangRgt));
      normalaxHub.rotateZ(tmpAngle);
      normalaxHub.transform(invTrgtRot);
      double deepAngleHub = tmpDirVec.angle(normalaxHub);

      double fixOverlapWheel = spokeRad*tan(deepAngleWheel);
      double fixOverlapHub = spokeRad*tan(deepAngleHub);
      tmpSpokeLength-=(fixOverlapWheel+fixOverlapHub);
      tmpMidPnt += 0.5*(fixOverlapHub-fixOverlapWheel)*tmpDirVec;

      TubsParams iSpokeParams( 0.0, spokeRad, 0.5*tmpSpokeLength);
      std::stringstream iSpokeName;
      iSpokeName<<"ProductionTargetSpokeRgt_"<<iSpk;

      VolumeInfo iSpokeInfo   = nestTubs( iSpokeName.str(),
                                          iSpokeParams,
                                          spokeMaterial,
                                          new CLHEP::HepRotation(tmpRotAxis,rotAngle),
                                          tmpMidPnt,
                                          prodTargetMotherInfo,
                                          0,
                                          prodTargetVisible,
                                          G4Colour::Gray(),
                                          prodTargetSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );
      ++iSpk;
    }

    std::map<double,CLHEP::Hep3Vector> const & anchoringPntsLft = tgt->anchoringPntsLft();

    iSpk=0;
    for (std::map<double,CLHEP::Hep3Vector>::const_iterator iPnt=anchoringPntsLft.begin(); iPnt!=anchoringPntsLft.end(); ++iPnt) {
      double tmpAngle = iPnt->first * CLHEP::deg;
      CLHEP::Hep3Vector normalax(cos(tmpAngle),sin(tmpAngle),0.0);
      CLHEP::Hep3Vector tmpAnchPnt = _supWheel_trgtPS_rIn*normalax;
      CLHEP::Hep3Vector tmpDirVec = tmpAnchPnt-iPnt->second;
      CLHEP::Hep3Vector tmpMidPnt(0.5*(tmpAnchPnt.x()+iPnt->second.x()),0.5*(tmpAnchPnt.y()+iPnt->second.y()),0.5*(tmpAnchPnt.z()+iPnt->second.z()));
      double tmpSpokeLength = tmpDirVec.mag();
      tmpDirVec = tmpDirVec.unit();
      CLHEP::Hep3Vector tmpRotAxis = tmpDirVec.cross(zax);
      double rotAngle = tmpDirVec.angle(zax);
      double deepAngleWheel = tmpDirVec.angle(normalax);

      CLHEP::Hep3Vector normalaxHub(sin(normalOverhangLft),0,cos(normalOverhangLft));
      normalaxHub.rotateZ(tmpAngle);
      normalaxHub.transform(invTrgtRot);
      double deepAngleHub = tmpDirVec.angle(normalaxHub);
      //CLHEP::HepRotation *tmpSpokeRot = new CLHEP::HepRotation(tmpRotAxis,rotAngle);

      double fixOverlapWheel = spokeRad*tan(deepAngleWheel);
      double fixOverlapHub = spokeRad*tan(deepAngleHub);
      tmpSpokeLength-=(fixOverlapWheel+fixOverlapHub);
      tmpMidPnt += 0.5*(fixOverlapHub-fixOverlapWheel)*tmpDirVec;

      TubsParams iSpokeParams( 0.0, spokeRad, 0.5*tmpSpokeLength);
      std::stringstream iSpokeName;
      iSpokeName<<"ProductionTargetSpokeLft_"<<iSpk;

      VolumeInfo iSpokeInfo   = nestTubs( iSpokeName.str(),
                                          iSpokeParams,
                                          spokeMaterial,
                                          new CLHEP::HepRotation(tmpRotAxis,rotAngle),
                                          tmpMidPnt,
                                          prodTargetMotherInfo,
                                          0,
                                          prodTargetVisible,
                                          G4Colour::Gray(),
                                          prodTargetSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );
      ++iSpk;
    }

  } // end Mu2eWorld::constructPS
}
