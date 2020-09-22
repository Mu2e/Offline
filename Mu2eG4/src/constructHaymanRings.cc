//
// Free function to create  Production Solenoid and Production Target.
//
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
#include "Mu2eG4/inc/constructHaymanRings.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/finishNesting.hh"

#include "ProductionSolenoidGeom/inc/PSVacuum.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Polycone.hh"
#include "G4SDManager.hh"
#include "G4Trd.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
using namespace std;

namespace mu2e {

  void constructHaymanRings(VolumeInfo const & parent, SimpleConfig const & _config) {

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

    std::cout << "production target mother info \n" <<
      tgt->position() << " " << parent.centerInMu2e() << std::endl;

    double const _supWheel_trgtPS_rIn         = _config.getDouble("supWheel_trgtPS_rIn");
    double const _supWheel_trgtPS_rOut        = _config.getDouble("supWheel_trgtPS_rOut");
    G4Material* _suppWheelMaterial = findMaterialOrThrow(_config.getString("supWheel_trgtPS_material"));

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "ProductionTargetSupportWheel", "supWheel_trgtPS" );

    G4ThreeVector _loclCenter(0.0,0.0,0.0);

    //    TubsParams suppWheelParams( _supWheel_trgtPS_rIn, _supWheel_trgtPS_rOut, _supWheel_trgtPS_halfLength);
    // making infinitely small support wheel
    TubsParams suppWheelParams( _supWheel_trgtPS_rIn, _supWheel_trgtPS_rOut, .0000001);
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

    std::cout << __func__ << " radius and halfLength =" << tgt->rOut() << " " << tgt->halfLength() << std::endl;
    std::cout << __func__ << " rotation angle        =" << tgt->productionTargetRotation() << std::endl;

    G4Material* prodTargetMaterial = findMaterialOrThrow(_config.getString("targetPS_materialName"));
    geomOptions->loadEntry( _config, "ProductionTarget", "targetPS" );

    //
    // this extra Core word because in making a sensitive detector, the code sees "ProductionTarget" in
    // "ProductionTarget*" and screws up.  Need something unique.  easier than fixing LV.name function 
    VolumeInfo prodTargetInfo   = nestTubs( "ProductionTargetCore",
                                            prodTargetParams,
                                            prodTargetMaterial,
                                            &tgt->productionTargetRotation(),
                                            _loclCenter,
                                            prodTargetMotherInfo,
                                            0,
                                            G4Colour::Yellow()
                                            );

    //
    // support rings, from Pushka email on 15 April 2019: from 15 to 17 mm, 2 mm half-length
    double innerRadiusRing = 15.;
    double outerRadiusRing = 17.;
    double halfLengthRing = 2.;
    TubsParams supportRingsParams(innerRadiusRing,outerRadiusRing,halfLengthRing);
    //    TubsParams supportRingsParams( 0., tgt->rOut(), tgt->halfLength());

    //
    //now move rings into position.  Since there is just a y-rotation it's easy.  Remember that the central z of the ring changes too; rotate about the y-axis,
    // about a point at the center, dx = rsin,dz = rcos

    /*
    //double deltaX = sin(tgt->productionTargetRotation().getTheta())*(tgt->halfLength()+halfLengthRing);
    //double deltaZ = cos(tgt->productionTargetRotation().getTheta())*(tgt->halfLength()+halfLengthRing);
    //    CLHEP::Hep3Vector SupportRingFarTSOffset (+deltaX, 0., +deltaZ);
    //CLHEP::Hep3Vector SupportRingNearTSOffset(-deltaX, 0., -deltaZ);
    CLHEP::Hep3Vector SupportRingFarTSOffset (+0., 0., +0.);
    CLHEP::Hep3Vector SupportRingNearTSOffset(+0., 0., +0.);
    VolumeInfo supportRingsFarTSInfo = nestTubs("SupportRingsFarTS",
						supportRingsParams,
						supportRingsMaterial,
						&tgt->productionTargetRotation(),
						_loclCenter+SupportRingFarTSOffset,
						prodTargetMotherInfo,
						0,
						G4Colour::Yellow()
						);

    VolumeInfo supportRingsNearTSInfo = nestTubs("SupportRingsNearTS",
						 supportRingsParams,
						 supportRingsMaterial,
						 &tgt->productionTargetRotation(),
						 _loclCenter+SupportRingNearTSOffset,
						 prodTargetMotherInfo,
						 0,
						 G4Colour::Yellow()
						 );
    std::cout << "local center = " << _loclCenter << std::endl;
   std::cout << "local center far = " << SupportRingFarTSOffset << std::endl;
   std::cout << "local center near = " << SupportRingNearTSOffset << std::endl;
    */
    // Now add fins for version 2
    // this is hayman target; first version has four fins, and looks very much like strawman 3.4 of doc-db 27281
    if ( tgt->version() > 1 ) {

      double finHalfLength = tgt->halfLength();
      //
      // in hayman, there are no hubs and in space no one can hear you scream.  Leave the code in for hooks later
      
      G4Trd * myTrd = new G4Trd("FinTrapezoid",
				//				finHalfLength, finHalfLengthOut,
				tgt->finThickness()/2.0, tgt->finThickness()/2.0,
				tgt->finHeight()/2.0,tgt->finHeight()/2.0,
				finHalfLength);


      G4Tubs* targetEndRingNearTS = new G4Tubs("targetEndRingFarTS",
					    innerRadiusRing,outerRadiusRing,halfLengthRing,0.,2.*M_PI);
      G4Tubs* targetEndRingFarTS = new G4Tubs("targetEndRingNearTS",
					    innerRadiusRing,outerRadiusRing,halfLengthRing,0.,2.*M_PI);

      VolumeInfo targetEndRingNearTSVol("targetEndRingNearTS",
					   _loclCenter,prodTargetMotherInfo.centerInWorld);
      VolumeInfo targetEndRingFarTSVol("targetEndRingFarTS", 
					   _loclCenter,prodTargetMotherInfo.centerInWorld);
      G4Material* endRingsMaterial = findMaterialOrThrow("G4_W");


     CLHEP::HepRotation* rotRingBase = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);

      CLHEP::HepRotation* rotRing = new CLHEP::HepRotation((*rotRingBase)*tgt->productionTargetRotation());

      // std::vector<double> finDims = {tgt->finThickness()/2.0,tgt->finHeight()/2.0,finHalfLength};
      double rToFin = tgt->rOut()+tgt->finHeight()/2.0;
      std::cout << "r variables " << rToFin << " " << tgt->rOut() << " " << tgt->finHeight()/2. << std::endl;
  
 
      CLHEP::HepRotation* rotFinBase = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
 
      std::cout << "rotfinbase = " << *rotFinBase << std::endl;
      CLHEP::HepRotation* rotFin1 = new CLHEP::HepRotation((*rotFinBase)*tgt->productionTargetRotation());
      CLHEP::HepRotation* rotFin2 = new CLHEP::HepRotation((*rotFinBase)*tgt->productionTargetRotation());
      CLHEP::HepRotation* rotFin3 = new CLHEP::HepRotation((*rotFinBase)*tgt->productionTargetRotation());
      CLHEP::HepRotation* rotFin4 = new CLHEP::HepRotation((*rotFinBase)*tgt->productionTargetRotation());


      rotFin1->rotateZ(-M_PI/4.);
      rotFin2->rotateZ(+M_PI/4.);
      rotFin3->rotateZ(+M_PI/4.);
      rotFin4->rotateZ(-M_PI/4.);

      CLHEP::Hep3Vector finOffset1(rToFin/sqrt(2.),-rToFin/sqrt(2.),0.);
      CLHEP::Hep3Vector finOffset2(-rToFin/sqrt(2.),-rToFin/sqrt(2.),0.);
      CLHEP::Hep3Vector finOffset3(rToFin/sqrt(2.),rToFin/sqrt(2.),0.);
      CLHEP::Hep3Vector finOffset4(-rToFin/sqrt(2.),rToFin/sqrt(2.),0.);
 
      // These are shifts in the unrotated frame in x and y.  But then when I apply
      // the rotation, since the fins are not centered on the z-axis, they pick up a z-shift that I must
      // take out.  rather than calculate it for the special case I'll do the matrix. If life gets more complicated
      // in the future (target rotated along and y) it's straightforward if tedious.
      //
      // it's complicated because the "finishNesting", rotates then shifts by what you give it; the above is wrong once the frame is 
      // rotated... 
      CLHEP::Hep3Vector fin1Shift(rToFin/sqrt(2.),-rToFin/sqrt(2.),0.);
      CLHEP::Hep3Vector fin2Shift(-rToFin/sqrt(2.),-rToFin/sqrt(2.),0.);
      CLHEP::Hep3Vector fin3Shift(rToFin/sqrt(2.),rToFin/sqrt(2.),0.);
      CLHEP::Hep3Vector fin4Shift(-rToFin/sqrt(2.),rToFin/sqrt(2.),0.);
      //      std::cout << "fin1shift before = " << fin1Shift << std::endl;
      //std::cout << "fin2shift before = " << fin2Shift << std::endl;
      //std::cout << "fin3shift before = " << fin3Shift << std::endl;
      //std::cout << "fin4shift before = " << fin4Shift << std::endl;

      fin1Shift.transform(*rotFin1);
      fin2Shift.transform(*rotFin2);
      fin3Shift.transform(*rotFin3);
      fin4Shift.transform(*rotFin4);

      //std::cout << "fin1shift after = " << fin1Shift << std::endl;
      //std::cout << "fin2shift after = " << fin2Shift << std::endl;
      //std::cout << "fin3shift after = " << fin3Shift << std::endl;
      //std::cout << "fin4shift after = " << fin4Shift << std::endl;

      finOffset1 += CLHEP::Hep3Vector(0.,0.,-fin1Shift.z());
      finOffset2 += CLHEP::Hep3Vector(0.,0.,-fin2Shift.z());
      finOffset3 += CLHEP::Hep3Vector(0.,0.,-fin3Shift.z());
      finOffset4 += CLHEP::Hep3Vector(0.,0.,-fin4Shift.z());

      VolumeInfo fin1Vol("ProductionTargetFin1",
			 _loclCenter,
			 prodTargetMotherInfo.centerInWorld);

      VolumeInfo fin2Vol("ProductionTargetFin2",
			 _loclCenter,
			 prodTargetMotherInfo.centerInWorld);

      VolumeInfo fin3Vol("ProductionTargetFin3",
			 _loclCenter,
			 prodTargetMotherInfo.centerInWorld);

      VolumeInfo fin4Vol("ProductionTargetFin4",
			 _loclCenter,
			 prodTargetMotherInfo.centerInWorld);

      fin1Vol.solid = myTrd;
      fin2Vol.solid = myTrd;
      fin3Vol.solid = myTrd;
      fin4Vol.solid = myTrd;

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
		    

      finishNesting( fin4Vol,
		     prodTargetMaterial,
		     rotFin4,
		     finOffset4,
		     prodTargetMotherInfo.logical,
		     0,
		     G4Colour::Magenta(),
		     "PS"
		     );
 

      targetEndRingNearTSVol.solid = targetEndRingNearTS;
      targetEndRingFarTSVol.solid = targetEndRingFarTS;

      //
      // put these just past edge of target in z, and shift them in x because target rotated about y axis.
      // rotation is about y, and xprime  = x cos - z sin, zprime = z cos + x sin.  x is zero, z is L/2 + whatever I need to push past core
      //but also recall 14deg from z-axis, not from x-axis...
      double lengthToEnd = tgt->halfLength()+(4.); //4. being half length of ring.  
      double deltaXRing =  sin(tgt->productionTargetRotation().getTheta());
      double deltaZRing =  cos(tgt->productionTargetRotation().getTheta());
      //
      // these offets are actual placements.  I just took absolute values for sin/cos above and deal with that in the offsets
      CLHEP::Hep3Vector farTSRingOffset (-deltaXRing*lengthToEnd,0.,-lengthToEnd*(deltaZRing));
      CLHEP::Hep3Vector nearTSRingOffset(+deltaXRing*lengthToEnd,0., lengthToEnd*(deltaZRing));
      //      std::cout << "angle of rotation in degrees = " << TMath::RadToDeg()*tgt->productionTargetRotation().getTheta() << std::endl;
      //std::cout << "length to end = " << lengthToEnd << std::endl;
      //std::cout << "ring offsets x,z = " << deltaXRing << " \n " << deltaZRing << std::endl;
     finishNesting(targetEndRingFarTSVol,
		    endRingsMaterial,
		    rotRing,
		    farTSRingOffset,
		    prodTargetMotherInfo.logical,
		    0,
		    G4Colour::Magenta(),
		    "PS"
		    );

      finishNesting(targetEndRingNearTSVol,
		    endRingsMaterial,
		    rotRing,
		    nearTSRingOffset,
		    prodTargetMotherInfo.logical,
		    0,
		    G4Colour::Magenta(),
		    "PS"
		    );


		    
    }



    // Using the old terms "right" and "left" to mean "downstream" (DS)
    // and "upstream" (US), respectively.

    // 
    //hayman has no hubs.  In this approximation the spokes will hang in space.  
    
    /*
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
    
    */
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
