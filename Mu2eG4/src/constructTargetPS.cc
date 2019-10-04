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
#include "GeomPrimitives/inc/Box.hh"
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
#include "GeometryService/inc/ProductionTargetMaker.hh"
#include "ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Polycone.hh"
#include "G4SDManager.hh"
#include "G4Trd.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/GenericFunctions/ASin.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"

#include "TMath.h"
using namespace std;

namespace mu2e {



  void constructTargetPS(VolumeInfo const & parent, SimpleConfig const & _config) {

    //
    // which target am I building?
    if (_config.getInt("targetPS_version") == ProductionTargetMaker::tier1){
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
      std::cout << "target position and hall origin = " << tgt->position() << "\n" <<
	parent.centerInMu2e() << std::endl;
      std::cout << "supWheel and envHalfLength = " << _clamp_supWheel_rOut  << "\n" <<
	envHalfLength << std::endl;

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
    }
    else if (_config.getInt("targetPS_version") == ProductionTargetMaker::hayman_v_2_0){
     int verbosityLevel                  = _config.getInt("PSHayman.verbosityLevel");
      verbosityLevel >0 &&
	cout << __func__ << " verbosityLevel on Hayman2.0              : " << verbosityLevel  << endl;
      G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();
      // Create variable to avoid multiple look-up

      G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
      
      bool prodTargetVisible   = geomOptions->isVisible( "ProductionTarget" );
      bool prodTargetSolid     = geomOptions->isSolid  ( "ProductionTarget" );
      bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible( "ProductionTarget" );
      bool placePV             = geomOptions->placePV( "ProductionTarget" );
      bool doSurfaceCheck      = geomOptions->doSurfaceCheck( "ProductionTarget" );
      

      // 

      // begin all names with ProductionTarget so when we build sensitive detectors we can make them all 
      // sensitive at once with LVname.find("ProductionTarget") !=std::string::npos
      //
      // Build the production target.
      GeomHandle<ProductionTarget> tgt;
      G4Material* prodTargetCoreMaterial = findMaterialOrThrow(tgt->targetCoreMaterial());
      G4Material* prodTargetFinMaterial  = findMaterialOrThrow(tgt->targetFinMaterial());
      TubsParams prodTargetMotherParams( 0.
      					 ,tgt->productionTargetMotherOuterRadius()
      					 ,tgt->productionTargetMotherHalfLength());

      G4ThreeVector _loclCenter(0.0,0.0,0.0);
      CLHEP::Hep3Vector zeroTranslation(0.,0.,0.);
      G4RotationMatrix* identityRotation = new G4RotationMatrix(CLHEP::HepRotation::IDENTITY);
      VolumeInfo prodTargetMotherInfo   = nestTubs( "ProductionTargetMother",
						    prodTargetMotherParams,
						    parent.logical->GetMaterial(),
						    0,
						    tgt->haymanProdTargetPosition() - _hallOriginInMu2e,
						    parent,
						    0,
						    prodTargetVisible,
						    G4Colour::Gray(),
						    prodTargetSolid,
						    forceAuxEdgeVisible,
						    placePV,
						    doSurfaceCheck
						    );

      std::cout << "target position and hall origin = " << tgt->haymanProdTargetPosition() << "\n" <<
	_hallOriginInMu2e << std::endl;
      if (verbosityLevel > 0){
	std::cout << __func__ << "created prodTargetMotherInfo " 
		  << tgt->productionTargetMotherOuterRadius() 
		  << " " <<tgt->productionTargetMotherHalfLength() << std::endl;
      }
      //
      // construct each section

      int numberOfSections = tgt->numberOfTargetSections();
      //
      // first build the core

      //
      // start at most negative end (technically could start at center but this is simpler since target is not symmetric about its center
      // this variable gets incremented every time we add anything
      double _currentZ = _loclCenter.z() - tgt->halfHaymanLength();
      if (verbosityLevel > 0){std::cout << __func__ << " current Z starts at " << _currentZ << std::endl;}
      //
      // running z locations for placing elements
      CLHEP::Hep3Vector _segmentCenter(0.,0.,0.);
      CLHEP::Hep3Vector _startingSegmentCenter(0.,0.,0.);
      double targetRadius = tgt->rOut();
      //
      // keeping track of these for sensitive volume creation
      int coreCopyNumber = 0;
      int finCopyNumber = 0;
      //     numberOfSections = 1;

      double angularSize = TMath::ASin((tgt->haymanFinThickness()/2.)/targetRadius);
      double dSphi = M_PI/2. - angularSize;
      double dPphi = 2.*angularSize;

      for (int ithSection = 0; ithSection < numberOfSections; ++ithSection){
	int numberOfSegments = tgt->numberOfSegmentsPerSection().at(ithSection);
	//	numberOfSegments = 1;
	std::string startName = "ProductionTargetStartingCoreSection_" + std::to_string(ithSection);
	//
	// first place starting segment

	TubsParams startingSegmentParams(0.,targetRadius,tgt->startingSectionThickness().at(ithSection)/2.);
	double currentHalfStartingSegment = tgt->startingSectionThickness().at(ithSection)/2.;
	//	
	//set currentZ to be in the center of the starting section
	_currentZ += tgt->startingSectionThickness().at(ithSection)/2.;
	if (verbosityLevel > 0){
	  std::cout << __func__ << " ithSection , startingSegmentCenter at " << ithSection << " " <<_currentZ << std::endl;
	  std::cout << "    and copy number for starting segment is now " << coreCopyNumber << std::endl;
	}
 	_startingSegmentCenter.setZ(_currentZ);

	CLHEP::Hep3Vector currentStartingSegmentCenter = tgt->productionTargetRotation().inverse()*_startingSegmentCenter;
	VolumeInfo startingSegment = nestTubs(startName,
					      startingSegmentParams,
					      prodTargetCoreMaterial,
					      &tgt->productionTargetRotation(),
					      currentStartingSegmentCenter,
					      prodTargetMotherInfo,
					      ithSection,
					      prodTargetVisible,
					      G4Colour::Yellow(),
					      prodTargetSolid,
					      forceAuxEdgeVisible,
					      placePV,
					      doSurfaceCheck
					      );

	//build fins around starting segment
	for (int ithFin = 0; ithFin < tgt->nHaymanFins(); ++ithFin)
	  //for (int ithFin = 0; ithFin < 1; ++ithFin)
	  {
	    //
	    // what we do is make a rectangular fin; cut out the appropriate core section with a subtraction volume;
	    // then that object is placed.
	    double currentFinAngle = tgt->finAngles().at(ithFin)*CLHEP::degree;
	    const std::string name = "ProductionTargetStartingSegmentSection_" + std::to_string(ithSection) + std::to_string(ithFin);
	    // divide by 2 to make box half-height since I gave it the "radius" to start with; arbitrary convention. 
	    double finHalfHeightAboveTarget = (tgt->finOuterRadius() - tgt->rOut())/2.;
	    G4Box* finBox = new G4Box("finBox",tgt->haymanFinThickness()/2.,finHalfHeightAboveTarget,currentHalfStartingSegment);
	    G4Tubs* tubCutout = new G4Tubs("finCutout",0.,tgt->rOut(),currentHalfStartingSegment+.0001,dSphi,dPphi);// need extra length for visualization tool
	    CLHEP::HepRotation* rotFin = new CLHEP::HepRotation(tgt->productionTargetRotation().inverse());
	    CLHEP::HepRotation* rotFinG4 = new CLHEP::HepRotation(tgt->productionTargetRotation());
	    rotFin->rotateZ(currentFinAngle);
	    rotFinG4->rotateZ(-currentFinAngle);
	    ++finCopyNumber;

	    G4ThreeVector finTranslation = currentStartingSegmentCenter + CLHEP::Hep3Vector(0.,finHalfHeightAboveTarget + tgt->rOut(),0.);
	    finTranslation = (*rotFin)*finTranslation;


	    CLHEP::Hep3Vector downshift = CLHEP::Hep3Vector(0.,-finHalfHeightAboveTarget - targetRadius*TMath::Cos(angularSize), 0.);
	    G4SubtractionSolid* finWithCutoutSolid   = new G4SubtractionSolid(name,finBox,tubCutout,identityRotation,downshift);
	    G4LogicalVolume*    finWithCutoutLogical = new G4LogicalVolume(finWithCutoutSolid,prodTargetFinMaterial,name);
	    G4VPhysicalVolume*  finWithCutout        = new G4PVPlacement(rotFinG4,finTranslation,finWithCutoutLogical,name,prodTargetMotherInfo.logical,0,finCopyNumber,false);

	    //
	    // this check also uses finWithCutout.  Indirectly it's a way of forcing an overlap check or it won't compile, since finWithCutout is never used otherwise.
	    if (doSurfaceCheck){checkForOverlaps(finWithCutout,_config,0);}
	  }

	//
	// set current z to be at the end of this starting segment before beginning loop on sections
	_currentZ += tgt->startingSectionThickness().at(ithSection)/2.;
	if (verbosityLevel > 0){std::cout << __func__ << " ithSection , startingSegment Z ends at " << ithSection << " " <<_currentZ << std::endl;}
	  double currentGap = tgt->thicknessOfGapPerSection().at(ithSection);
	  double currentHalfSegment = tgt->thicknessOfSegmentPerSection().at(ithSection)/2.;

	for (int ithSegment = 0; ithSegment < numberOfSegments; ++ithSegment){
	  std::string name = "ProductionTargetCoreSection_" + std::to_string(ithSection) 
	    + "_Segment_" + std::to_string(ithSegment);
	  if (verbosityLevel > 0) {
	    std::cout << __func__ << "name = " << name << " and core copy number = " << coreCopyNumber << std::endl;
	    std::cout << __func__ << "gap, and segment half = " << currentGap << " " << currentHalfSegment << std::endl;
	  }
	  TubsParams segmentParams(0.,targetRadius,currentHalfSegment);
	  _currentZ += currentHalfSegment + currentGap;
	  if (verbosityLevel > 0){
	    std::cout << __func__ << " ithSection , current Z for segment center is at " << ithSection << " " 
		      << ithSegment << " "  <<_currentZ << std::endl;
	  }
	  std::cout << "segment center before anything = " << _segmentCenter << std::endl;
	  _segmentCenter.setZ(_currentZ);
	  CLHEP::Hep3Vector currentSegmentCenter = tgt->productionTargetRotation().inverse()*_segmentCenter;
	  VolumeInfo coreSegment = nestTubs(name,
					    segmentParams,
					    prodTargetCoreMaterial,
					    &tgt->productionTargetRotation(),
					    currentSegmentCenter,
					    prodTargetMotherInfo,
					    coreCopyNumber,
					    prodTargetVisible,
					    G4Colour::Yellow(),
					    prodTargetSolid,
					    forceAuxEdgeVisible,
					    placePV,
					    doSurfaceCheck
					    );
	  //
	  //now add fins surrounding this core segment.  Build them as simple rectangles with a G4 subtraction volume for the core.
	  
	  	  for (int ithFin = 0; ithFin < tgt->nHaymanFins(); ++ithFin)
	  //	  for (int ithFin = 0; ithFin < 1; ++ithFin)
	    {
	      //
	      // what we do is make a rectangular fin; cut out the appropriate core section with a subtraction volume;
	      // then that object is placed.
	      double currentFinAngle = tgt->finAngles().at(ithFin)*CLHEP::degree;
	      const std::string name = "ProductionTargetFinSection_" + std::to_string(ithSection) + std::to_string(ithSegment) + std::to_string(ithFin);
	      if (verbosityLevel > 1)
		{std::cout << "Fin Section and Angle = " << name << " " << currentFinAngle << " fin copy number = " << finCopyNumber << std::endl;}
	      // divide by 2 to make box half-height since I gave it the "radius" to start with; arbitrary convention. 
	      double finHalfHeightAboveTarget = (tgt->finOuterRadius() - tgt->rOut())/2.;
	      std::cout << "finHeightAboveTarget = " << tgt->finOuterRadius() << " " << tgt->rOut() << " " << finHalfHeightAboveTarget << std::endl;
	      G4Box* finBox = new G4Box("finBox",tgt->haymanFinThickness()/2.,finHalfHeightAboveTarget,currentHalfSegment);
	      // 
	      // this tubs is in the xy plane with an extent along z.  So phi goes in the xy plane.
	      G4Tubs* tubCutout = new G4Tubs("finCutout",0.,tgt->rOut(),currentHalfSegment+.0001,dSphi,dPphi);// need extra length for visualization tool
	      //
	      // I now have two G4Solids, a box and a cutout.  Combine them to make a logical volume. the box is 
	      // oriented so that it is "z" thick and "radius" high.  the tubs is made to be the same, so no rotation
	      // and also no translation needed 
	      //
	      // and this is confusing but needed.  I move the fin where it needs to go using CLHEP, which requires the inverse matrix 
	      // -- but when I feed it to G4 I need the regular form.  Just to add, the .inverse is active and was inverted already. sigh...
	      // I'm sure there are better ways to do this but this has the advantage of my understanding it
              CLHEP::HepRotation* rotFin = new CLHEP::HepRotation(tgt->productionTargetRotation().inverse());
              CLHEP::HepRotation* rotFinG4 = new CLHEP::HepRotation(tgt->productionTargetRotation());
	      rotFin->rotateZ(currentFinAngle);
	      rotFinG4->rotateZ(-currentFinAngle);
	      ++finCopyNumber;

	      G4ThreeVector finTranslation = currentSegmentCenter + CLHEP::Hep3Vector(0.,finHalfHeightAboveTarget + tgt->rOut(),0.);
	      finTranslation = (*rotFin)*finTranslation;


	      if (verbosityLevel > 1){
		std::cout << "ithfin; current fin angle; rotFin; finTranslation = " << ithFin << " " << currentFinAngle << *rotFin << " " << *rotFinG4 << " " << finTranslation << std::endl;
	      }
	      //
	      //couple of good debugging tricks: make the subtraction volume a union volume, then as you move it around you can see it.  Offset in whatever
	      //direction you need to line up, and still see everything.
	      //moving cutout to "bottom" of fin, then up by target radius so that center of cutout arc is at center of target
	      CLHEP::Hep3Vector downshift = CLHEP::Hep3Vector(0.,-finHalfHeightAboveTarget - targetRadius*TMath::Cos(angularSize), 0.);
	      G4SubtractionSolid* finWithCutoutSolid   = new G4SubtractionSolid(name,finBox,tubCutout,identityRotation,downshift);
	      G4LogicalVolume*    finWithCutoutLogical = new G4LogicalVolume(finWithCutoutSolid,prodTargetFinMaterial,name);
	      G4VPhysicalVolume*  finWithCutout        = new G4PVPlacement(rotFinG4,finTranslation,finWithCutoutLogical,name,prodTargetMotherInfo.logical,0,finCopyNumber,false);

	      //
	      // this check also uses finWithCutout.  Indirectly it's a way of forcing an overlap check or it won't compile, since finWithCutout is never used otherwise.
	      if (doSurfaceCheck){checkForOverlaps(finWithCutout,_config,0);}
	    }
	  _currentZ += currentHalfSegment;
	  //
	  // if this is the last segment, add another gap before next section starts!
	  if (ithSegment == numberOfSegments-1) {
	    if (verbosityLevel > 0){std::cout << " z before = " << _currentZ << std::endl;}
	    _currentZ += currentGap;
	    if (verbosityLevel > 0){std::cout << " z after = " << _currentZ << std::endl;}
	  }
	  if (verbosityLevel > 0){std::cout << __func__ << " ending at z=" << _currentZ << std::endl;}
	  //	  prodTargetCoreInfoBySegment.push_back(coreSegment);
	}
	  ++coreCopyNumber;
	//
	// if this is the last section, add on one more beginning block. decided not to extend vectors and have zero segments, just ugly
	if (ithSection == (numberOfSections -1) ){
	  startName = "ProductionTargetStartingCoreSection_" + std::to_string(ithSection+1);
	  TubsParams startingSegmentParamsEnd(0.,targetRadius,tgt->startingSectionThickness().at(ithSection)/2.);

	  //	
	  //set currentZ to be in the center of the starting section
	  _currentZ += tgt->startingSectionThickness().at(ithSection)/2.;
	  if (verbosityLevel > 0){std::cout << __func__ << " ithSection+1 , startingSegmentCenter at " << ithSection+1 << " " <<_currentZ << std::endl;}
	  _startingSegmentCenter.setZ(_currentZ);
	  CLHEP::Hep3Vector currentStartingSegmentCenter = tgt->productionTargetRotation().inverse()*_startingSegmentCenter;
	  startingSegment = nestTubs(startName,
				     startingSegmentParamsEnd,
				     prodTargetCoreMaterial,
				     &tgt->productionTargetRotation(),
				     currentStartingSegmentCenter,
				     prodTargetMotherInfo,
				     ithSection+1,
				     prodTargetVisible,
				     G4Colour::Yellow(),
				     prodTargetSolid,
				     forceAuxEdgeVisible,
				     placePV,
				     doSurfaceCheck
				     );
	//build fins around starting segment
	for (int ithFin = 0; ithFin < tgt->nHaymanFins(); ++ithFin)
	  //for (int ithFin = 0; ithFin < 1; ++ithFin)
	  {
	    //
	    // what we do is make a rectangular fin; cut out the appropriate core section with a subtraction volume;
	    // then that object is placed.
	    double currentFinAngle = tgt->finAngles().at(ithFin)*CLHEP::degree;
	    const std::string name = "ProductionTargetStartingSegmentSection_" + std::to_string(ithSection) + std::to_string(ithFin);
	    // divide by 2 to make box half-height since I gave it the "radius" to start with; arbitrary convention. 
	    double finHalfHeightAboveTarget = (tgt->finOuterRadius() - tgt->rOut())/2.;
	    G4Box* finBox = new G4Box("finBox",tgt->haymanFinThickness()/2.,finHalfHeightAboveTarget,currentHalfStartingSegment);
	    G4Tubs* tubCutout = new G4Tubs("finCutout",0.,tgt->rOut(),currentHalfStartingSegment+.0001,dSphi,dPphi);// need extra length for visualization tool
	    CLHEP::HepRotation* rotFin = new CLHEP::HepRotation(tgt->productionTargetRotation().inverse());
	    CLHEP::HepRotation* rotFinG4 = new CLHEP::HepRotation(tgt->productionTargetRotation());
	    rotFin->rotateZ(currentFinAngle);
	    rotFinG4->rotateZ(-currentFinAngle);
	    ++finCopyNumber;

	    G4ThreeVector finTranslation = currentStartingSegmentCenter + CLHEP::Hep3Vector(0.,finHalfHeightAboveTarget + tgt->rOut(),0.);
	    finTranslation = (*rotFin)*finTranslation;


	    CLHEP::Hep3Vector downshift = CLHEP::Hep3Vector(0.,-finHalfHeightAboveTarget - targetRadius*TMath::Cos(angularSize), 0.);
	    G4SubtractionSolid* finWithCutoutSolid   = new G4SubtractionSolid(name,finBox,tubCutout,identityRotation,downshift);
	    G4LogicalVolume*    finWithCutoutLogical = new G4LogicalVolume(finWithCutoutSolid,prodTargetFinMaterial,name);
	    G4VPhysicalVolume*  finWithCutout        = new G4PVPlacement(rotFinG4,finTranslation,finWithCutoutLogical,name,prodTargetMotherInfo.logical,0,finCopyNumber,false);

	    //
	    // this check also uses finWithCutout.  Indirectly it's a way of forcing an overlap check or it won't compile, since finWithCutout is never used otherwise.
	    if (doSurfaceCheck){checkForOverlaps(finWithCutout,_config,0);}
	  }

	  //
	  // set current z to be at the end of this starting segment as a final check
	  _currentZ += tgt->startingSectionThickness().at(ithSection)/2.;
	  if (verbosityLevel > 0){
	    std::cout << __func__ << " ithSection+1 , startingSegment Z ends at " << ithSection+1 << " " <<_currentZ << std::endl;
	    std::cout << "             with copy number for last segment  = " << ithSection+1 << std::endl;
	  }
	}
      }
    }
  }
}

 // end Mu2eWorld::constructPS

