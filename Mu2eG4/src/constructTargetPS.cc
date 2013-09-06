//
// Free function to create  Production Solenoid and Production Target.
//
// $Id: constructTargetPS.cc,v 1.1 2013/09/06 19:39:18 tassiell Exp $
// $Author: tassiell $
// $Date: 2013/09/06 19:39:18 $
//
// Original author KLG based on Mu2eWorld constructPS
//
// Notes:
// Construct the PS. Parent volume is the air inside of the hall.

// C++ includes
#include <iostream>
#include <sstream>

// Mu2e includes.
#include "BeamlineGeom/inc/Beamline.hh"
#include "ProductionSolenoidGeom/inc/ProductionSolenoid.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/Polycone.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/constructTargetPS.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

#include "ProductionSolenoidGeom/inc/PSVacuum.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Polycone.hh"
#include "G4SDManager.hh"

#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  void constructTargetPS(VolumeInfo const & parent, SimpleConfig const & _config) {
    
    ProductionSolenoid const & psgh = *(GeomHandle<ProductionSolenoid>());

    //Tube const & psVacVesselInnerParams     = *psgh.getVacVesselInnerParamsPtr();
    //Tube const & psVacVesselOuterParams     = *psgh.getVacVesselOuterParamsPtr();

    // Extract some information from the config file.

    int verbosityLevel                  = _config.getInt("PS.verbosityLevel");
    bool const _psVisible               = _config.getBool("PS.visible");
    bool const _psSolid                 = _config.getBool("PS.solid");
                                                
    //G4Material* psVacVesselMaterial = findMaterialOrThrow(psVacVesselInnerParams.materialName());

    verbosityLevel >0 && 
      cout << __func__ << " verbosityLevel                   : " << verbosityLevel  << endl;

    bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

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
                                            _psVisible,
                                            G4Colour::Blue(),
                                            _psSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                           );

    double const _supWheel_trgtPS_rIn         = _config.getDouble("supWheel_trgtPS_rIn");
    double const _supWheel_trgtPS_rOut        = _config.getDouble("supWheel_trgtPS_rOut");
    G4Material* _suppWheelMaterial = findMaterialOrThrow(_config.getString("supWheel_trgtPS_material"));
    bool const _supWheel_trgtPS_visible       = _config.getBool("supWheel_trgtPS.visible");
    bool const _supWheel_trgtPS_solid         = _config.getBool("supWheel_trgtPS.solid");

    G4ThreeVector _loclCenter(0.0,0.0,0.0);

    TubsParams suppWheelParams( _supWheel_trgtPS_rIn, _supWheel_trgtPS_rOut, _supWheel_trgtPS_halfLength);
    VolumeInfo suppWheelInfo   = nestTubs( "ProductionTargetSupportWheel",
                                            suppWheelParams,
                                            _suppWheelMaterial,
                                            0,
                                            _loclCenter,
                                            prodTargetMotherInfo,
                                            0,
                                            _supWheel_trgtPS_visible,
                                            G4Colour::Gray(),
                                            _supWheel_trgtPS_solid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                           );

    double const _clamp_supWheel_rIn         = _config.getDouble("clamp_supWheel_rIn");
    G4Material* _clampSpWheelMaterial = findMaterialOrThrow(_config.getString("clamp_supWheel_material"));
    bool const _clamp_supWheel_visible       = _config.getBool("clamp_supWheel.visible");
    bool const _clamp_supWheel_solid         = _config.getBool("clamp_supWheel.solid");

    G4ThreeVector _clampPosR(0.0,0.0,_supWheel_trgtPS_halfLength+_clamp_supWheel_halfLength);

    TubsParams clampSpWheelParams( _clamp_supWheel_rIn, _clamp_supWheel_rOut, _clamp_supWheel_halfLength);
    VolumeInfo clampSpWheelInfoR   = nestTubs( "ClampSupportWheel_R",
                                            clampSpWheelParams,
                                            _clampSpWheelMaterial,
                                            0,
                                            _clampPosR,
                                            prodTargetMotherInfo,
                                            0,
                                            _clamp_supWheel_visible,
                                            G4Colour::Gray(),
                                            _clamp_supWheel_solid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                           );

//    G4ThreeVector _clampPosL(0.0,0.0,-(_supWheel_trgtPS_halfLength+_clamp_supWheel_halfLength));
//
//    VolumeInfo clampSpWheelInfoL   = nestTubs( "ClampSupportWheel_L",
//                                            clampSpWheelParams,
//                                            _clampSpWheelMaterial,
//                                            0,
//                                            _clampPosL,
//                                            prodTargetMotherInfo,
//                                            0,
//                                            _clamp_supWheel_visible,
//                                            G4Colour::Gray(),
//                                            _clamp_supWheel_solid,
//                                            forceAuxEdgeVisible,
//                                            placePV,
//                                            doSurfaceCheck
//                                           );

    TubsParams prodTargetParams( 0., tgt->rOut(), tgt->halfLength());

    G4Material* prodTargetMaterial = findMaterialOrThrow(_config.getString("targetPS_materialName"));
    
    bool prodTargetVisible = _config.getBool("targetPS.visible",true);
    bool prodTargetSolid   = _config.getBool("targetPS.solid",true);

    VolumeInfo prodTargetInfo   = nestTubs( "ProductionTarget",
                                            prodTargetParams,
                                            prodTargetMaterial,
                                            &tgt->productionTargetRotation(),
                                            _loclCenter,
                                            prodTargetMotherInfo,
                                            0,
                                            prodTargetVisible,
                                            G4Colour::Magenta(),
                                            prodTargetSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );

    Polycone const & pHubRgtParams = *tgt->getHubsRgtPtr();
    VolumeInfo prodTargetHubRgtInfo  = nestPolycone("ProductionTargetHubRgt",
                                    pHubRgtParams.getPolyconsParams(),
                                    findMaterialOrThrow(pHubRgtParams.materialName()),
                                    &tgt->productionTargetRotation(),
                                    pHubRgtParams.originInMu2e(),
                                    prodTargetMotherInfo,
                                    0,
                                    prodTargetVisible,
                                    G4Colour::Magenta(),
                                    prodTargetSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    Polycone const & pHubLftParams = *tgt->getHubsLftPtr();
    VolumeInfo prodTargetHubLftInfo  = nestPolycone("ProductionTargetHubLft",
                                    pHubLftParams.getPolyconsParams(),
                                    findMaterialOrThrow(pHubLftParams.materialName()),
                                    &tgt->productionTargetRotation(),
                                    pHubLftParams.originInMu2e(),
                                    prodTargetMotherInfo,
                                    0,
                                    prodTargetVisible,
                                    G4Colour::Magenta(),
                                    prodTargetSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
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
