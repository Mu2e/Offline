//
// Free function to create Production Target.
//

// C++ includes
#include <iostream>
#include <sstream>
#include <cmath>

// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.
#include "BeamlineGeom/inc/Beamline.hh"
#include "ProductionSolenoidGeom/inc/ProductionSolenoid.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/Polycone.hh"
#include "GeomPrimitives/inc/Box.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
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
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4Trd.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/GenericFunctions/ASin.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4PVPlacement.hh"

#include "Geant4/G4UnionSolid.hh"
#include "Geant4/G4IntersectionSolid.hh"
#include <cmath>
using namespace std;

namespace mu2e {

  void placeSubtractionVolumeBoxTubs(
                                     const std::string& name //volume info for main object
                                     ,const CLHEP::Hep3Vector& translation
                                     ,VolumeInfo const& parent
                                     ,G4Box*  boxCutout //action volume object
                                     ,G4Tubs* tubsCutout
                                     ,const CLHEP::HepRotation* rotAngle
                                     ,G4Material* material// finish nesting args
                                     ,G4RotationMatrix* rotOverall
                                     ,G4ThreeVector const& offset
                                     ,int const copyNo
                                     ,const G4Colour color
                                     ,const std::string& lookupToken
                                     ,const bool verbose = false
                                     )
  {
    VolumeInfo boxWithTubsCutoutInfo(name,translation,parent.centerInWorld);
    boxWithTubsCutoutInfo.solid = new G4SubtractionSolid(name,boxCutout,tubsCutout, rotOverall, offset);
    finishNesting(boxWithTubsCutoutInfo,material,rotAngle,translation,parent.logical,copyNo,color,lookupToken,verbose);
  }

  void constructTargetPS(VolumeInfo const & parent, SimpleConfig const & _config) {

    Mu2eG4Helper    & _helper = *(art::ServiceHandle<Mu2eG4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

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
      if(verbosityLevel > 0) {
        G4cout << "inside clause 1 of " << __func__ << G4endl;
        G4cout << "target position and hall origin = " << tgt->position() << "\n" << _hallOriginInMu2e << " " <<
          parent.centerInMu2e() << G4endl;
        G4cout << "supWheel and envHalfLength = " << _clamp_supWheel_rOut  << "\n" <<
          envHalfLength << G4endl;
      }

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
      if(verbosityLevel > 0)
        G4cout << "_loclCenter for Tier 1 target = " << _loclCenter << G4endl;
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
        CLHEP::HepRotation* rotFinBase = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
        rotFinBase->rotateX(90.0*CLHEP::degree);
        rotFinBase->rotateZ(90.0*CLHEP::degree);

        CLHEP::HepRotation* rotFin1 = reg.add(CLHEP::HepRotation((*rotFinBase)*tgt->productionTargetRotation()));
        CLHEP::HepRotation* rotFin2 = reg.add(CLHEP::HepRotation((*rotFinBase)*tgt->productionTargetRotation()));
        CLHEP::HepRotation* rotFin3 = reg.add(CLHEP::HepRotation((*rotFinBase)*tgt->productionTargetRotation()));
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
                                            reg.add(CLHEP::HepRotation(tmpRotAxis,rotAngle)),
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
        //CLHEP::HepRotation *tmpSpokeRot = reg.add(CLHEP::HepRotation(tmpRotAxis,rotAngle));

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
                                            reg.add(CLHEP::HepRotation(tmpRotAxis,rotAngle)),
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
      G4Material* prodTargetSupportRingMaterial = findMaterialOrThrow(tgt->supportRingMaterial());

      TubsParams prodTargetMotherParams( 0.
                                         ,tgt->productionTargetMotherOuterRadius()
                                         ,tgt->productionTargetMotherHalfLength());

      G4ThreeVector _loclCenter(0.0,0.0,0.0);
      G4ThreeVector zeroTranslation(0.,0.,0.);
      G4RotationMatrix* targetRotation = reg.add(G4RotationMatrix(tgt->productionTargetRotation().inverse()));
      if (verbosityLevel > 2){G4cout << __PRETTY_FUNCTION__ << "target rotation  = " << *targetRotation << G4endl;}
      VolumeInfo prodTargetMotherInfo   = nestTubs( "ProductionTargetMother",
                                                    prodTargetMotherParams,
                                                    parent.logical->GetMaterial(),
                                                    0,
                                                    tgt->haymanProdTargetPosition() - parent.centerInMu2e(),
                                                    parent,
                                                    0,
                                                    G4Colour::Blue(),
                                                    "PS"
                                                    );

      if (verbosityLevel > 0){
        G4cout << "target position and hall origin = "
                  << tgt->haymanProdTargetPosition() << "\n" <<
          _hallOriginInMu2e << " " << parent.centerInMu2e() << G4endl;
      }
      if (verbosityLevel > 2){
        G4cout << __PRETTY_FUNCTION__ << "created prodTargetMotherInfo "
                  << tgt->productionTargetMotherOuterRadius()
                  << " " <<tgt->productionTargetMotherHalfLength() << G4endl;
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
      if (verbosityLevel > 2){G4cout << __PRETTY_FUNCTION__ << " current Z starts at " << _currentZ << G4endl;}
      //
      // running z locations for placing elements
      CLHEP::Hep3Vector _segmentCenter(0.,0.,0.);
      CLHEP::Hep3Vector _startingSegmentCenter(0.,0.,0.);
      double targetRadius = tgt->rOut();
      //
      // keeping track of these for sensitive volume creation
      int coreCopyNumber = 0;
      int finCopyNumber = 0;
      int finTopCopyNumber = 0;

      double angularSize = asin((tgt->haymanFinThickness()/2.)/targetRadius);
      double dSphi = M_PI/2. - angularSize;
      double dPphi = 2.*angularSize;

      //
      //  need this for a) drawing subtraction volumes and b) rectangular cutouts out of circular support rings, etc. otherwise visualizer gets
      //  confused and how can G4 know which of two volumes some plane, for example, belongs to?
      constexpr double magicOffset = 0.0001;

      //
      //other constants
      double cutoutHeightAboveTarget = tgt->supportRingInnerRadius() + tgt->supportRingCutoutThickness()/2. - magicOffset;

      CLHEP::HepRotation* rotRing = reg.add(CLHEP::HepRotation(tgt->productionTargetRotation()));

      //
      // we're going to use these over and over; save time, memory, make cleaner and put in registry.  save typing with currentFinAngles
      // fins are constructed at 90^o relative to boxes and it's just cleaner to do it this way
      std::vector<CLHEP::HepRotation*> rotFinsG4;
      std::vector<CLHEP::HepRotation*> rotBoxesG4;
      std::vector<double> currentFinAngles;
      for (int ithFin = 0; ithFin < tgt->nHaymanFins(); ++ithFin){
        rotFinsG4.emplace_back(reg.add(CLHEP::HepRotation(tgt->productionTargetRotation())));
        rotBoxesG4.emplace_back(reg.add(CLHEP::HepRotation::IDENTITY));
        //
        // here I have a little confusing fix.  The fin is built along the y-axis.  But the rotation is given wrt the x axis.  Hence I need to
        // subtract off that 90^o . Then the - sign is CLHEP vs G4
        currentFinAngles.emplace_back(tgt->finAngles().at(ithFin));
        rotFinsG4.at(ithFin)->rotateZ(-(-M_PI/2. + currentFinAngles.at(ithFin)));
        rotBoxesG4.at(ithFin)->rotateZ(-(-M_PI/2. + currentFinAngles.at(ithFin)));
        if (verbosityLevel > 1){G4cout << "rotFin " << ithFin << " " << *rotFinsG4.at(ithFin) << G4endl;}
      }

      G4Box* cutoutBox = reg.add(G4Box("cutoutBox",tgt->supportRingCutoutThickness()/2.+ magicOffset,
                                        tgt->haymanFinThickness() + magicOffset,tgt->supportRingCutoutLength()/2.+ magicOffset));


      // get real
      for (int ithSection = 0; ithSection < numberOfSections; ++ithSection){
        int numberOfSegments = tgt->numberOfSegmentsPerSection().at(ithSection);
        std::string startName = "ProductionTargetStartingCoreSection_" + std::to_string(ithSection+1);
        //
        // first place starting segment

        TubsParams startingSegmentParams(0.,targetRadius,tgt->startingSectionThickness().at(ithSection)/2.);
        double currentHalfStartingSegment = tgt->startingSectionThickness().at(ithSection)/2.;
        //
        //set currentZ to be in the center of the starting section
        _currentZ += tgt->startingSectionThickness().at(ithSection)/2.;
        if (verbosityLevel > 3){
          G4cout << __PRETTY_FUNCTION__ << " ithSection , startingSegmentCenter at " << ithSection << " " <<_currentZ << G4endl;
          G4cout << "    and copy number for starting segment is now " << coreCopyNumber << G4endl;
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
          {
            const std::string name = "ProductionTargetFinStartingSection_" + std::to_string(ithSection+1) + "_Fin_" + std::to_string(ithFin);
            if (verbosityLevel > 1)
              {G4cout << "Fin Section and Angle = " << name << " " << currentFinAngles.at(ithFin) << " fin copy number = " << finCopyNumber << G4endl;}

            double finHalfHeightAboveTarget = (tgt->finOuterRadius() - tgt->rOut())/2.;
            G4Box* finBox = new G4Box("finBox",tgt->haymanFinThickness()/2.,finHalfHeightAboveTarget,currentHalfStartingSegment);

            //
            // this tubs is in the xy plane with an extent along z.  So phi goes in the xy plane.
            // need extra length for visualization tool and overlap avoidance
            G4Tubs* tubCutout = new G4Tubs("finCutout",0.,tgt->rOut(),currentHalfStartingSegment + magicOffset,dSphi,dPphi);
            ++finCopyNumber;

            G4ThreeVector finTranslation = currentStartingSegmentCenter;
            //
            // G4 translates and then rotates.  I want to offset the fin to its final location relative to the core and then let G4 perform rotations,
            // both of the fin about the z axis and then the entire production target rotation angle
            double distanceFromCenter = finHalfHeightAboveTarget + tgt->rOut();
            CLHEP::Hep3Vector offsetRelativeToCore(distanceFromCenter*cos(currentFinAngles.at(ithFin))
                                                   ,distanceFromCenter*sin(currentFinAngles.at(ithFin)),0.);
            //
            // finTranslation was already put in the target rotated frame.  I have to rotate the offset vector as well
            offsetRelativeToCore *= tgt->productionTargetRotation().inverse();
            finTranslation = finTranslation + offsetRelativeToCore;

            CLHEP::Hep3Vector downshift = CLHEP::Hep3Vector(0.,-finHalfHeightAboveTarget - targetRadius*cos(angularSize), 0.);
            placeSubtractionVolumeBoxTubs(name
                                          ,finTranslation
                                          ,prodTargetMotherInfo
                                          ,finBox
                                          ,tubCutout
                                          ,rotFinsG4.at(ithFin)
                                          ,prodTargetFinMaterial
                                          ,nullptr
                                          ,downshift
                                          ,finCopyNumber
                                          ,G4Colour::Blue()
                                          ,"PS"
                                          ,verbosityLevel>1);

            //
            // for each fin, build the little section at the top that joins it to the next fin.  This will be  a G4Box with a cutout.
            // Here make the cutout a little tubs without the fancy angular business, since I don't have fins overlapping fins up at the top.
            const std::string nameTop = "ProductionTargetFinTopStartingSection_" + std::to_string(ithSection+1) + "_Fin_" + std::to_string(ithFin);

            double gapRadius = tgt->thicknessOfGapPerSection().at(ithSection)/2.;
            double finTopHalfHeight =( tgt->finOuterRadius() - tgt->heightOfRectangularGapPerSection().at(ithSection))/2.;
            if (verbosityLevel > 2){G4cout << "finTopHalfHeight = " << finTopHalfHeight <<G4endl;}
            G4Box* finTopBox = new G4Box("finTopBox",tgt->haymanFinThickness()/2.,finTopHalfHeight,tgt->thicknessOfGapPerSection().at(ithSection)/2.);
            G4Tubs* finTopCutout = new G4Tubs("finTopCutout",0.,gapRadius + magicOffset,tgt->haymanFinThickness()/2. + magicOffset,0.,2.*M_PI);
            ++finTopCopyNumber;

            if (verbosityLevel > 2){G4cout << "current starting SegmentCenter, current HalfSegment, gap radius = "
                                              << currentStartingSegmentCenter << " " << currentHalfStartingSegment
                                              << " " << gapRadius << G4endl;}
            G4ThreeVector finTopLocation = currentStartingSegmentCenter
                      + tgt->productionTargetRotation().inverse()*CLHEP::Hep3Vector(0.,0.,currentHalfStartingSegment + gapRadius);
            //
            // for debugging
            CLHEP::Hep3Vector offsetTop(0.,0.,0.);
            //
            // rotate the offset vector
            G4ThreeVector finTopTranslation = finTopLocation + offsetTop;
            double distanceTopFromCenter = tgt->finOuterRadius() - finTopHalfHeight;
            CLHEP::Hep3Vector offsetTopRelativeToCore(distanceTopFromCenter*cos(currentFinAngles.at(ithFin))
                                                      ,distanceTopFromCenter*sin(currentFinAngles.at(ithFin)),0.);
            //
            // finTranslation was already put in the target rotated frame.  I have to rotate the offset vector as well
            offsetTopRelativeToCore *= tgt->productionTargetRotation().inverse();
            finTopTranslation = finTopTranslation + offsetTopRelativeToCore;
            if (verbosityLevel > 2){G4cout << "finTopTranslation = " << finTopTranslation << G4endl;}
            CLHEP::Hep3Vector downshiftTopBox = CLHEP::Hep3Vector(0.,-finTopHalfHeight, 0.);
            //
            // these are oriented with gap radius along z.  Recall still in mother volume so this axis is OK
            G4RotationMatrix* tubTopRotation = reg.add(G4RotationMatrix(CLHEP::HepRotation::IDENTITY));

            tubTopRotation->rotateY(M_PI/2.);

            VolumeInfo finTopWithCutoutInfo(nameTop,finTopTranslation,prodTargetMotherInfo.centerInWorld);
            finTopWithCutoutInfo.solid = new G4SubtractionSolid(nameTop,finTopBox,finTopCutout,tubTopRotation,downshiftTopBox);
            finishNesting(finTopWithCutoutInfo
                          ,prodTargetFinMaterial
                          ,rotFinsG4.at(ithFin)
                          ,finTopTranslation
                          ,prodTargetMotherInfo.logical
                          ,finTopCopyNumber
                          ,G4Colour::White()
                          ,"ProductionTarget"
                          ,verbosityLevel>1);
            //
            // this check also uses finWithCutout.  Indirectly it's a way of forcing an overlap check or it won't compile, since finWithCutout is never used otherwise.
            //      if (doSurfaceCheck){checkForOverlaps(finTopWithCutout,_config,0);}

          }

        //
        // for very first starting section, we need to build the support ring
        //
        // build support ring at this end.  I need to build an annulus with cutouts for the fins.
        int boxCopyNumber = 0;
        if (ithSection==0){
          std::string name = "ProductionTargetNegativeEndRing";
          G4Tubs* supportRing = new G4Tubs(name,
                                           tgt->supportRingInnerRadius(),
                                           tgt->supportRingOuterRadius(),
                                           tgt->supportRingLength()/2.,
                                           0.,
                                           2.*M_PI);

          if (verbosityLevel > 0){
            G4cout << __PRETTY_FUNCTION__ << ": \n " << name.c_str() << ":\n  (rin, rout, half length) = ("
                   << tgt->supportRingInnerRadius() << ", " << tgt->supportRingOuterRadius()
                   << ", " << tgt->supportRingLength()/2. << ")" << G4endl;
          }
          // move ring to end of target and then back by cutout size (signs for negative side)
          G4ThreeVector ringTranslation = currentStartingSegmentCenter
            - tgt->productionTargetRotation().inverse()*CLHEP::Hep3Vector(0.,0.,tgt->supportRingLength()/2.);

          if (verbosityLevel > 0){
            G4cout << "  center = (0., 0., -" << (currentStartingSegmentCenter.mag() + tgt->supportRingLength()/2.) << ")" << G4endl;
            G4cout << "  center (wrt target mother) = " << ringTranslation << G4endl;
          }
          // the way we handle the cutout is a little different; the tubs is centered on zero, not offset from the z axis
          // build one cutout box for each fin.  This ignores the extra cutouts (easy enough to put in later) and the mounting mechanisms,
          // which won't matter for muon yield, heat deposition, etc -- just looks like tungsten for those purposes

          //
          // this is all so I can add multiple cutouts to the ring and only do one placement.
          //      std::vector<G4SubtractionSolid*> ringWithCutoutSolid;
          std::vector<G4SubtractionSolid*> ringWithCutoutSolid;

          std::string nameRing = "ProductionTargetNegativeRingCutout";

          for (int ithFin = 0; ithFin < tgt->nHaymanFins(); ++ithFin){
            double currentFinAngle = tgt->finAngles().at(ithFin);
            const std::string name = "ProductionTargetNegativeRingCutout_" + std::to_string(ithSection)
              + "_Segment_" +"_Fin_" + std::to_string(ithFin);
            ++boxCopyNumber;

            if (verbosityLevel > 6)
              {G4cout << "Fin Section and Angle = " << name << " " << currentFinAngle << " fin copy number = " << finCopyNumber << G4endl;}
            if (verbosityLevel > 2){G4cout << " cutoutHeightAboveTarget = " << tgt->supportRingInnerRadius() << " "
                                              << " " << tgt->supportRingCutoutThickness() << " " << cutoutHeightAboveTarget << G4endl;}

            //
            // G4 translates and then rotates.  I want to offset the fin to its final location relative to the core and then let G4 perform rotations,
            // both of the fin about the z axis and then the entire production target rotation angle
            double distanceFromCenter = cutoutHeightAboveTarget;
            CLHEP::Hep3Vector boxShift(distanceFromCenter*cos(currentFinAngles.at(ithFin)),distanceFromCenter*sin(currentFinAngles.at(ithFin)),0.);
            CLHEP::Hep3Vector offsetRelativeToCore(distanceFromCenter*cos(currentFinAngles.at(ithFin)),distanceFromCenter*sin(currentFinAngles.at(ithFin)),
                                                           tgt->supportRingCutoutLength()/2. + magicOffset);

            if (verbosityLevel > 3){
              G4cout << "production target rotation angle = " << tgt->productionTargetRotation().getTheta() << G4endl;
              G4cout << "fin dump: " << currentFinAngles.at(ithFin) << " " << distanceFromCenter << " " << offsetRelativeToCore << G4endl;
              G4cout << "ithfin; current fin angle; rotFin; finTranslation = " << ithFin << " "
                     << currentFinAngles.at(ithFin)  << " " << *rotFinsG4.at(ithFin) << " " << ringTranslation << G4endl;
              G4cout << "production target rotation = " << tgt->productionTargetRotation() << G4endl;
            }

            if (ithFin == 0){
                      ringWithCutoutSolid.push_back(new G4SubtractionSolid(name,supportRing,cutoutBox,rotBoxesG4.at(ithFin),offsetRelativeToCore));
            }else {
                      ringWithCutoutSolid.push_back(new G4SubtractionSolid(name,ringWithCutoutSolid.at(ithFin-1),cutoutBox,rotBoxesG4.at(ithFin),offsetRelativeToCore));
            }
          }


          if (verbosityLevel > 0){
            G4cout << "  center (at end) = (0., 0., " << ringTranslation.mag() << ")" <<G4endl;
            G4cout << "  center (wrt target mother at end) = " << ringTranslation << G4endl;

            G4cout << "  rotation = " << rotRing << G4endl;
          }

          VolumeInfo ringWithCutoutNegative(name,ringTranslation,prodTargetMotherInfo.centerInWorld);
          ringWithCutoutNegative.solid = ringWithCutoutSolid.at(tgt->nHaymanFins() - 1);
          finishNesting(ringWithCutoutNegative
                        ,prodTargetSupportRingMaterial
                        ,rotRing
                        ,ringTranslation
                        ,prodTargetMotherInfo.logical
                        ,finCopyNumber
                        ,G4Colour::White()
                        ,"ProductionTarget"
                        ,verbosityLevel>1);


        }


        /*    -------------------end of support ring, first starting section  ---------------------*/

        //
        // set current z to be at the end of this starting segment before beginning loop on sections
        _currentZ += tgt->startingSectionThickness().at(ithSection)/2.;
        if (verbosityLevel > 0){G4cout << __PRETTY_FUNCTION__ << " ithSection , startingSegment Z ends at " << ithSection << " " <<_currentZ << G4endl;}
        double currentGap = tgt->thicknessOfGapPerSection().at(ithSection);
        double currentHalfSegment = tgt->thicknessOfSegmentPerSection().at(ithSection)/2.;

        for (int ithSegment = 0; ithSegment < numberOfSegments; ++ithSegment){
          std::string name = "ProductionTargetCoreSection_" + std::to_string(ithSection+1)
            + "_Segment_" + std::to_string(ithSegment);
          if (verbosityLevel > 0) {
            G4cout << __PRETTY_FUNCTION__ << "name = " << name << " and core copy number = " << coreCopyNumber << G4endl;
            G4cout << __PRETTY_FUNCTION__ << "gap, and segment half = " << currentGap << " " << currentHalfSegment << G4endl;
          }
          TubsParams segmentParams(0.,targetRadius,currentHalfSegment);
          _currentZ += currentHalfSegment + currentGap;
          if (verbosityLevel > 0){
            G4cout << __PRETTY_FUNCTION__ << " ithSection , current Z for segment center is at " << ithSection << " "
                      << ithSegment << " "  <<_currentZ << G4endl;
          }
          _segmentCenter.setZ(_currentZ);
          G4ThreeVector currentSegmentCenter = tgt->productionTargetRotation().inverse()*_segmentCenter;
          //CLHEP::Hep3Vector currentSegmentCenter = tgt->productionTargetRotation()*_segmentCenter;
          //      CLHEP::Hep3Vector currentSegmentCenter = _segmentCenter;
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
          if (verbosityLevel > 2){G4cout << "segment center after volumeinfo = " << _segmentCenter << G4endl;}

          //
          //now add fins surrounding this core segment.  Build them as simple rectangles with a G4 subtraction volume for the core.

          for (int ithFin = 0; ithFin < tgt->nHaymanFins(); ++ithFin)
            //                                    for (int ithFin = 0; ithFin < 1; ++ithFin)
            {
              const std::string name = "ProductionTargetFinSection_" + std::to_string(ithSection+1) + "_Segment_" + std::to_string(ithSegment) + "_Fin_" + std::to_string(ithFin);
              if (verbosityLevel > 1)
                {G4cout << "Fin Section and Angle = " << name << " " << currentFinAngles.at(ithFin) << " fin copy number = " << finCopyNumber << G4endl;}
              // divide by 2 to make box half-height since I gave it the "radius" to start with; arbitrary convention.
              double finHalfHeightAboveTarget = (tgt->finOuterRadius() - tgt->rOut())/2.;
              if (verbosityLevel > 2){G4cout << "finHeightAboveTarget = " << tgt->finOuterRadius() << " "
                                                 << tgt->rOut() << " " << finHalfHeightAboveTarget << G4endl;}
              G4Box* finBox = new G4Box("finBox",tgt->haymanFinThickness()/2.,finHalfHeightAboveTarget,currentHalfSegment);
              //
              // this tubs is in the xy plane with an extent along z.  So phi goes in the xy plane.
              G4Tubs* tubCutout = new G4Tubs("finCutout",0.,tgt->rOut(),currentHalfSegment + magicOffset,dSphi,dPphi);// need extra length for visualization tool
              //
              // I now have two G4Solids, a box and a cutout.  Combine them to make a logical volume. the box is
              // oriented so that it is "z" thick and "radius" high.  the tubs is made to be the same, so no rotation
              // and also no translation needed

              ++finCopyNumber;

              //
              //for debugging
              CLHEP::Hep3Vector offset(0.,0.,0.);
              G4ThreeVector finTranslation = currentSegmentCenter;

              double distanceFromCenter = finHalfHeightAboveTarget + tgt->rOut();
              CLHEP::Hep3Vector offsetRelativeToCore(distanceFromCenter*cos(currentFinAngles.at(ithFin)),distanceFromCenter*sin(currentFinAngles.at(ithFin)),0.);
              //
              // finTranslation was already put in the target rotated frame.  I have to rotate the offset vector as well
              offsetRelativeToCore *= tgt->productionTargetRotation().inverse();
              finTranslation = finTranslation + offsetRelativeToCore;


              CLHEP::Hep3Vector downshift = CLHEP::Hep3Vector(0.,-finHalfHeightAboveTarget - targetRadius*cos(angularSize), 0.);
              placeSubtractionVolumeBoxTubs(name
                                            ,finTranslation
                                            ,prodTargetMotherInfo
                                            ,finBox
                                            ,tubCutout
                                            ,rotFinsG4.at(ithFin)
                                            ,prodTargetFinMaterial
                                            ,nullptr
                                            ,downshift
                                            ,finCopyNumber
                                            ,G4Colour::Blue()
                                            ,"PS"
                                            ,verbosityLevel>1);
            //
            // for each fin, build the little section at the top that joins it to the next fin.  This will be  a G4Box with a cutout.
            // Here make the cutout a little tubs without the fancy angular business, since I don't have fins overlapping fins up at the top.
            const std::string nameTop = "ProductionTargetFinTopSection_" + std::to_string(ithSection+1)
              + "_Segment_" + std::to_string(ithSegment) + "_Fin_" + std::to_string(ithFin);

            //Much easier to visualize but the code is equally yucky
            double gapRadius = tgt->thicknessOfGapPerSection().at(ithSection)/2.;
            double finTopHalfHeight =( tgt->finOuterRadius() - tgt->heightOfRectangularGapPerSection().at(ithSection))/2.;
            if (verbosityLevel > 2){G4cout << "finTopHalfHeight = " << finTopHalfHeight <<G4endl;}
            G4Box* finTopBox = new G4Box("finTopBox",tgt->haymanFinThickness()/2.,finTopHalfHeight,tgt->thicknessOfGapPerSection().at(ithSection)/2.);
            G4Tubs* finTopCutout = new G4Tubs("finTopCutout",0.,gapRadius + magicOffset,tgt->haymanFinThickness()/2. + magicOffset,0.,2.*M_PI);
            ++finTopCopyNumber;

            if (verbosityLevel > 2){G4cout << "current starting SegmentCenter, current HalfSegment, gap radius = "
                                              << currentStartingSegmentCenter << " " << currentHalfStartingSegment << " "
                                              << gapRadius << G4endl;}
            G4ThreeVector finTopBoxLocation = currentSegmentCenter
                      + tgt->productionTargetRotation().inverse()*CLHEP::Hep3Vector(0.,0.,currentHalfSegment + gapRadius);
            //
            // for debugging
            CLHEP::Hep3Vector offsetTop(0.,0.,0.);
            //
            // rotate the offset vector
            G4ThreeVector finTopTranslation = finTopBoxLocation + offsetTop;
            //
            // G4 translates and then rotates.  I want to offset the fin to its final location relative to the core and then let G4 perform rotations,
            // both of the fin about the z axis and then the entire production target rotation angle
            double distanceTopFromCenter = tgt->finOuterRadius() - finTopHalfHeight;

            CLHEP::Hep3Vector offsetTopRelativeToCore(distanceTopFromCenter*cos(currentFinAngles.at(ithFin)),distanceTopFromCenter*sin(currentFinAngles.at(ithFin)),0.);
            //
            // finTranslation was already put in the target rotated frame.  I have to rotate the offset vector as well
            offsetTopRelativeToCore *= tgt->productionTargetRotation().inverse();
            finTopTranslation = finTopTranslation + offsetTopRelativeToCore;
            if (verbosityLevel > 2){G4cout << "finTopTranslation = " << finTopTranslation << G4endl;}
            CLHEP::Hep3Vector downshiftTopBox = CLHEP::Hep3Vector(0.,-finTopHalfHeight, 0.);

            //
            // these are oriented with gap radius along z.  Recall still in mother volume so this axis is OK
            G4RotationMatrix* tubTopRotation = reg.add(G4RotationMatrix(CLHEP::HepRotation::IDENTITY));
            tubTopRotation->rotateY(M_PI/2.);

            VolumeInfo finTopWithCutoutInfo(nameTop,finTopTranslation,prodTargetMotherInfo.centerInWorld);
            finTopWithCutoutInfo.solid = new G4SubtractionSolid(nameTop,finTopBox,finTopCutout,tubTopRotation,downshiftTopBox);
            finishNesting(finTopWithCutoutInfo
                          ,prodTargetFinMaterial
                          ,rotFinsG4.at(ithFin)
                          ,finTopTranslation
                          ,prodTargetMotherInfo.logical
                          ,finTopCopyNumber
                          ,G4Colour::White()
                          ,"ProductionTarget"
                          ,verbosityLevel>1);

            }
          _currentZ += currentHalfSegment;
          //
          // if this is the last segment, add another gap before next section starts!
          if (ithSegment == numberOfSegments-1) {
            if (verbosityLevel > 0){G4cout << " z before = " << _currentZ << G4endl;}
            _currentZ += currentGap;
            if (verbosityLevel > 0){G4cout << " z after = " << _currentZ << G4endl;}
          }
          if (verbosityLevel > 0){G4cout << __PRETTY_FUNCTION__ << " ending at z=" << _currentZ << G4endl;}
        }
          ++coreCopyNumber;
        //
        // if this is the last section, add on one more beginning block. decided not to extend vectors and have zero segments, just ugly
        if (ithSection == (numberOfSections -1) ){
          startName = "ProductionTargetStartingCoreSection_" + std::to_string(ithSection+1+1);
          TubsParams startingSegmentParamsEnd(0.,targetRadius,tgt->startingSectionThickness().at(ithSection)/2.);

          //
          //set currentZ to be in the center of the starting section
          _currentZ += tgt->startingSectionThickness().at(ithSection)/2.;
          if (verbosityLevel > 0){G4cout << __PRETTY_FUNCTION__ << " ithSection+1 , startingSegmentCenter at " << ithSection+1 << " " <<_currentZ << G4endl;}
          _startingSegmentCenter.setZ(_currentZ);
          CLHEP::Hep3Vector currentStartingSegmentCenter = tgt->productionTargetRotation().inverse()*_startingSegmentCenter;
          double currentHalfStartingSegment = tgt->startingSectionThickness().at(ithSection)/2.;
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
          for (int ithFin = 0; ithFin < tgt->nHaymanFins(); ++ithFin){

            double currentFinAngle = tgt->finAngles().at(ithFin);
            const std::string name = "ProductionTargetFinStartingSection_" + std::to_string(ithSection+1+1) + "_Fin_" + std::to_string(ithFin);

            if (verbosityLevel > 1)
              {G4cout << "Fin Section and Angle = " << name << " " << currentFinAngle << " fin copy number = " << finCopyNumber << G4endl;}

            double finHalfHeightAboveTarget = (tgt->finOuterRadius() - tgt->rOut())/2.;
            G4Box* finBox = new G4Box("finBox",tgt->haymanFinThickness()/2.,finHalfHeightAboveTarget,currentHalfStartingSegment);
            //
            // this tubs is in the xy plane with an extent along z.  So phi goes in the xy plane.
            // need extra length for visualization tool and overlap avoidance
            G4Tubs* tubCutout = new G4Tubs("finCutout",0.,tgt->rOut(),currentHalfStartingSegment + magicOffset,dSphi,dPphi);

            CLHEP::HepRotation* rotFin = reg.add(CLHEP::HepRotation(tgt->productionTargetRotation().inverse()));
            CLHEP::HepRotation* rotFinG4 = reg.add(CLHEP::HepRotation(tgt->productionTargetRotation()));
            //
            // ok here I have a little confusing fix.  The fin is built along the y-axis.  But the rotation is given wrt the x axis.  Hence I need to
            // subtract off that 90^o
            double rotAdj = -M_PI/2 + currentFinAngle;
            rotFin->rotateZ(rotAdj);
            //
            // g4 version.
            rotFinG4->rotateZ(-rotAdj);
            ++finCopyNumber;

            //
            // for debugging
            CLHEP::Hep3Vector offset(0.,0.,0.);
              //
              // rotate the offset vector
            G4ThreeVector finTranslation = currentStartingSegmentCenter;
            //
            // G4 translates and then rotates.  I want to offset the fin to its final location relative to the core and then let G4 perform rotations,
            // both of the fin about the z axis and then the entire production target rotation angle
            double distanceFromCenter = finHalfHeightAboveTarget + tgt->rOut();
            CLHEP::Hep3Vector offsetRelativeToCore(distanceFromCenter*cos(currentFinAngle),distanceFromCenter*sin(currentFinAngle),0.);
            //
            // finTranslation was already put in the target rotated frame.  I have to rotate the offset vector as well
            offsetRelativeToCore *= tgt->productionTargetRotation().inverse();
            finTranslation = finTranslation + offsetRelativeToCore;
            if (verbosityLevel > 2){G4cout << "finTranslation = " << finTranslation << G4endl;}
            CLHEP::Hep3Vector downshift = CLHEP::Hep3Vector(0.,-finHalfHeightAboveTarget - targetRadius*cos(angularSize), 0.);

            VolumeInfo finWithCutoutInfo(name,finTranslation,prodTargetMotherInfo.centerInWorld);
            finWithCutoutInfo.solid = new G4SubtractionSolid(name,finBox,tubCutout,nullptr,downshift);
            finishNesting(finWithCutoutInfo
                          ,prodTargetFinMaterial
                          ,rotFinsG4.at(ithFin)
                          ,finTranslation
                          ,prodTargetMotherInfo.logical
                          ,finCopyNumber
                          ,G4Colour::White()
                          ,"ProductionTarget"
                          ,verbosityLevel>1);

          }
          //
          // set current z to be at the end of this starting segment as a final check
          _currentZ += tgt->startingSectionThickness().at(ithSection)/2.;
          if (verbosityLevel > 0){
            G4cout << __PRETTY_FUNCTION__ << " ithSection+1 , startingSegment Z ends at " << ithSection+1 << " " <<_currentZ << G4endl;
            G4cout << "             with copy number for last segment  = " << ithSection+1 << G4endl;
          }
          /***********and now the final end ring  ***********/
          std::string name = "ProductionTargetPositiveEndRing";
          G4Tubs* supportRing = new G4Tubs(name,
                                           tgt->supportRingInnerRadius(),
                                           tgt->supportRingOuterRadius(),
                                           tgt->supportRingLength()/2.,
                                           0.,
                                           2.*M_PI);
          if (verbosityLevel > 0){
            G4cout << __PRETTY_FUNCTION__ << ": \n " << name.c_str() << ":\n  (rin, rout, half length) = ("
                   << tgt->supportRingInnerRadius() << ", " << tgt->supportRingOuterRadius()
                   << ", " << tgt->supportRingLength()/2. << ")" << G4endl;
          }
          // move ring to end of target and then back by cutout size (signs for positive side)
          G4ThreeVector ringTranslation = currentStartingSegmentCenter
            + tgt->productionTargetRotation().inverse()*CLHEP::Hep3Vector(0.,0.,tgt->supportRingLength()/2.);

          if (verbosityLevel > 0){
            G4cout << "  center = (0., 0., " << (currentStartingSegmentCenter.mag() + tgt->supportRingLength()/2.) << ")" << G4endl;
            G4cout << "  center (wrt target mother) = " << ringTranslation << G4endl;
          }
          std::vector<G4SubtractionSolid*> ringWithCutoutSolid;
          std::string nameRing = "ProductionTargetPositiveRingCutout";

          for (int ithFin = 0; ithFin < tgt->nHaymanFins(); ++ithFin){
            double currentFinAngle = tgt->finAngles().at(ithFin);
            const std::string name = "ProductionTargetPositiveRingCutout_" + std::to_string(ithSection)
              + "_Segment_" +"_Fin_" + std::to_string(ithFin);
            ++boxCopyNumber;


            //
            // G4 translates and then rotates.  I want to offset the fin to its final location relative to the core and then let G4 perform rotations,
            // both of the fin about the z axis and then the entire production target rotation angle
            double distanceFromCenter = cutoutHeightAboveTarget;

            CLHEP::Hep3Vector boxShift(distanceFromCenter*cos(currentFinAngle),distanceFromCenter*sin(currentFinAngle),0.);
            CLHEP::Hep3Vector offsetRelativeToCore(distanceFromCenter*cos(currentFinAngle),distanceFromCenter*sin(currentFinAngle),
                                                   -tgt->supportRingCutoutLength()/2. - magicOffset);// note flipped sign from other end of target

            CLHEP::Hep3Vector downshift = CLHEP::Hep3Vector(0.,0.,0.);//for debugging
            offsetRelativeToCore = offsetRelativeToCore + downshift;
            if (ithFin == 0){
              ringWithCutoutSolid.push_back(new G4SubtractionSolid(name,supportRing,cutoutBox,rotBoxesG4.at(ithFin),offsetRelativeToCore));
            }else {
              ringWithCutoutSolid.push_back(new G4SubtractionSolid(name,ringWithCutoutSolid.at(ithFin-1),cutoutBox,rotBoxesG4.at(ithFin),offsetRelativeToCore));
            }
          }


          if (verbosityLevel > 0){
            G4cout << "  center (at end) = (0., 0., " << ringTranslation.mag() << ")" <<G4endl;
            G4cout << "  center (wrt target mother at end) = " << ringTranslation << G4endl;

            G4cout << "  rotation = " << rotRing << G4endl;
          }
          VolumeInfo ringWithCutoutPositive(name,ringTranslation,prodTargetMotherInfo.centerInWorld);
          ringWithCutoutPositive.solid = ringWithCutoutSolid.at(tgt->nHaymanFins() - 1); // what does this =  mean?
          finishNesting(ringWithCutoutPositive
                        ,prodTargetSupportRingMaterial
                        ,rotRing
                        ,ringTranslation
                        ,prodTargetMotherInfo.logical
                        ,finCopyNumber
                        ,G4Colour::White()
                        ,"ProductionTarget"
                        ,verbosityLevel>1);

        } //end if(ithSection == (numberOfSections -1) )
      } //end for(int ithSection...)

      //Add support structures for the production target
      if(tgt->supportsBuild()) {
        G4Material* suppWheelMaterial = findMaterialOrThrow(tgt->supportWheelMaterial());
        G4ThreeVector localWheelCenter(0.0,0.0,0.0); //no offset
        double suppWheelParams[] = {tgt->supportWheelRIn(), tgt->supportWheelROut(), tgt->supportWheelHL()};
        //create the volume info for the support wheel+rods
        VolumeInfo suppWheelInfo( "ProductionTargetSupportWheel", localWheelCenter, prodTargetMotherInfo.centerInMu2e());
        suppWheelInfo.solid = new G4Tubs("ProductionTargetSupportWheel_wheel", suppWheelParams[0], suppWheelParams[1],
                                         suppWheelParams[2], 0., CLHEP::twopi);
                                               // suppWheelParams,
                                               // suppWheelMaterial,
                                               // 0,
                                               // localWheelCenter,
                                               // prodTargetMotherInfo,
                                               // 0,
                                               // G4Colour::Gray(),
                                               // "PS"
                                               // );

        // add spokes //

        //spoke info
        const int nspokesperside = tgt->nSpokesPerSide();
        G4Material* spokeMaterial = findMaterialOrThrow(tgt->spokeMaterial());
        //target info
        double rTarget = tgt->supportRingOuterRadius(); //radius of the support ring to attach to
        double zTarget = tgt->halfHaymanLength(); //where along the target to attach
        double smallGap = 0.001; //for adding small offsets to avoid overlaps due to precision
        //initialize parameter vectors
        //features on wheel
        const vector<double> supportWheelFeatureAngles = tgt->supportWheelFeatureAngles();
        const vector<double> supportWheelFeatureArcs   = tgt->supportWheelFeatureArcs  ();
        const vector<double> supportWheelFeatureRIns   = tgt->supportWheelFeatureRIns  ();
        //support rods in wheel
        const vector<double> supportWheelRodHL           = tgt->supportWheelRodHL          ();
        const vector<double> supportWheelRodOffset       = tgt->supportWheelRodOffset      ();
        const vector<double> supportWheelRodRadius       = tgt->supportWheelRodRadius      ();
        const vector<double> supportWheelRodRadialOffset = tgt->supportWheelRodRadialOffset();
        const vector<double> supportWheelRodWireOffsetD  = tgt->supportWheelRodWireOffsetD ();
        const vector<double> supportWheelRodWireOffsetU  = tgt->supportWheelRodWireOffsetU ();
        const vector<double> supportWheelRodAngles       = tgt->supportWheelRodAngles      ();
        //spoke (support wire) angles
        const vector<double> spokeTargetAnglesD = tgt->spokeTargetAnglesD();
        const vector<double> spokeTargetAnglesU = tgt->spokeTargetAnglesU();
        if(verbosityLevel > 0)
          std::cout << "Printing information about production target supports:\n";

        const double targetAngle = tgt->rotHaymanY(); //assume target angle is only in the x-z plane for supports
        CLHEP::HepRotation* rodRot = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
        rodRot->rotateY(-1.*targetAngle);

        for(int istream = 0; istream < 2; ++istream) {
          for(int ispoke = 0; ispoke < nspokesperside; ++ispoke) {
            const double wheelAngle =  supportWheelRodAngles[ispoke]*CLHEP::degree;
            //get angle of the support rod on the wheel and the angle on the target the wire connects to
            const double targetWireAngle = (istream == 0) ? spokeTargetAnglesD[ispoke]*CLHEP::degree
              : spokeTargetAnglesU[ispoke]*CLHEP::degree;
            double rWheel = supportWheelRodRadialOffset[ispoke]; // radius of the wheel to attach to
            CLHEP::Hep3Vector rodCenter(rWheel*cos(wheelAngle), rWheel*sin(wheelAngle), 0.);
            const double rodOffset = supportWheelRodOffset[ispoke];
            rodCenter += CLHEP::Hep3Vector(sin(targetAngle)*rodOffset, 0., cos(targetAngle)*rodOffset);
            if(istream == 0) { //only do once
              //add the features near the support rods in the bicycle wheel
              const double featureAngle = supportWheelFeatureAngles[ispoke]*CLHEP::degree; //angle of feature center
              const double featureArc   = supportWheelFeatureArcs[ispoke]*CLHEP::degree; //width in angle
              const double featureRIn   = supportWheelFeatureRIns[ispoke]; //inner radius of feature
              const double featureROut = tgt->supportWheelRIn() + smallGap; //ensure they overlap for union
              // double featureR = (featureRIn + featureROut)/2.; //radius of feature center
              // CLHEP::Hep3Vector featureCenter(featureR*cos(featureAngle), featureR*sin(featureAngle), 0.);
              CLHEP::Hep3Vector featureCenter(localWheelCenter); //center is wheel center
              double featureParams[] = {featureRIn, featureROut, tgt->supportWheelHL(), featureAngle - featureArc/2. /*phi0*/, featureArc /*dphi*/};
              G4Tubs* featureTubs = new G4Tubs("ProductionTargetSupportFeature_" +std::to_string(ispoke),
                                               featureParams[0], featureParams[1], featureParams[2], featureParams[3], featureParams[4]);
              suppWheelInfo.solid = new G4UnionSolid("ProductionTargetSupportWheelFeature_union_"+std::to_string(ispoke),
                                                     suppWheelInfo.solid, featureTubs, 0, featureCenter);
              //add the support rod to the wheel
              G4Tubs* rodTubs = new G4Tubs("ProductionTargetSupportRod_" +std::to_string(ispoke),
                                           0., supportWheelRodRadius[ispoke], supportWheelRodHL[ispoke], 0., CLHEP::twopi);
              suppWheelInfo.solid = new G4UnionSolid("ProductionTargetSupportWheelRod_union_"+std::to_string(ispoke),
                                                     suppWheelInfo.solid, rodTubs, rodRot, rodCenter);
            }
            const int side = (1-2*istream); //+1 or -1
            //info about wire connection
            //get end of the rod on this side
            CLHEP::Hep3Vector rodAxis(sin(targetAngle), 0., cos(targetAngle));
            CLHEP::Hep3Vector wheelPos(rodCenter);
            wheelPos += side*supportWheelRodHL[ispoke]*rodAxis;
            //translate from rod center to edge
            CLHEP::Hep3Vector rodCenterToWire(cos(wheelAngle)*cos(targetAngle),
                                              sin(wheelAngle)*cos(targetAngle),
                                              -cos(wheelAngle)*sin(targetAngle));
            wheelPos -= supportWheelRodRadius[ispoke]*rodCenterToWire;
            double zWireRodOffset = (istream == 0) ? supportWheelRodWireOffsetD[ispoke] : supportWheelRodWireOffsetU[ispoke];
            wheelPos -= side*zWireRodOffset*rodAxis;

            //get wire position on target
            CLHEP::Hep3Vector targetPos(rTarget*cos(targetWireAngle), rTarget*sin(targetWireAngle), side*zTarget);
            targetPos = tgt->productionTargetRotation().inverse()*targetPos; //rotate from target frame to mother frame
            CLHEP::Hep3Vector spokeAxis((wheelPos-targetPos).unit());
            CLHEP::Hep3Vector targetAxis(0.,0.,side);
            targetAxis = tgt->productionTargetRotation().inverse()*targetAxis;
            CLHEP::Hep3Vector zax(0.,0.,1.);
            if(verbosityLevel > 0)
              std::cout << "istream " << istream << " ispoke " << ispoke << std::endl
                        << "target pos " << targetPos << "\nwheel pos " << wheelPos << std::endl
                        << "Target axis " << targetAxis << "\nSpoke axis " << spokeAxis << std::endl
                        << "Rod axis " << rodAxis << "\nRod center to wire axis " << rodCenterToWire << std::endl;
            //to remove overlaps where the wire connects, need angle of wire and surface connecting to
            //remove overlap at target
            double wireTargetAngle = targetAxis.angle(-1.*spokeAxis);
            double deltaLength = (abs(tan(wireTargetAngle)) > 1.e-6) ? abs(tgt->spokeRadius()/tan(wireTargetAngle)) : 0.; //give up if ~paralle
            targetPos += (deltaLength+0.1)*spokeAxis; //subtract off the length
            if(verbosityLevel > 0)
              std::cout << "wire target angle " << wireTargetAngle << " delta L " << deltaLength
                        << " target pos " << targetPos <<std::endl;

            //next remove overlap at rod
            //
            double wireRodAngle = abs((rodCenterToWire).angle(spokeAxis));
            deltaLength = abs(tan(wireRodAngle)/tgt->spokeRadius());
            wheelPos -= (deltaLength+1.)*spokeAxis;
            if(verbosityLevel > 0)
              std::cout << "wire rod angle " << wireRodAngle << " delta L " << deltaLength
                        << " wheel pos " << wheelPos <<std::endl;

            CLHEP::Hep3Vector spokeCenter((wheelPos+targetPos)/2.);
            double spokeLength = abs((wheelPos-targetPos).mag());
            TubsParams spokeParams(0., tgt->spokeRadius(), 0.5*spokeLength);
            CLHEP::HepRotation* spokeRot = reg.add(CLHEP::HepRotation(spokeAxis.cross(zax), spokeAxis.angle(zax)));
            std::stringstream spokeName;
            spokeName << "ProductionTargetSpokeWire_" ;
            if(istream == 0)
              spokeName << "Downstream_";
            else
              spokeName << "Upstream_";
            spokeName << ispoke;

            VolumeInfo spokeInfo   = nestTubs( spokeName.str(),
                                               spokeParams,
                                               spokeMaterial,
                                               spokeRot,
                                               spokeCenter,
                                               prodTargetMotherInfo,
                                               0,
                                               G4Colour::Gray(),
                                               "PS"
                                               );

          } //end spokes loop
        } //end stream loop
        finishNesting(suppWheelInfo,
                      suppWheelMaterial,
                      0,
                      localWheelCenter,
                      prodTargetMotherInfo.logical,
                      0,
                      G4Colour::Gray(),
                      "PS"
                      );

      } //end adding support structures
    } //end ProductionTargetMaker::hayman_v_2_0
  } //end constructTargetPS
} //end namespace mu2e

 // end Mu2eWorld::constructPS
