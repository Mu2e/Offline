/*
 * ITrackerBuilder.cpp
 *
 *  Created on: Feb 11, 2010
 *      Author: tassiell
 */

#include "Mu2eG4/inc/ITrackerBuilder.hh"

#include "G4Hype.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4Sphere.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VisAttributes.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/ITGasLayerSD.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "cetlib/exception.h"
#include "cetlib/pow.h"
#include "globals.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <boost/regex.hpp>
#include <iostream>
#include <sstream>
//#include "Mu2eG4/inc/ITGasLayerSD_ExtWireData.hh"
//#include "Mu2eG4/inc/ITGasLayerSD_v2.hh"
//#include "Mu2eG4/inc/ITGasLayerSD_v3.hh"

using namespace std;

using cet::pow;

namespace mu2e {

bool checkOverlap, detailedCheck;

VolumeInfo ITrackerBuilder::constructTracker( G4LogicalVolume* mother, double zOff ){

        // Master geometry for the ITracker.
        GeomHandle<ITracker> itracker;
        art::ServiceHandle<GeometryService> geom;
        SimpleConfig const& config  = geom->config();

        VolumeInfo trackerInfo;

        // Make the mother volume for the ITracker.
        string trackerName("TrackerMother");

        double z0    = CLHEP::mm * itracker->z0();
        G4ThreeVector trackerOffset(0.,0.,z0-zOff);

        checkOverlap = config.getBool("g4.doSurfaceCheck",false);
        detailedCheck = checkOverlap&&config.getBool("itracker.doDetailedSurfCheck",false);

        /*      if ( itracker->rOut() >= ((G4Tubs *)mother->GetSolid())->GetOuterRadius() )
          throw cet::exception("GEOM") <<"The ITracker doesn't fit inside the DS, please check the external radius of the ITracker\n"
                                       <<"rout ITracker="<<itracker->rOut()<<" DS outer radius="<<((G4Tubs *)mother->GetSolid())->GetOuterRadius()<<"\n"
                                       <<(*(G4Tubs *)mother->GetSolid());
        */
        if (itracker->isExternal()) {
                throw cet::exception("GEOM") <<"This GDML file option is temporarily disabled\n";
                //            G4GDMLParser parser;
                //            parser.Read(itracker->extFile().c_str());
                //            trackerInfo.logical = parser.GetWorldVolume()->GetLogicalVolume();
                //            trackerInfo.logical->SetName(trackerName.c_str());
                //        trackerInfo.physical =  new G4PVPlacement( 0,
                //                                                   trackerOffset,
                //                                                   trackerInfo.logical,
                //                                                   trackerName,
                //                                                   mother,
                //                                                   0,
                //                                                   0);
                //
                //        // Visualization attributes of the the mother volume.
                //        {
                //          G4VisAttributes* visAtt = new G4VisAttributes(true, G4Colour::Red() );
                //          visAtt->SetForceSolid(true);
                //          visAtt->SetForceAuxEdgeVisible (false);
                //          visAtt->SetVisibility(true);
                //          visAtt->SetDaughtersInvisible(true);
                //          trackerInfo.logical->SetVisAttributes(visAtt);
                //        }
                //
                //        int nDaugVol = trackerInfo.logical->GetNoDaughters();
                //        string iVolName;
                //        VolumeInfo gasLayerInfo;
                //
                //        boost::regex vaolSyntax("^[gw]volS[0-9]{2}R[0-9]{2}_");
                //        boost::match_results<std::string::const_iterator> matchingStr;
                //
                //            G4SDManager* SDman   = G4SDManager::GetSDMpointer();
                //            std::vector<int> submatches;
                //            submatches.push_back(-1);
                //        for ( int iVol = 0;        iVol<nDaugVol; ++iVol) {
                //                gasLayerInfo.logical = trackerInfo.logical->GetDaughter(iVol)->GetLogicalVolume();
                //                iVolName=gasLayerInfo.logical->GetName();
                //
                //                /*
                //                boost::sregex_token_iterator i(iVolName.begin(), iVolName.end(), vaolSyntax, -1);
                //            boost::sregex_token_iterator j;
                //
                //            unsigned count = 0;
                //            while(i != j)
                //            {
                //               cout << *i++ << endl;
                //               count++;
                //            }
                //            cout << "There were " << count << " tokens found." << endl;
                //            */
                //
                //            if ( boost::regex_search(iVolName,matchingStr,vaolSyntax) ){
                //                        G4String cellSDname = matchingStr.str(0).substr(0,matchingStr.str(0).size()-1);
                //                         cout<<cellSDname<<endl;
                //                        ITGasLayerSD* glSD     = new ITGasLayerSD_ExtWireData( cellSDname );
                //                        SDman->AddNewDetector( glSD );
                //                        gasLayerInfo.logical->SetSensitiveDetector( glSD );
                //                }
                //        }
        } else {

                G4VisAttributes* visAtt = new G4VisAttributes(true, G4Colour::White() );
                visAtt->SetForceSolid(true);
                visAtt->SetForceAuxEdgeVisible (false);
                visAtt->SetVisibility(false);
                visAtt->SetDaughtersInvisible(false);

                char shape[30], vol[30], shape_name_FD[30], shape_name_SD[30], vol_name_FD[30], vol_name_SD[30], wire_name[40];
                sprintf(shape_name_FD,"tube_Field");
                sprintf(shape_name_SD,"tube_Sense");

                int superlayer =0,iring=0;
                //G4SDManager* SDman   = G4SDManager::GetSDMpointer();

                G4ThreeVector positionTracker = G4ThreeVector(0,0,0);
                //        boost::shared_array<SuperLayer> SLayers = itracker->getSuperLayersArray();

                G4VisAttributes* visAttWall = new G4VisAttributes(true, G4Colour::Red() );
                visAttWall->SetForceSolid(true);
                visAttWall->SetForceAuxEdgeVisible (false);
                visAttWall->SetVisibility(true);
                visAttWall->SetDaughtersInvisible(true);

                G4VisAttributes* visAttGas = new G4VisAttributes(true, G4Colour::Cyan() );
                visAttGas->SetForceSolid(true);
                visAttGas->SetForceAuxEdgeVisible (false);
                visAttGas->SetForceLineSegmentsPerCircle(12);
                visAttGas->SetVisibility(itracker->displayGasLayer());
                visAttGas->SetDaughtersInvisible(!itracker->displayWires());

                G4VisAttributes* visAttFw = new G4VisAttributes(true, G4Colour::Grey() );
                visAttFw->SetForceSolid(true);
                visAttFw->SetForceAuxEdgeVisible (false);
                visAttFw->SetVisibility(itracker->displayWires());
                visAttFw->SetDaughtersInvisible(true);

                G4VisAttributes* visAttSw = new G4VisAttributes(true, G4Colour::Yellow() );
                visAttSw->SetForceSolid(true);
                visAttSw->SetForceAuxEdgeVisible (false);
                visAttSw->SetVisibility(itracker->displayWires());
                visAttSw->SetDaughtersInvisible(true);

//                for (int iWall = 0; iWall < itracker->getNWalls(); ++iWall) {
//                        VolumeInfo WallInfo;
//                        boost::shared_ptr<Wall> iwall = itracker->getWall(iWall);
//                        WallInfo = buildWall(iwall.get(),itracker->endcapType());
//                        WallInfo.logical->SetVisAttributes(visAttWall);
//                        WallInfo.physical = new G4PVPlacement(iwall->getPos(),
//                                        WallInfo.logical,                // its logical volume
//                                        iwall->getName(),                // its name
//                                        WallInfo.logical,                // its mother  volume
//                                        false,                                        // no boolean operations
//                                        0);                                                // copy number
//                }

                double outerWallInnerRadius=0.0;
                double innerWallOuterRadius=0.0;
                double planeEndCapThikness=0.0;
                bool   reshapeWall=false;
                double innerWallHaltLength=0.0;
                double outerWallInnerRadiusWithElec=0.0;

                multimap<Wall::Walltype,boost::shared_ptr<Wall> >::iterator walls_it;
                for ( walls_it=itracker->getWalls()->begin() ; walls_it != itracker->getWalls()->end(); walls_it++ ) {
                        boost::shared_ptr<Wall> iwall = walls_it->second;
                        if (iwall->getType() == Wall::outer) {
                                outerWallInnerRadius = iwall->getRmin();
                                outerWallInnerRadiusWithElec = outerWallInnerRadius;
                        }
                        else if (iwall->getType() == Wall::inner) {
                                innerWallOuterRadius = iwall->getRmax();

                                Wall *wall = iwall.get();
                                int nSub=wall->getNShells();
                                for (int ishell=0; ishell<nSub; ishell++){
                                        G4String iWallShellMat = wall->getMaterialsName()->at(ishell);
                                        if ( iWallShellMat.contains("ITGas") ) {
                                                innerWallOuterRadius-=wall->getThicknesses()->at(ishell);
                                        }
                                }
                                innerWallHaltLength = iwall->_pDz;
                        }
                        else if (iwall->getType() == Wall::endcap) {
                                if ( itracker->endcapType() == ITracker::Plane ) {
                                        reshapeWall=true;
                                        planeEndCapThikness=2.0*iwall->getDz();
                                }
                        }
                }

                double oringWidth  = config.getDouble("itracker.oringWidth",0.0);
                double oringHeight = config.getDouble("itracker.oringHeight",0.0);
                G4Material* oringMat = findMaterialOrThrow( config.getString("itracker.oringMaterial","WAGVacuum") );
                double oringCenter = innerWallHaltLength;
                double travkVolRmax = config.getDouble("itracker.elctContRmax",itracker->rOut());
                bool addElctronicBox = false;
                if (travkVolRmax<itracker->rOut()) { throw cet::exception("GEOM") <<"The tracker volume maximum radius is less than ITracker R out\n"; }
                if (travkVolRmax>itracker->rOut()) { addElctronicBox = true; }

                G4Material* Vacuum = findMaterialOrThrow( "WAGVacuum" );
                if (oringWidth>0.0 && oringWidth>0.0) {
                        trackerInfo.solid = new G4Tubs("Itracker", innerWallOuterRadius-oringHeight-0.001,travkVolRmax,itracker->maxEndCapDim(),0.0,360.0*CLHEP::degree);
                        trackerInfo.logical = new G4LogicalVolume(trackerInfo.solid , Vacuum, trackerName,0,0,0);
                        trackerInfo.logical->SetVisAttributes(visAtt);
                        oringWidth*=0.5;
                        oringCenter += oringWidth;
                        G4VSolid* oringInnWallShape = new G4Tubs("oringInnerWall", innerWallOuterRadius-oringHeight,innerWallOuterRadius,oringWidth,0.0,360.0*CLHEP::degree);
                        G4LogicalVolume*   oringInnWallVol = new G4LogicalVolume(oringInnWallShape,oringMat,"oringInnerWallVol",0,0,0);
                        G4VPhysicalVolume *oringInnWallPos_R = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,oringCenter),
                                        oringInnWallVol,      // its logical volume
                                        "oringInnerWallVol",  // its name
                                        trackerInfo.logical,  // its mother  volume
                                        false,                // no boolean operations
                                        0,                    // copy number
                                        detailedCheck);
                        oringInnWallPos_R->GetCopyNo(); //just to remove the warning during compiling
                        G4VPhysicalVolume *oringInnWallPos_L = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,-oringCenter),
                                        oringInnWallVol,      // its logical volume
                                        "oringInnerWallVol",  // its name
                                        trackerInfo.logical,  // its mother  volume
                                        false,                // no boolean operations
                                        1,                    // copy number
                                        detailedCheck);
                        oringInnWallPos_L->GetCopyNo(); //just to remove the warning during compiling
                } else {
                        trackerInfo.solid = new G4Tubs("Itracker", itracker->r0()-0.001,travkVolRmax,itracker->maxEndCapDim(),0.0,360.0*CLHEP::degree);
                        trackerInfo.logical = new G4LogicalVolume(trackerInfo.solid , Vacuum, trackerName,0,0,0);
                        trackerInfo.logical->SetVisAttributes(visAtt);
                }

                if (reshapeWall) {
                        for ( walls_it=itracker->getWalls()->begin() ; walls_it != itracker->getWalls()->end(); walls_it++ ) {
                                boost::shared_ptr<Wall> iwall = walls_it->second;
                                if (iwall->getType() == Wall::outer /*|| iwall->getType() == Wall::inner*/) {
                                        iwall->_pDz+=planeEndCapThikness;
                                        if (addElctronicBox) {
                                                iwall->_pRmax = travkVolRmax;
                                                iwall->_pRmin = travkVolRmax-iwall->getTotalThickness();
                                                outerWallInnerRadiusWithElec = iwall->_pRmin;
                                        }
                                }
                                else if (iwall->getType() == Wall::endcap) {
                                        iwall->_pRmin=innerWallOuterRadius;
                                        iwall->_pRmax=outerWallInnerRadius;
                                }
                        }
                }
                for ( walls_it=itracker->getWalls()->begin() ; walls_it != itracker->getWalls()->end(); walls_it++ ) {
                        VolumeInfo WallInfo;
                        boost::shared_ptr<Wall> iwall = walls_it->second;
                        WallInfo = buildWall(iwall.get(),itracker->endcapType());
                        WallInfo.logical->SetVisAttributes(visAttWall);
                        WallInfo.physical = new G4PVPlacement(iwall->getPos(),
                                        WallInfo.logical,                // its logical volume
                                        iwall->getName(),                // its name
                                        trackerInfo.logical,             // its mother  volume
                                        false,                           // no boolean operations
                                        0,                               // copy number
                                        checkOverlap);
                        //if (iwall->getType() == Wall::outer) outerWallInnerRadius = iwall->getRmin();
                }
                if (itracker->endcapType() == ITracker::Plane ) {
                        char tShapeName[50], tVolName[50];
                        double oringZEnd = oringCenter+oringWidth;
                        double extraInnWallLength = itracker->maxEndCapDim() - oringZEnd;
                        extraInnWallLength *= 0.5;
                        G4VSolid* extraInnerWallShape = new G4Tubs("extraInnerWall", ((G4Tubs*)trackerInfo.solid)->GetRMin(),innerWallOuterRadius,extraInnWallLength,0.0,360.0*CLHEP::degree);
                        G4LogicalVolume*   extraInnerWallVol = new G4LogicalVolume(extraInnerWallShape,Vacuum,"extraInnerWallVol",0,0,0);

                        for ( walls_it=itracker->getWalls()->begin() ; walls_it != itracker->getWalls()->end(); walls_it++ ) {
                                boost::shared_ptr<Wall> iwall = walls_it->second;
                                if (iwall->getType() == Wall::endcap) {
                                        if ( itracker->endcapType() == ITracker::Plane ) {
                                                Wall *wall = iwall.get();
                                                int nSub=wall->getNShells();
                                                int iSub=0;
                                                double iExtWallShellRmax = innerWallOuterRadius, iExtWallShellRmin;
                                                for (int ishell=nSub-config.getInt("itracker.structShellsAreLast",0); ishell<nSub; ishell++){

                                                        sprintf(tShapeName,"extraInnerWall_sub%i",iSub);
                                                        sprintf(tVolName,"extraInnerWallVol_sub%i",iSub);
                                                        iExtWallShellRmin = iExtWallShellRmax - wall->getThicknesses()->at(ishell);
                                                        G4Material *extInnWallShellMat = findMaterialOrThrow(wall->getMaterialsName()->at(ishell).c_str());

                                                        G4VSolid* extraInnerWallShellShape = new G4Tubs(tShapeName, iExtWallShellRmin,iExtWallShellRmax,extraInnWallLength,0.0,360.0*CLHEP::degree);
                                                        G4LogicalVolume*   extraInnerWallShellVol = new G4LogicalVolume(extraInnerWallShellShape,extInnWallShellMat,tVolName,0,0,0);
                                                        G4VPhysicalVolume *extraInnerWallShellPos = new G4PVPlacement(0,
                                                                        G4ThreeVector(0,0,0),
                                                                        extraInnerWallShellVol,   // its logical volume
                                                                        tVolName,                 // its name
                                                                        extraInnerWallVol,        // its mother  volume
                                                                        false,                    // no boolean operations
                                                                        0,                        // copy number
                                                                        detailedCheck);
                                                        extraInnerWallShellPos->GetCopyNo(); //just to remove the warning during compiling

                                                        iExtWallShellRmax = iExtWallShellRmin;
                                                        ++iSub;
                                                }
                                        }
                                }
                        }

                        double extraInnWallCenter = oringZEnd+extraInnWallLength;
                        G4VPhysicalVolume *extraInnerWallPos_R = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,extraInnWallCenter),
                                        extraInnerWallVol,        // its logical volume
                                        "extraInnerWallVol_R",      // its name
                                        trackerInfo.logical,      // its mother  volume
                                        false,                    // no boolean operations
                                        0,                        // copy number
                                        checkOverlap);
                        extraInnerWallPos_R->GetCopyNo(); //just to remove the warning during compiling
                        G4VPhysicalVolume *extraInnerWallPos_L = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,-extraInnWallCenter),
                                        extraInnerWallVol,        // its logical volume
                                        "extraInnerWallVol_L",      // its name
                                        trackerInfo.logical,      // its mother  volume
                                        false,                    // no boolean operations
                                        1,                        // copy number
                                        checkOverlap);
                        extraInnerWallPos_L->GetCopyNo(); //just to remove the warning during compiling

                }

                if ( config.getBool("itracker.detailedWireSupport",false) && addElctronicBox ) {
                        double  elctContLength = itracker->maxEndCapDim()-innerWallHaltLength;
                        double  elctContWallThick = config.getDouble("itracker.elctContWallThick");
                        G4Material *elctContWallMat = findMaterialOrThrow( config.getString("itracker.elctContWallMat") );
                        G4Material *elctContFillMat = findMaterialOrThrow( config.getString("itracker.elctContFillMat") );

                        elctContLength*=0.5;
                        G4VSolid* elctContShape = new G4Tubs("elctCont", outerWallInnerRadius,outerWallInnerRadiusWithElec,elctContLength,0.0,360.0*CLHEP::degree);
                        G4LogicalVolume*   elctContVol = new G4LogicalVolume(elctContShape,elctContFillMat,"elctContVol",0,0,0);

                        elctContWallThick*=0.5;
                        G4VSolid* elctContWallUpStShape = new G4Tubs("elctContWallUpSt", outerWallInnerRadius,outerWallInnerRadiusWithElec,elctContWallThick,0.0,360.0*CLHEP::degree);
                        G4LogicalVolume*   elctContWallUpStVol = new G4LogicalVolume(elctContWallUpStShape,elctContWallMat,"elctContWallUpStVol",0,0,0);
                        G4VPhysicalVolume *elctContWallUpStPos = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,-elctContLength+elctContWallThick),
                                        elctContWallUpStVol,     // its logical volume
                                        "elctContWallUpStVol",   // its name
                                        elctContVol,             // its mother  volume
                                        false,                   // no boolean operations
                                        0,                       // copy number
                                        detailedCheck);
                        elctContWallUpStPos->GetCopyNo(); //just to remove the warning during compiling

                        G4VSolid* elctContWallDwnStShape = new G4Tubs("elctContWallDwnSt", outerWallInnerRadius,outerWallInnerRadiusWithElec,elctContWallThick,0.0,360.0*CLHEP::degree);
                        G4LogicalVolume*   elctContWallDwnStVol = new G4LogicalVolume(elctContWallDwnStShape,elctContWallMat,"elctContWallDwnStVol",0,0,0);
                        G4VPhysicalVolume *elctContWallDwnStPos = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,elctContLength-elctContWallThick),
                                        elctContWallDwnStVol,    // its logical volume
                                        "elctContWallDwnStVol",  // its name
                                        elctContVol,             // its mother  volume
                                        false,                   // no boolean operations
                                        0,                       // copy number
                                        detailedCheck);
                        elctContWallDwnStPos->GetCopyNo(); //just to remove the warning during compiling

                        double elctContWallBottLeng = elctContLength-2.0*elctContWallThick;
                        double elctContWallBottRmax = outerWallInnerRadius+2.0*elctContWallThick;
                        G4VSolid* elctContWallBottShape = new G4Tubs("elctContWallBott", outerWallInnerRadius,elctContWallBottRmax,elctContWallBottLeng,0.0,360.0*CLHEP::degree);
                        G4LogicalVolume*   elctContWallBottVol = new G4LogicalVolume(elctContWallBottShape,elctContWallMat,"elctContWallBottVol",0,0,0);
                        G4VPhysicalVolume *elctContWallBottPos = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,0),
                                        elctContWallBottVol,     // its logical volume
                                        "elctContWallBottVol",   // its name
                                        elctContVol,             // its mother  volume
                                        false,                   // no boolean operations
                                        0,                       // copy number
                                        detailedCheck);
                        elctContWallBottPos->GetCopyNo(); //just to remove the warning during compiling

                        double elctContHeight = outerWallInnerRadiusWithElec-elctContWallBottRmax;

                        double  electCardLength = config.getDouble("itracker.electCardLength");
                        double  electCardThick  = config.getDouble("itracker.electCardThick");
                        double  electCardHeight = config.getDouble("itracker.electCardHeight");
                        G4Material *electCardMaterial = findMaterialOrThrow( config.getString("itracker.electCardMaterial") );
                        electCardLength*=0.5;

                        if ( electCardLength>elctContWallBottLeng || electCardHeight> elctContHeight ) { throw cet::exception("GEOM") <<"The electronic boards exceed their mother volume dim\n"; }

                        double electCardAngle = electCardThick/elctContWallBottRmax;
                        G4VSolid* electCardShape = new G4Tubs("electCard",elctContWallBottRmax,elctContWallBottRmax+electCardHeight,electCardLength,-0.5*electCardAngle,electCardAngle);
                        G4LogicalVolume*   electCardVol = new G4LogicalVolume(electCardShape,electCardMaterial,"electCardVol",0,0,0);
                        double electCardsNumber = itracker->nSuperLayers()*itracker->nRing()*itracker->nSWire()/config.getInt("itracker.electCardChanPerBord");
                        double electCardAngleStep = CLHEP::twopi/electCardsNumber;
                        for (int iElCard=0; iElCard<electCardsNumber; ++iElCard) {
                                G4VPhysicalVolume *electCardPos = new G4PVPlacement(
                                                HepGeom::Transform3D( HepGeom::RotateZ3D( ((double)iElCard)*electCardAngleStep ) ),
                                                electCardVol,     // its logical volume
                                                "electCardVol",   // its name
                                                elctContVol,      // its mother  volume
                                                false,            // no boolean operations
                                                iElCard,          // copy number
                                                detailedCheck);
                                electCardPos->GetCopyNo(); //just to remove the warning during compiling
                        }



                        G4VPhysicalVolume *elctContPos_R = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,itracker->maxEndCapDim()-elctContLength),
                                        elctContVol,             // its logical volume
                                        "elctContVol",           // its name
                                        trackerInfo.logical,     // its mother  volume
                                        false,                   // no boolean operations
                                        0,                       // copy number
                                        detailedCheck);
                        elctContPos_R->GetCopyNo(); //just to remove the warning during compiling

                        G4VPhysicalVolume *elctContPos_L = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,-(itracker->maxEndCapDim()-elctContLength)),
                                        elctContVol,             // its logical volume
                                        "elctContVol",           // its name
                                        trackerInfo.logical,     // its mother  volume
                                        false,                   // no boolean operations
                                        1,                       // copy number
                                        detailedCheck);
                        elctContPos_L->GetCopyNo(); //just to remove the warning during compiling


                        G4Material *leftFillGasFillMat = findMaterialOrThrow( config.getString("itracker.fillMaterial") );
                        G4VSolid* leftFillGasShape = new G4Tubs("leftFillGas", outerWallInnerRadius,outerWallInnerRadiusWithElec,innerWallHaltLength,0.0,360.0*CLHEP::degree);
                        G4LogicalVolume*   leftFillGasVol = new G4LogicalVolume(leftFillGasShape,leftFillGasFillMat,"leftFillGasVol",0,0,0);
                        G4VPhysicalVolume *leftFillGasPos = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,0),
                                        leftFillGasVol,           // its logical volume
                                        "leftFillGasVol",         // its name
                                        trackerInfo.logical,      // its mother  volume
                                        false,                    // no boolean operations
                                        0,                        // copy number
                                        detailedCheck);
                        leftFillGasPos->GetCopyNo(); //just to remove the warning during compiling
                }

                bool activeWireSD = config.getBool("itracker.ActiveWiresSD",false);
                for (int iSl = 0; iSl < itracker->nSuperLayers(); iSl++){

                        SuperLayer *SLayer = itracker->getSuperLayer(iSl);

                        for (int iLy=0; iLy < SLayer->nLayers(); iLy++ ){
                                //                   for (int iLy=0; iLy < SLayers[iSl].nLayers(); iLy++ ){

                                VolumeInfo LayerInfo;
                                boost::shared_ptr<ITLayer> ily = SLayer->getLayer(iLy);
                                //                        const ITLayer *ily = SLayers[iSl].getLayer(iLy);
                                superlayer = ily->Id().getSuperLayer();
                                iring = ily->Id().getLayer();
                                if (ily->getLayerType() == ITLayer::wire) {
                                        sprintf(shape,"wS%dR%d",superlayer,iring);
                                        sprintf(vol,"wvolS%02dR%02d",superlayer,iring);
                                } else if (ily->getLayerType() == ITLayer::gas) {
                                        sprintf(shape,"gS%dR%d",superlayer,iring);
                                        sprintf(vol,"gvolS%02dR%02d",superlayer,iring);
                                } else {
                                        sprintf(shape,"S%dR%d",superlayer,iring);
                                        sprintf(vol,"volS%02dR%02d",superlayer,iring);
                                }

                                //                     cout<<ily->Id()<<" IR "<<ily->getDetail()->centerInnerRadiusRing()<<" OR "<<
                                //                                     ily->getDetail()->centerOuterRadiusRing()<<" SI "<<ily->getDetail()->stereoAngleInnerRing()<<" SO "<<
                                //                                     ily->getDetail()->stereoAngleOuterRing()<<" HL "<<ily->getDetail()->halfLength()<<endl;

                                //cout<<ily->getDetail()->centerOuterRadiusRing()<<" "<<sqrt( square(ily->getDetail()->centerOuterRadiusRing()) +
                                //                                                square(ily->getDetail()->halfLength()*tan(ily->getDetail()->stereoAngleOuterRing())) ) <<" "<< outerWallInnerRadius<<endl;
                                //if ( sqrt( square(ily->getDetail()->centerOuterRadiusRing()) +
                                //                square(ily->getDetail()->halfLength()*tan(ily->getDetail()->stereoAngleOuterRing())) ) > outerWallInnerRadius )
                                //        throw cet::exception("GEOM") <<"The ITracker layer "<<ily->Id()<<" doesn't fit inside the ITracker outer wall\n";

                                LayerInfo.solid = new G4Hype(shape,ily->getDetail()->centerInnerRadiusRing(),
                                                ily->getDetail()->centerOuterRadiusRing(),ily->getDetail()->stereoAngleInnerRing(),
                                                ily->getDetail()->stereoAngleOuterRing(),ily->getDetail()->halfLength());
                                G4Material* GasMix = findMaterialOrThrow( ily->getDetail()->materialName() );
                                LayerInfo.logical = new G4LogicalVolume(LayerInfo.solid,GasMix,vol,0,0,0);
                                if (ily->voxelizationFactor()==0.0) {
                                        LayerInfo.logical->SetOptimisation(false);
                                }
                                else{
                                        LayerInfo.logical->SetSmartless(ily->voxelizationFactor());
                                }
                                LayerInfo.logical->SetVisAttributes(visAttGas);

                                LayerInfo.physical = new G4PVPlacement(0,               // no rotation
                                                positionTracker,         // at (x,y,z)
                                                LayerInfo.logical,       // its logical volume
                                                vol,                     // its name
                                                trackerInfo.logical,     // its mother  volume
                                                false,                   // no boolean operations
                                                0,                       // copy number
                                                checkOverlap);

                                if (ily->getLayerType() != ITLayer::undefined) {
//                                        ITGasLayerSD* glSD;
//                                        if ( itracker->geomType()==ITracker::Hexagonal ) glSD = new ITGasLayerSD( vol );//ITGasLayerSD_v2( vol );
//                                        else if ( itracker->geomType()==ITracker::Square ) glSD = new ITGasLayerSD( vol );//ITGasLayerSD_v3( vol );
//                                        SDman->AddNewDetector( glSD );
//                                        LayerInfo.logical->SetSensitiveDetector( glSD );
                                  G4VSensitiveDetector *sd = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::TrackerGas());
                                  if(sd) LayerInfo.logical->SetSensitiveDetector(sd);

                                }

                                VolumeInfo tmpFieldWireInfo;
                                for ( int iFw=0; iFw < ily->nFieldWires(); iFw++){
                                        VolumeInfo FieldWireInfo;
                                        boost::shared_ptr<Wire> iwire = ily->getFWire(iFw);
                                        boost::shared_ptr<WireDetail> wdet = iwire->getDetail();
                                        if (iFw==0) {
                                                sprintf(vol_name_FD,"tubeFD_%d_%d",superlayer,iring);
                                                tmpFieldWireInfo = buildWire(wdet->outerRadius(),wdet->halfLength(),shape_name_FD,vol_name_FD,wdet->materialNames(),wdet->shellsThicknesses(),activeWireSD);
                                                tmpFieldWireInfo.logical->SetVisAttributes(visAttFw);
                                        }
                                        FieldWireInfo.solid = tmpFieldWireInfo.solid;
                                        FieldWireInfo.logical = tmpFieldWireInfo.logical;
                                        sprintf(wire_name,"%s_%i",vol_name_FD,iwire->Id().getUId());
                                        FieldWireInfo.name = wire_name;
                                        FieldWireInfo.physical = new G4PVPlacement(iwire->get3DTransfrom(),
                                                        FieldWireInfo.logical,         // its logical volume
                                                        wire_name,                     // its name
                                                        LayerInfo.logical,             // its mother  volume
                                                        false,                         // no boolean operations
                                                        iwire->Id().getUId(),          // copy number
                                                        detailedCheck);
                                }

                                VolumeInfo tmpSenseWireInfo;
                                for ( int iSw=0; iSw < ily->nCells(); iSw++){
                                        VolumeInfo SenseWireInfo;
                                        boost::shared_ptr<Wire> iwire = ily->getCell(iSw)->getWire();
                                        boost::shared_ptr<WireDetail> wdet = iwire->getDetail();
                                        if (iSw==0) {
                                                sprintf(vol_name_SD,"tubeSD_%d_%d",superlayer,iring);
                                                tmpSenseWireInfo = buildWire(wdet->outerRadius(),wdet->halfLength(),shape_name_SD,vol_name_SD,wdet->materialNames(),wdet->shellsThicknesses(),activeWireSD,true);
                                                tmpSenseWireInfo.logical->SetVisAttributes(visAttSw);
                                        }
                                        SenseWireInfo.solid = tmpSenseWireInfo.solid;
                                        SenseWireInfo.logical = tmpSenseWireInfo.logical;
                                        sprintf(wire_name,"%s_%i",vol_name_SD,iwire->Id().getUId());
                                        SenseWireInfo.name = wire_name;
                                        SenseWireInfo.physical = new G4PVPlacement(iwire->get3DTransfrom(),
                                                        SenseWireInfo.logical,         // its logical volume
                                                        wire_name,                     // its name
                                                        LayerInfo.logical,             // its mother  volume
                                                        false,                         // no boolean operations
                                                        iwire->Id().getUId(),          // copy number
                                                        detailedCheck);
                                }

                        }

                }




                trackerInfo.physical =  new G4PVPlacement( 0,
                                trackerOffset,
                                trackerInfo.logical,
                                trackerName,
                                mother,
                                0,
                                0,
                                checkOverlap);

                if ( checkOverlap ) { cout<<"IT Overlap Checking "<<trackerInfo.physical->CheckOverlaps(100000,0.0001,true)<<endl; }
        }

        return trackerInfo;

}

VolumeInfo ITrackerBuilder::buildWire(float radius, float length, char *shapeName, char *volName, const std::vector<std::string> &materialName, const std::vector<double> &thicknesses, bool activeWireSD, bool isSense){

        VolumeInfo wire;
        wire.solid = new G4Tubs(shapeName,0.0,radius,length,0.0,360.0*CLHEP::degree);
        int nSub=materialName.size();
        if (nSub==1){
                wire.logical = new G4LogicalVolume(wire.solid,findMaterialOrThrow( materialName.at(0).c_str() ),volName,0,0,0);
        }
        else {
                wire.logical = new G4LogicalVolume(wire.solid,findMaterialOrThrow( "WAGVacuum" ),volName,0,0,0);
                char tShapeName[50], tVolName[50];
                double oldRadius = 0.0;
                double iRadius = 0.0;

                for (int ishell=0; ishell<nSub; ishell++){
                        sprintf(tShapeName,"%s_sub%i",shapeName,ishell);
                        sprintf(tVolName,"%s_sub%i",volName,ishell);
                        iRadius+=thicknesses.at(ishell);
                        G4Tubs *tswire = new G4Tubs(tShapeName,oldRadius,iRadius,length,0.0,360.0*CLHEP::degree);
                        //          cout<<tShapeName<<" "<<oldRadius<<" "<<iRadius<<" "<<length<<" "<<materialName.at(ishell)<<endl;
                        oldRadius=iRadius;

                        G4LogicalVolume *tlogicWire = new G4LogicalVolume(tswire,findMaterialOrThrow(materialName.at(ishell).c_str()),tVolName,0,0,0);
                        if (activeWireSD) {
                                if (isSense) {
                                        tlogicWire->SetSensitiveDetector( G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::TrackerSWires()) );
                                } else {
                                        tlogicWire->SetSensitiveDetector( G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::ITrackerFWires()) );
                                }
                        }

                        G4VPhysicalVolume *tphysWire = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,0),
                                        tlogicWire,       // its logical volume
                                        tVolName,         // its name
                                        wire.logical,     // its mother  volume
                                        false,            // no boolean operations
                                        0);               // copy number
                        tphysWire->GetCopyNo(); //just to remove the warning during compiling
                }
        }
        return wire;
}

VolumeInfo ITrackerBuilder::buildWall(Wall *wall, ITracker::EnCapType endcapType){

        art::ServiceHandle<GeometryService> geom;
        SimpleConfig const& config  = geom->config();
        bool hasDetailedSupport = config.getBool("itracker.detailedWireSupport",false);

        int skipSub(-1);

        VolumeInfo wallInfo;
        char volName[50], shapeName[50];
        sprintf(volName,"vol_%s",wall->getName().c_str());
        sprintf(shapeName,"shape_%s",wall->getName().c_str());
        if (wall->getType()==Wall::endcap && endcapType == ITracker::Spherical ) {
                wallInfo.solid = new G4Sphere( shapeName,wall->getRmin(),wall->getRmax(),wall->getSPhi(),wall->getDPhi(),wall->getSTheta(),wall->getDTheta() );
        } else {
                double rMax(wall->getRmax());
                if (wall->getType()==Wall::inner) {
                        int nSub=wall->getNShells();
                        for (int ishell=0; ishell<nSub; ishell++){
                                G4String iWallShellMat = wall->getMaterialsName()->at(ishell);
                                if ( iWallShellMat.contains("ITGas") ) {
                                        rMax-=wall->getThicknesses()->at(ishell);
                                        skipSub=ishell;
                                }
                        }
                }
                wallInfo.solid = new G4Tubs( shapeName,wall->getRmin(),rMax,wall->getDz(),wall->getSPhi(),wall->getDPhi() );
        }
        int nSub=wall->getNShells();
        if (nSub==1){
                wallInfo.logical = new G4LogicalVolume(wallInfo.solid,findMaterialOrThrow( wall->getMaterialsName()->at(0).c_str() ),volName,0,0,0);
        }
        else {
                wallInfo.logical = new G4LogicalVolume(wallInfo.solid,findMaterialOrThrow( "WAGVacuum" ),volName,0,0,0);
                char tShapeName[50], tVolName[50];
                double oldRadius = wall->getRmin();
                double iRadius = oldRadius;
                double iZpos=-wall->getDz();
                double iHalfThickness=0.0;
                double spdWebBaseExcess=0.0;
                for (int ishell=0; ishell<nSub; ishell++){
                        if (ishell==skipSub) { continue; }
                        sprintf(tShapeName,"%s_sub%i",shapeName,ishell);
                        sprintf(tVolName,"%s_sub%i",volName,ishell);
                        iRadius+=wall->getThicknesses()->at(ishell);
                        G4VSolid *tswall=0x0;
                        if (wall->getType()==Wall::endcap ){
                                if ( endcapType == ITracker::Spherical ) {
                                        tswall = new G4Sphere( tShapeName,oldRadius,iRadius,wall->getSPhi(),wall->getDPhi(),wall->getSTheta(),wall->getDTheta() );
                                        //          cout<<tShapeName<<" "<<oldRadius<<" "<<iRadius<<" "<<length<<" "<<wall->getMaterialsName()->at(ishell)<<endl;
                                        oldRadius=iRadius;
                                        iZpos=0.0;
                                } else if ( endcapType == ITracker::Plane ) {
                                        iHalfThickness=wall->getThicknesses()->at(ishell)*0.5;
                                        iZpos+=iHalfThickness;
                                        tswall = new G4Tubs(tShapeName,wall->getRmin(),wall->getRmax(),iHalfThickness,wall->getSPhi(),wall->getDPhi() );
                                }
                        } else {
                                tswall = new G4Tubs(tShapeName,oldRadius,iRadius,wall->getDz(),wall->getSPhi(),wall->getDPhi() );
                                oldRadius=iRadius;
                                iZpos=0.0;
                        }


                        G4LogicalVolume *tlogicwall = new G4LogicalVolume(tswall,findMaterialOrThrow(wall->getMaterialsName()->at(ishell).c_str()),tVolName,0,0,0);
                        if ( hasDetailedSupport && wall->getType()==Wall::endcap && endcapType == ITracker::Plane ) {
                                switch (ishell) {
                                        case 0:
                                                spdWebBaseExcess=constructSpiderWeb(tlogicwall,config);
                                                break;
                                        case 1:
                                                constructWireAnchoring(tlogicwall,config,spdWebBaseExcess);
                                                break;
                                        case 2:
                                                constructSignalCables(tlogicwall,config);
                                                constructHvCables(tlogicwall,config);
                                                break;
                                        default:
                                                break;
                                }
                        }
                        G4VPhysicalVolume *tphyswall = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,iZpos),
                                        tlogicwall,       // its logical volume
                                        tVolName,         // its name
                                        wallInfo.logical, // its mother  volume
                                        false,            // no boolean operations
                                        0,                // copy number
                                        detailedCheck);
                        tphyswall->GetCopyNo(); //just to remove the warning during compiling
                        iZpos+=iHalfThickness;
                }

        }
        return wallInfo;
}

double ITrackerBuilder::constructSpiderWeb(G4LogicalVolume* localMother, SimpleConfig const& config) {
        GeomHandle<ITracker> itracker;
        double minR, maxR;
        char tShapeName[50], tVolName[50];
        double motherDz = ((G4Tubs *)localMother->GetSolid())->GetZHalfLength();

        minR = ((G4Tubs *)localMother->GetSolid())->GetInnerRadius();
        maxR = ((G4Tubs *)localMother->GetSolid())->GetOuterRadius();

        G4Material* spdWebMat =  findMaterialOrThrow( config.getString("itracker.spdWebSpokeMaterial") );

        int spdWebSpokesNumber = config.getInt("itracker.spdWebSpokesNumber");
        vector<double> spdWebSpokeFacePntsX, spdWebSpokeFacePntsY;
        config.getVectorDouble("itracker.spdWebSpokeFacePntsX",spdWebSpokeFacePntsX);
        config.getVectorDouble("itracker.spdWebSpokeFacePntsY",spdWebSpokeFacePntsY);
        if (spdWebSpokeFacePntsX.size()!=spdWebSpokeFacePntsY.size()) { throw cet::exception("GEOM") <<"The number of points X is different from that of Y for the vertexes coordinates of the spoke shape\n"; }
        double spdWebBaseWidth = config.getDouble("itracker.spdWebBaseWidth");
        double spdWebBaseThickness = config.getDouble("itracker.spdWebBaseThickness");

        double effSpdWebBaseWidth = 0.5*spdWebBaseWidth;
        double spdWebBaseWidthExcess = 0.0;
        if (effSpdWebBaseWidth>motherDz) { effSpdWebBaseWidth=motherDz; spdWebBaseWidthExcess=spdWebBaseWidth-2.0*motherDz; }

        sprintf(tShapeName,"spdWebBase");
        sprintf(tVolName,"spdWebBaseVol");
        G4VSolid* spdWebBaseShape = new G4Tubs(tShapeName, minR,minR+spdWebBaseThickness,effSpdWebBaseWidth,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   spdWebBaseVol = new G4LogicalVolume(spdWebBaseShape,spdWebMat,tVolName,0,0,0);
        G4VPhysicalVolume *spdWebBasePos = new G4PVPlacement(0,
                                               G4ThreeVector(0,0,0),
                                               spdWebBaseVol,    // its logical volume
                                               tVolName,         // its name
                                               localMother,      // its mother  volume
                                               false,            // no boolean operations
                                               0,                // copy number
                                               detailedCheck);
        spdWebBasePos->GetCopyNo(); //just to remove the warning during compiling

        sprintf(tShapeName,"spdWebSpoke");
        sprintf(tVolName,"spdWebSpokeVol");

        double spdWebSpokeHeight = maxR-minR-spdWebBaseThickness;
        spdWebSpokeHeight -= (maxR - std::sqrt( maxR*maxR - 100.0 )); //safety factor assuming a max front spoke size of 20 mm
        std::vector<G4TwoVector> spdWebSpokeFacePnts;
        for (unsigned int ipnts=0; ipnts<spdWebSpokeFacePntsX.size(); ++ipnts) {
                spdWebSpokeFacePnts.push_back( G4TwoVector( spdWebSpokeFacePntsX.at(ipnts),spdWebSpokeFacePntsY.at(ipnts) ) );
        }
        G4VSolid* spdWebSpokeShape;
        if (config.getBool("itracker.useSimplefiedSpoke",false)) {
                double simplefiedSpokeWidth = 0.5*config.getDouble("itracker.simplefiedSpokeWidth");
                spdWebSpokeShape = new G4Box(tShapeName,simplefiedSpokeWidth,simplefiedSpokeWidth,0.5*spdWebSpokeHeight);
        } else {
                spdWebSpokeShape = new G4ExtrudedSolid(tShapeName, spdWebSpokeFacePnts, 0.5*spdWebSpokeHeight,G4TwoVector(0.0,0.0),1.0,G4TwoVector(0.0,0.0),1.0);
        }
        G4LogicalVolume*   spdWebSpokeVol = new G4LogicalVolume(spdWebSpokeShape,spdWebMat,tVolName,0,0,0);
        HepGeom::Transform3D baseSpokePos ( HepGeom::TranslateX3D(0.5*spdWebSpokeHeight+minR+spdWebBaseThickness) /** HepGeom::RotateZ3D(90.0*CLHEP::degree)*/ * HepGeom::RotateY3D(90.0*CLHEP::degree) );
        double spokeAngleStep = CLHEP::twopi/((double)spdWebSpokesNumber);
        for (int iSpokes=0; iSpokes<spdWebSpokesNumber; ++iSpokes) {
                G4VPhysicalVolume *spdWebSpokePos = new G4PVPlacement(
                                HepGeom::Transform3D( HepGeom::RotateZ3D( ((double)iSpokes)*spokeAngleStep ) * baseSpokePos ),
                                spdWebSpokeVol,   // its logical volume
                                tVolName,         // its name
                                localMother,      // its mother  volume
                                false,            // no boolean operations
                                iSpokes,          // copy number
                                detailedCheck);
                spdWebSpokePos->GetCopyNo(); //just to remove the warning during compiling
        }

        return spdWebBaseWidthExcess;
}

void ITrackerBuilder::constructWireAnchoring(G4LogicalVolume* localMother, SimpleConfig const& config, double spdWebBaseExcess) {
        GeomHandle<ITracker> itracker;

        bool isDownStream = localMother->GetName().contains("_R");
        string endCapSide;
        if (isDownStream) {
                endCapSide = "_R";
        } else {
                endCapSide = "_L";
        }

        double minR, maxR, swR, tmpZincr;
        maxR=0.0;
        char tShapeName[50], tVolName[50];
        char ancrVolName[50], spacerVolName[50], componContVolName[50];
        double motherDz = ((G4Tubs *)localMother->GetSolid())->GetZHalfLength();
        G4Material* motherMat = localMother->GetMaterial();
        int superlayer =0,iring=0;
        //cout<<"vol name "<< localMother->GetName()<<" half z "<<motherDz<<" with matrial "<<*(localMother->GetMaterial())<<endl;

        //parameters for Field wire Boards
        double FwBoardWidth, FwBoardThick;
        FwBoardWidth = config.getDouble("itracker.fwBoardWidth");
        FwBoardThick = config.getDouble("itracker.fwBoardThickness");
        G4Material* FwBoardMat =  findMaterialOrThrow( config.getString("itracker.fwBoardMaterial") );
        FwBoardWidth*=0.5;
        if (FwBoardWidth>motherDz) { throw cet::exception("GEOM") <<"The wire anchoring Fw board exceeds its mother volume dim\n"; }
        G4ThreeVector fwRelPos(0,0,-motherDz+FwBoardWidth);

        //parameters for Sense wire Boards
        double SwBoardWidth, SwBoardThick;
        SwBoardWidth = config.getDouble("itracker.swBoardWidth");
        SwBoardThick = config.getDouble("itracker.swBoardThickness");
        G4Material* SwBoardMat =  findMaterialOrThrow( config.getString("itracker.swBoardMaterial") );
        SwBoardWidth*=0.5;
        if (FwBoardWidth>motherDz) { throw cet::exception("GEOM") <<"The wire anchoring Sw board exceeds its mother volume dim\n"; }
        G4ThreeVector swRelPos(0,0,-motherDz+SwBoardWidth);

        //parameters for spacers
        double startingSpacerRmin, startingSpacerRmax, finishingSpacerRmin, finishingSpacerRmax;
        startingSpacerRmax=0.0;
        double spacerRmin, spacerRmax, spacerBaseWidth, spacerBaseThick, spacerCoreThick;
        double spacerCoreSurfFrac, spokeAngle;
        int    spacerNumCores, spacerNumSpokesPerCore;
        spacerBaseWidth = config.getDouble("itracker.spacerBaseWidth");
        spacerBaseThick = config.getDouble("itracker.spacerBaseThickness");
        spacerCoreThick = config.getDouble("itracker.spacerCoreThickness");
        spacerCoreSurfFrac = config.getDouble("itracker.spacerCoreSurfFrac");
        spacerNumCores = config.getInt("itracker.spacerNumCores");
        spacerNumSpokesPerCore = config.getInt("itracker.spacerNumSpokesPerCore");
        G4Material* spacerMat =  findMaterialOrThrow( config.getString("itracker.spacerMaterial") );
        spacerBaseWidth*=0.5;
        spacerCoreThick*=0.5;
        if (spacerBaseWidth>motherDz || spacerCoreThick>spacerBaseWidth) { throw cet::exception("GEOM") <<"The wire anchoring Spacer exceeds its mother volume dim\n"; }
        G4ThreeVector spacerRelPos(0,0,-motherDz+spacerBaseWidth);
        int spacerCoreNumSpokes = spacerNumSpokesPerCore*spacerNumCores;
        double spetSpokeRot = CLHEP::twopi/((double)spacerCoreNumSpokes);
        std::vector<HepGeom::RotateZ3D *> spacerSpokesRot;// = new HepGeom::RotateZ3D[spacerCoreNumSpokes];
        for (int iSpSpoke=0; iSpSpoke<spacerCoreNumSpokes; ++iSpSpoke) {
                spacerSpokesRot.push_back( new HepGeom::RotateZ3D( ((double)iSpSpoke)*spetSpokeRot ) );
                //spacerSpokesRot[iSpSpoke].set( CLHEP::HepRotationZ( ((double)iSpSpoke)*spetSpokeRot ) );
        }
        spokeAngle = spetSpokeRot*spacerCoreSurfFrac;


        double nCellPerLayer;
        double componContWidth = motherDz-spacerBaseWidth;
        double compomContHeight, componContRmin, componContRmax, compStepAngle;
        G4ThreeVector componContRelPos(0,0,motherDz-componContWidth);
        //parameters for  HV Capacitance
        G4Material* hvCapMat =  findMaterialOrThrow( config.getString("itracker.hvCapMaterial") );
        double hvCapLength = 0.5*config.getDouble("itracker.hvCapLength");
        double hvCapWidth  = config.getDouble("itracker.hvCapWidth");
        double hvCapHeight = config.getDouble("itracker.hvCapHeight");
        double hvCapAngle;

        if (hvCapLength>componContWidth) { throw cet::exception("GEOM") <<"The HV Capacitance exceeds its mother volume dim\n"; }
        compomContHeight = hvCapHeight;

        //parameters for  Termination Resistance
        G4Material* termResMat =  findMaterialOrThrow( config.getString("itracker.trmResMaterial") );
        double termResLength = 0.5*config.getDouble("itracker.trmResLength");
        double termResWidth  = config.getDouble("itracker.trmResWidth");
        double termResHeight = config.getDouble("itracker.trmResHeight");
        double termResAngle;

        if (termResLength>componContWidth) { throw cet::exception("GEOM") <<"The termination Resistance exceeds its mother volume dim\n"; }
        if (termResHeight>compomContHeight) { compomContHeight = termResHeight; }

        //parameters for  Hv Resistance
        G4Material* hvResMat =  findMaterialOrThrow( config.getString("itracker.hvResMaterial") );
        double hvResLength = 0.5*config.getDouble("itracker.hvResLength");
        double hvResWidth  = config.getDouble("itracker.hvResWidth");
        double hvResHeight = config.getDouble("itracker.hvResHeight");
        double hvResAngle;

        if (hvResLength>componContWidth) { throw cet::exception("GEOM") <<"The HV Resistance exceeds its mother volume dim\n"; }
        if (hvResHeight>compomContHeight) { compomContHeight = hvResHeight; }

        bool firstLayer=true;
        for (int iSl = 0; iSl < itracker->nSuperLayers(); iSl++){

                SuperLayer *SLayer = itracker->getSuperLayer(iSl);

                for (int iLy=0; iLy < SLayer->nLayers(); iLy++ ){

                        VolumeInfo LayerInfo;
                        boost::shared_ptr<ITLayer> ily = SLayer->getLayer(iLy);
                        superlayer = ily->Id().getSuperLayer();
                        iring = ily->Id().getLayer();
                        if (ily->nCells()>0) {
                                /*
                                        swR = ily->getDetail()->centerInnerRadiusRing();
                                        tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleInnerRing());
                                        swR = std::sqrt( swR*swR + tmpZincr*tmpZincr);

                                        boost::shared_ptr<ITLayer> ilyDwn = SLayer->getLayer(iLy-1);
                                        minR = ilyDwn->getDetail()->centerInnerRadiusRing();
                                        tmpZincr = ilyDwn->getDetail()->halfLength()*std::tan(ilyDwn->getDetail()->stereoAngleInnerRing());
                                        minR = std::sqrt( minR*minR + tmpZincr*tmpZincr);

                                        boost::shared_ptr<ITLayer> ilyUp  = SLayer->getLayer(iLy+1);
                                        maxR = ilyUp->getDetail()->centerOuterRadiusRing();
                                        tmpZincr = ilyUp->getDetail()->halfLength()*std::tan(ilyUp->getDetail()->stereoAngleOuterRing());
                                        maxR = std::sqrt( maxR*maxR + tmpZincr*tmpZincr);
                                 */

                                minR = ily->getDetail()->centerInnerRadiusRing();
                                tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleInnerRing());
                                minR = std::sqrt( minR*minR + tmpZincr*tmpZincr);
                                if(firstLayer) {
                                        startingSpacerRmax = minR;
                                        firstLayer = false;
                                }

                                maxR = ily->getDetail()->centerOuterRadiusRing();
                                tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleOuterRing());
                                maxR = std::sqrt( maxR*maxR + tmpZincr*tmpZincr);

                                swR = (maxR+minR)*0.5;

                                sprintf(tShapeName,"anchorS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(ancrVolName,"anchorVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* anchorShape = new G4Tubs(tShapeName, minR,maxR,motherDz,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   anchorVol = new G4LogicalVolume(anchorShape,motherMat,ancrVolName,0,0,0);


                                //------------- start wires boards -------------
                                sprintf(tShapeName,"brdFldDwnS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"brdFldDwnVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* brdFldDwnShape = new G4Tubs(tShapeName, minR,minR+FwBoardThick,FwBoardWidth,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   brdFldDwnVol = new G4LogicalVolume(brdFldDwnShape,FwBoardMat,tVolName,0,0,0);
                                G4VPhysicalVolume *brdFldDwnPos = new G4PVPlacement(0,
                                                fwRelPos,
                                                brdFldDwnVol,     // its logical volume
                                                tVolName,         // its name
                                                anchorVol,        // its mother  volume
                                                false,            // no boolean operations
                                                0,                // copy number
                                                detailedCheck);
                                brdFldDwnPos->GetCopyNo(); //just to remove the warning during compiling

                                sprintf(tShapeName,"brdFldUpS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"brdFldUpVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* brdFldUpShape = new G4Tubs(tShapeName, maxR-FwBoardThick,maxR,FwBoardWidth,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   brdFldUpVol = new G4LogicalVolume(brdFldUpShape,FwBoardMat,tVolName,0,0,0);
                                G4VPhysicalVolume *brdFldUpPos = new G4PVPlacement(0,
                                                fwRelPos,
                                                brdFldUpVol,      // its logical volume
                                                tVolName,         // its name
                                                anchorVol,        // its mother  volume
                                                false,            // no boolean operations
                                                0,                // copy number
                                                detailedCheck);
                                brdFldUpPos->GetCopyNo(); //just to remove the warning during compiling

                                sprintf(tShapeName,"brdSncS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"brdSncVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* brdSncShape = new G4Tubs(tShapeName, swR-SwBoardThick,swR,SwBoardWidth,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   brdSncVol = new G4LogicalVolume(brdSncShape,SwBoardMat,tVolName,0,0,0);
                                G4VPhysicalVolume *brdSncPos = new G4PVPlacement(0,
                                                swRelPos,
                                                brdSncVol,        // its logical volume
                                                tVolName,         // its name
                                                anchorVol,        // its mother  volume
                                                false,            // no boolean operations
                                                0,                // copy number
                                                detailedCheck);
                                brdSncPos->GetCopyNo(); //just to remove the warning during compiling
                                //------------- end wires boards -------------

                                //------------- start spacers -------------
                                //start down spacers
                                spacerRmin = minR+FwBoardThick;
                                spacerRmax = swR-SwBoardThick;

                                sprintf(tShapeName,"spacerDwnS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(spacerVolName,"spacerDwnVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* spacerDwnShape = new G4Tubs(tShapeName, spacerRmin,spacerRmax,spacerBaseWidth,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   spacerDwnVol = new G4LogicalVolume(spacerDwnShape,motherMat,spacerVolName,0,0,0);

                                sprintf(tShapeName,"spacerDwnBotS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"spacerDwnBotVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* spacerDwnBotShape = new G4Tubs(tShapeName, spacerRmin,spacerRmin+spacerBaseThick,spacerBaseWidth,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   spacerDwnBotVol = new G4LogicalVolume(spacerDwnBotShape,spacerMat,tVolName,0,0,0);
                                G4VPhysicalVolume *spacerDwnBotPos = new G4PVPlacement(0,
                                                G4ThreeVector(0,0,0),
                                                spacerDwnBotVol,  // its logical volume
                                                tVolName,         // its name
                                                spacerDwnVol,     // its mother  volume
                                                false,            // no boolean operations
                                                0,                // copy number
                                                detailedCheck);
                                spacerDwnBotPos->GetCopyNo(); //just to remove the warning during compiling

                                sprintf(tShapeName,"spacerDwnTopS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"spacerDwnTopVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* spacerDwnTopShape = new G4Tubs(tShapeName, spacerRmax-spacerBaseThick,spacerRmax,spacerBaseWidth,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   spacerDwnTopVol = new G4LogicalVolume(spacerDwnTopShape,spacerMat,tVolName,0,0,0);
                                G4VPhysicalVolume *spacerDwnTopPos = new G4PVPlacement(0,
                                                G4ThreeVector(0,0,0),
                                                spacerDwnTopVol,  // its logical volume
                                                tVolName,         // its name
                                                spacerDwnVol,     // its mother  volume
                                                false,            // no boolean operations
                                                0,                // copy number
                                                detailedCheck);
                                spacerDwnTopPos->GetCopyNo(); //just to remove the warning during compiling

                                sprintf(tShapeName,"spacerDwnSpokeS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"spacerDwnSpokeVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* spacerDwnSpokeShape = new G4Tubs(tShapeName, spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,spacerCoreThick,-0.5*spokeAngle,spokeAngle);
                                G4LogicalVolume*   spacerDwnSpokeVol = new G4LogicalVolume(spacerDwnSpokeShape,spacerMat,tVolName,0,0,0);
                                for (int iSpSpoke=0; iSpSpoke<spacerCoreNumSpokes; ++iSpSpoke){
                                        G4VPhysicalVolume *spacerDwnSpokePos = new G4PVPlacement(*spacerSpokesRot.at(iSpSpoke),
                                                        spacerDwnSpokeVol,  // its logical volume
                                                        tVolName,           // its name
                                                        spacerDwnVol,       // its mother  volume
                                                        false,              // no boolean operations
                                                        iSpSpoke,           // copy number
                                                        detailedCheck);
                                        spacerDwnSpokePos->GetCopyNo(); //just to remove the warning during compiling
                                }

                                G4VPhysicalVolume *spacerDwnPos = new G4PVPlacement(0,
                                                spacerRelPos,
                                                spacerDwnVol,     // its logical volume
                                                spacerVolName,    // its name
                                                anchorVol,        // its mother  volume
                                                false,            // no boolean operations
                                                0,                // copy number
                                                detailedCheck);
                                spacerDwnPos->GetCopyNo(); //just to remove the warning during compiling
                                //end down spacers

                                //start up spacers
                                spacerRmin = swR;
                                spacerRmax = maxR-FwBoardThick;

                                sprintf(tShapeName,"spacerUpS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(spacerVolName,"spacerUpVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* spacerUpShape = new G4Tubs(tShapeName, spacerRmin,spacerRmax,spacerBaseWidth,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   spacerUpVol = new G4LogicalVolume(spacerUpShape,motherMat,spacerVolName,0,0,0);

                                sprintf(tShapeName,"spacerUpBotS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"spacerUpBotVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* spacerUpBotShape = new G4Tubs(tShapeName, spacerRmin,spacerRmin+spacerBaseThick,spacerBaseWidth,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   spacerUpBotVol = new G4LogicalVolume(spacerUpBotShape,spacerMat,tVolName,0,0,0);
                                G4VPhysicalVolume *spacerUpBotPos = new G4PVPlacement(0,
                                                G4ThreeVector(0,0,0),
                                                spacerUpBotVol,   // its logical volume
                                                tVolName,         // its name
                                                spacerUpVol,      // its mother  volume
                                                false,            // no boolean operations
                                                0,                // copy number
                                                detailedCheck);
                                spacerUpBotPos->GetCopyNo(); //just to remove the warning during compiling

                                sprintf(tShapeName,"spacerUpTopS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"spacerUpTopVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* spacerUpTopShape = new G4Tubs(tShapeName, spacerRmax-spacerBaseThick,spacerRmax,spacerBaseWidth,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   spacerUpTopVol = new G4LogicalVolume(spacerUpTopShape,spacerMat,tVolName,0,0,0);
                                G4VPhysicalVolume *spacerUpTopPos = new G4PVPlacement(0,
                                                G4ThreeVector(0,0,0),
                                                spacerUpTopVol,   // its logical volume
                                                tVolName,         // its name
                                                spacerUpVol,      // its mother  volume
                                                false,            // no boolean operations
                                                0,                // copy number
                                                detailedCheck);
                                spacerUpTopPos->GetCopyNo(); //just to remove the warning during compiling

                                sprintf(tShapeName,"spacerUpSpokeS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"spacerUpSpokeVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* spacerUpSpokeShape = new G4Tubs(tShapeName, spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,spacerCoreThick,-0.5*spokeAngle,spokeAngle);
                                G4LogicalVolume*   spacerUpSpokeVol = new G4LogicalVolume(spacerUpSpokeShape,spacerMat,tVolName,0,0,0);
                                for (int iSpSpoke=0; iSpSpoke<spacerCoreNumSpokes; ++iSpSpoke){
                                        G4VPhysicalVolume *spacerUpSpokePos = new G4PVPlacement(*spacerSpokesRot.at(iSpSpoke),
                                                        spacerUpSpokeVol,   // its logical volume
                                                        tVolName,           // its name
                                                        spacerUpVol,        // its mother  volume
                                                        false,              // no boolean operations
                                                        iSpSpoke,           // copy number
                                                        detailedCheck);
                                        spacerUpSpokePos->GetCopyNo(); //just to remove the warning during compiling
                                }

                                G4VPhysicalVolume *spacerUpPos = new G4PVPlacement(0,
                                                spacerRelPos,
                                                spacerUpVol,      // its logical volume
                                                spacerVolName,    // its name
                                                anchorVol,        // its mother  volume
                                                false,            // no boolean operations
                                                0,                // copy number
                                                detailedCheck);
                                spacerUpPos->GetCopyNo(); //just to remove the warning during compiling
                                //end up spacers
                                //------------- end spacers -------------

                                nCellPerLayer = ily->nCells();
                                compStepAngle = CLHEP::twopi/nCellPerLayer;
                                //------------- start Term components -------------
                                componContRmin = swR;
                                componContRmax = componContRmin+compomContHeight;

                                sprintf(tShapeName,"componContTermS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(componContVolName,"componContTermVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                G4VSolid* componContTermShape = new G4Tubs(tShapeName, componContRmin,componContRmax,componContWidth,0.0,360.0*CLHEP::degree);
                                G4LogicalVolume*   componContTermVol = new G4LogicalVolume(componContTermShape,motherMat,componContVolName,0,0,0);

                                //start HV capacitance
                                sprintf(tShapeName,"hvCapS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"hvCapVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                hvCapAngle = hvCapWidth/componContRmin;
                                G4VSolid* hvCapShape = new G4Tubs(tShapeName, componContRmin,componContRmin+hvCapHeight,hvCapLength,-0.5*hvCapAngle,hvCapAngle);
                                G4LogicalVolume*   hvCapVol = new G4LogicalVolume(hvCapShape,hvCapMat,tVolName,0,0,0);
                                for (int iNcell=0; iNcell<nCellPerLayer; ++iNcell){
                                        G4VPhysicalVolume *hvCapPos = new G4PVPlacement( HepGeom::RotateZ3D( ((double)iNcell)*compStepAngle ),
                                                        hvCapVol,           // its logical volume
                                                        tVolName,           // its name
                                                        componContTermVol,  // its mother  volume
                                                        false,              // no boolean operations
                                                        iNcell,             // copy number
                                                        detailedCheck);
                                        hvCapPos->GetCopyNo(); //just to remove the warning during compiling
                                }

                                //end HV capacitance

                                //start Term Resistance
                                sprintf(tShapeName,"termResS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                sprintf(tVolName,"termResVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                termResAngle = termResWidth/componContRmin;
                                G4VSolid* termResShape = new G4Tubs(tShapeName, componContRmin,componContRmin+termResHeight,termResLength,-0.5*termResAngle,termResAngle);
                                G4LogicalVolume*   termResVol = new G4LogicalVolume(termResShape,termResMat,tVolName,0,0,0);
                                for (int iNcell=0; iNcell<nCellPerLayer; ++iNcell){
                                        G4VPhysicalVolume *termResPos = new G4PVPlacement( HepGeom::RotateZ3D( (((double)iNcell)+0.5)*compStepAngle ),
                                                        termResVol,         // its logical volume
                                                        tVolName,           // its name
                                                        componContTermVol,  // its mother  volume
                                                        false,              // no boolean operations
                                                        iNcell,             // copy number
                                                        detailedCheck);
                                        termResPos->GetCopyNo(); //just to remove the warning during compiling
                                }
                                //end Term Resistance

                                G4VPhysicalVolume *componContTermPos = new G4PVPlacement(0,
                                                componContRelPos,
                                                componContTermVol, // its logical volume
                                                componContVolName, // its name
                                                anchorVol,         // its mother  volume
                                                false,             // no boolean operations
                                                0,                 // copy number
                                                detailedCheck);
                                componContTermPos->GetCopyNo(); //just to remove the warning during compiling
                                //------------- end Term components -------------

                                //------------- start HV Resistance components -------------
                                if (isDownStream) {
                                        componContRmax = swR-SwBoardThick;
                                        componContRmin = componContRmax-compomContHeight;

                                        sprintf(tShapeName,"componContHvS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                        sprintf(componContVolName,"componContHvVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                        G4VSolid* componContTermShape = new G4Tubs(tShapeName, componContRmin,componContRmax,componContWidth,0.0,360.0*CLHEP::degree);
                                        G4LogicalVolume*   componContTermVol = new G4LogicalVolume(componContTermShape,motherMat,componContVolName,0,0,0);

                                        sprintf(tShapeName,"hvResS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                        sprintf(tVolName,"hvResVolS%dR%d%s",superlayer,iring,endCapSide.c_str());
                                        hvResAngle = hvResWidth/componContRmin;
                                        G4VSolid* hvResShape = new G4Tubs(tShapeName, componContRmax-hvResHeight,componContRmax,hvResLength,-0.5*hvResAngle,hvResAngle);
                                        G4LogicalVolume*   hvResVol = new G4LogicalVolume(hvResShape,hvResMat,tVolName,0,0,0);
                                        for (int iNcell=0; iNcell<nCellPerLayer; ++iNcell){
                                                G4VPhysicalVolume *hvResPos = new G4PVPlacement( HepGeom::RotateZ3D( (((double)iNcell)+0.5)*compStepAngle ),
                                                                hvResVol,           // its logical volume
                                                                tVolName,           // its name
                                                                componContTermVol,  // its mother  volume
                                                                false,              // no boolean operations
                                                                iNcell,             // copy number
                                                                detailedCheck);
                                                hvResPos->GetCopyNo(); //just to remove the warning during compiling
                                        }

                                        G4VPhysicalVolume *componContTermPos = new G4PVPlacement(0,
                                                        componContRelPos,
                                                        componContTermVol, // its logical volume
                                                        componContVolName, // its name
                                                        anchorVol,         // its mother  volume
                                                        false,             // no boolean operations
                                                        0,                 // copy number
                                                        detailedCheck);
                                        componContTermPos->GetCopyNo(); //just to remove the warning during compiling
                                }
                                //------------- end HV Resistance components -------------

                                G4VPhysicalVolume *anchorPos = new G4PVPlacement(0,
                                                G4ThreeVector(0,0,0),
                                                anchorVol,        // its logical volume
                                                ancrVolName,      // its name
                                                localMother,      // its mother  volume
                                                false,            // no boolean operations
                                                0,                // copy number
                                                detailedCheck);
                                anchorPos->GetCopyNo(); //just to remove the warning during compiling

                        }
                }

        }

        minR = ((G4Tubs *)localMother->GetSolid())->GetInnerRadius();
        startingSpacerRmin = minR;
        //remaining part of the base of the spider web
        if (spdWebBaseExcess>0.0) {
                G4Material* spdWebMat =  findMaterialOrThrow( config.getString("itracker.spdWebSpokeMaterial") );
                double spdWebBaseThickness = config.getDouble("itracker.spdWebBaseThickness");
                double effSpdWebBaseWidth = 0.5*spdWebBaseExcess;
                if (effSpdWebBaseWidth>motherDz) { throw cet::exception("GEOM") <<"The exceeding part of the Spider Web Base exceeds its mother volume dim\n"; }
                startingSpacerRmin = minR+spdWebBaseThickness;

                sprintf(tShapeName,"spdWebBase");
                sprintf(tVolName,"spdWebBaseVol");
                G4VSolid* spdWebBaseShape = new G4Tubs(tShapeName, minR, startingSpacerRmin, effSpdWebBaseWidth,0.0,360.0*CLHEP::degree);
                G4LogicalVolume*   spdWebBaseVol = new G4LogicalVolume(spdWebBaseShape,spdWebMat,tVolName,0,0,0);
                G4VPhysicalVolume *spdWebBasePos = new G4PVPlacement(0,
                                                       G4ThreeVector(0,0,-motherDz+effSpdWebBaseWidth),
                                                       spdWebBaseVol,    // its logical volume
                                                       tVolName,         // its name
                                                       localMother,      // its mother  volume
                                                       false,            // no boolean operations
                                                       0,                // copy number
                                                       detailedCheck);
                spdWebBasePos->GetCopyNo(); //just to remove the warning during compiling
        }

        //------------- start starting spacer and guard field wires -------------
        sprintf(tShapeName,"anchorS%dR%d",0,-1);
        sprintf(ancrVolName,"anchorVolS%dR%d",0,-1);
        G4VSolid* anchorShape = new G4Tubs(tShapeName, startingSpacerRmin,startingSpacerRmax,motherDz,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   anchorVol = new G4LogicalVolume(anchorShape,motherMat,ancrVolName,0,0,0);

        sprintf(tShapeName,"brdFldUpS%dR%d",0,-1);
        sprintf(tVolName,"brdFldUpVolS%dR%d",0,-1);
        G4VSolid* brdFldUpShape = new G4Tubs(tShapeName, startingSpacerRmax-FwBoardThick,startingSpacerRmax,FwBoardWidth,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   brdFldUpVol = new G4LogicalVolume(brdFldUpShape,FwBoardMat,tVolName,0,0,0);
        G4VPhysicalVolume *brdFldUpPos = new G4PVPlacement(0,
                        fwRelPos,
                        brdFldUpVol,      // its logical volume
                        tVolName,         // its name
                        anchorVol,        // its mother  volume
                        false,            // no boolean operations
                        0,                // copy number
                        detailedCheck);
        brdFldUpPos->GetCopyNo(); //just to remove the warning during compiling

        //start down spacers
        spacerRmin = startingSpacerRmin;
        spacerRmax = startingSpacerRmax-FwBoardThick;

        sprintf(tShapeName,"spacerDwnS%dR%d",0,-1);
        sprintf(spacerVolName,"spacerDwnVolS%dR%d",0,-1);
        G4VSolid* spacerDwnShape = new G4Tubs(tShapeName, spacerRmin,spacerRmax,spacerBaseWidth,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   spacerDwnVol = new G4LogicalVolume(spacerDwnShape,motherMat,spacerVolName,0,0,0);

        sprintf(tShapeName,"spacerDwnBotS%dR%d",0,-1);
        sprintf(tVolName,"spacerDwnBotVolS%dR%d",0,-1);
        G4VSolid* spacerDwnBotShape = new G4Tubs(tShapeName, spacerRmin,spacerRmin+spacerBaseThick,spacerBaseWidth,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   spacerDwnBotVol = new G4LogicalVolume(spacerDwnBotShape,spacerMat,tVolName,0,0,0);
        G4VPhysicalVolume *spacerDwnBotPos = new G4PVPlacement(0,
                        G4ThreeVector(0,0,0),
                        spacerDwnBotVol,  // its logical volume
                        tVolName,         // its name
                        spacerDwnVol,     // its mother  volume
                        false,            // no boolean operations
                        0,                // copy number
                        detailedCheck);
        spacerDwnBotPos->GetCopyNo(); //just to remove the warning during compiling

        sprintf(tShapeName,"spacerDwnTopS%dR%d",0,-1);
        sprintf(tVolName,"spacerDwnTopVolS%dR%d",0,-1);
        G4VSolid* spacerDwnTopShape = new G4Tubs(tShapeName, spacerRmax-spacerBaseThick,spacerRmax,spacerBaseWidth,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   spacerDwnTopVol = new G4LogicalVolume(spacerDwnTopShape,spacerMat,tVolName,0,0,0);
        G4VPhysicalVolume *spacerDwnTopPos = new G4PVPlacement(0,
                        G4ThreeVector(0,0,0),
                        spacerDwnTopVol,  // its logical volume
                        tVolName,         // its name
                        spacerDwnVol,     // its mother  volume
                        false,            // no boolean operations
                        0,                // copy number
                        detailedCheck);
        spacerDwnTopPos->GetCopyNo(); //just to remove the warning during compiling

        sprintf(tShapeName,"spacerDwnSpokeS%dR%d",0,-1);
        sprintf(tVolName,"spacerDwnSpokeVolS%dR%d",0,-1);
        G4VSolid* spacerDwnSpokeShape = new G4Tubs(tShapeName, spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,spacerCoreThick,-0.5*spokeAngle,spokeAngle);
        G4LogicalVolume*   spacerDwnSpokeVol = new G4LogicalVolume(spacerDwnSpokeShape,spacerMat,tVolName,0,0,0);
        for (int iSpSpoke=0; iSpSpoke<spacerCoreNumSpokes; ++iSpSpoke){
                G4VPhysicalVolume *spacerDwnSpokePos = new G4PVPlacement(*spacerSpokesRot.at(iSpSpoke),
                                spacerDwnSpokeVol,  // its logical volume
                                tVolName,           // its name
                                spacerDwnVol,       // its mother  volume
                                false,              // no boolean operations
                                iSpSpoke,           // copy number
                                detailedCheck);
                spacerDwnSpokePos->GetCopyNo(); //just to remove the warning during compiling
        }

        G4VPhysicalVolume *spacerDwnPos = new G4PVPlacement(0,
                        spacerRelPos,
                        spacerDwnVol,     // its logical volume
                        spacerVolName,    // its name
                        anchorVol,        // its mother  volume
                        false,            // no boolean operations
                        0,                // copy number
                        detailedCheck);
        spacerDwnPos->GetCopyNo(); //just to remove the warning during compiling
        //end down spacers

        G4VPhysicalVolume *anchorPos = new G4PVPlacement(0,
                        G4ThreeVector(0,0,0),
                        anchorVol,        // its logical volume
                        ancrVolName,      // its name
                        localMother,      // its mother  volume
                        false,            // no boolean operations
                        0,                // copy number
                        detailedCheck);
        anchorPos->GetCopyNo(); //just to remove the warning during compiling
       //------------- end starting spacer and guard field wires -------------

        //------------- start finishing spacer and guard field wires -------------
        finishingSpacerRmin = maxR;
        finishingSpacerRmax = ((G4Tubs *)localMother->GetSolid())->GetOuterRadius();

        sprintf(tShapeName,"anchorS%dR%d",itracker->nSuperLayers(),-1);
        sprintf(ancrVolName,"anchorVolUpS%dR%d",itracker->nSuperLayers(),-1);
        G4VSolid* anchorShapeUp = new G4Tubs(tShapeName, finishingSpacerRmin,finishingSpacerRmax,motherDz,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   anchorVolUp = new G4LogicalVolume(anchorShapeUp,motherMat,ancrVolName,0,0,0);

        sprintf(tShapeName,"brdFldDwnS%dR%d",itracker->nSuperLayers(),-1);
        sprintf(tVolName,"brdFldDwnVolS%dR%d",itracker->nSuperLayers(),-1);
        G4VSolid* brdFldDwnShape = new G4Tubs(tShapeName, finishingSpacerRmin,finishingSpacerRmin+FwBoardThick,FwBoardWidth,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   brdFldDwnVol = new G4LogicalVolume(brdFldDwnShape,FwBoardMat,tVolName,0,0,0);
        G4VPhysicalVolume *brdFldDwnPos = new G4PVPlacement(0,
                        fwRelPos,
                        brdFldDwnVol,     // its logical volume
                        tVolName,         // its name
                        anchorVolUp,      // its mother  volume
                        false,            // no boolean operations
                        0,                // copy number
                        detailedCheck);
        brdFldDwnPos->GetCopyNo(); //just to remove the warning during compiling

        //start up spacer
        spacerRmin = finishingSpacerRmin+FwBoardThick;
        spacerRmax = finishingSpacerRmax;

        sprintf(tShapeName,"spacerUpS%dR%d",itracker->nSuperLayers(),-1);
        sprintf(spacerVolName,"spacerUpVolS%dR%d",itracker->nSuperLayers(),-1);
        G4VSolid* spacerUpShape = new G4Tubs(tShapeName, spacerRmin,spacerRmax,spacerBaseWidth,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   spacerUpVol = new G4LogicalVolume(spacerUpShape,motherMat,spacerVolName,0,0,0);

        sprintf(tShapeName,"spacerUpBotS%dR%d",itracker->nSuperLayers(),-1);
        sprintf(tVolName,"spacerUpBotVolS%dR%d",itracker->nSuperLayers(),-1);
        G4VSolid* spacerUpBotShape = new G4Tubs(tShapeName, spacerRmin,spacerRmin+spacerBaseThick,spacerBaseWidth,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   spacerUpBotVol = new G4LogicalVolume(spacerUpBotShape,spacerMat,tVolName,0,0,0);
        G4VPhysicalVolume *spacerUpBotPos = new G4PVPlacement(0,
                        G4ThreeVector(0,0,0),
                        spacerUpBotVol,   // its logical volume
                        tVolName,         // its name
                        spacerUpVol,      // its mother  volume
                        false,            // no boolean operations
                        0,                // copy number
                        detailedCheck);
        spacerUpBotPos->GetCopyNo(); //just to remove the warning during compiling

        sprintf(tShapeName,"spacerUpTopS%dR%d",itracker->nSuperLayers(),-1);
        sprintf(tVolName,"spacerUpTopVolS%dR%d",itracker->nSuperLayers(),-1);
        G4VSolid* spacerUpTopShape = new G4Tubs(tShapeName, spacerRmax-spacerBaseThick,spacerRmax,spacerBaseWidth,0.0,360.0*CLHEP::degree);
        G4LogicalVolume*   spacerUpTopVol = new G4LogicalVolume(spacerUpTopShape,spacerMat,tVolName,0,0,0);
        G4VPhysicalVolume *spacerUpTopPos = new G4PVPlacement(0,
                        G4ThreeVector(0,0,0),
                        spacerUpTopVol,   // its logical volume
                        tVolName,         // its name
                        spacerUpVol,      // its mother  volume
                        false,            // no boolean operations
                        0,                // copy number
                        detailedCheck);
        spacerUpTopPos->GetCopyNo(); //just to remove the warning during compiling

        sprintf(tShapeName,"spacerUpSpokeS%dR%d",itracker->nSuperLayers(),-1);
        sprintf(tVolName,"spacerUpSpokeVolS%dR%d",itracker->nSuperLayers(),-1);
        G4VSolid* spacerUpSpokeShape = new G4Tubs(tShapeName, spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,spacerCoreThick,-0.5*spokeAngle,spokeAngle);
        G4LogicalVolume*   spacerUpSpokeVol = new G4LogicalVolume(spacerUpSpokeShape,spacerMat,tVolName,0,0,0);
        for (int iSpSpoke=0; iSpSpoke<spacerCoreNumSpokes; ++iSpSpoke){
                G4VPhysicalVolume *spacerUpSpokePos = new G4PVPlacement(*spacerSpokesRot.at(iSpSpoke),
                                spacerUpSpokeVol,   // its logical volume
                                tVolName,           // its name
                                spacerUpVol,        // its mother  volume
                                false,              // no boolean operations
                                iSpSpoke,           // copy number
                                detailedCheck);
                spacerUpSpokePos->GetCopyNo(); //just to remove the warning during compiling
        }

        G4VPhysicalVolume *spacerUpPos = new G4PVPlacement(0,
                        spacerRelPos,
                        spacerUpVol,      // its logical volume
                        spacerVolName,    // its name
                        anchorVolUp,      // its mother  volume
                        false,            // no boolean operations
                        0,                // copy number
                        detailedCheck);
        spacerUpPos->GetCopyNo(); //just to remove the warning during compiling
        //end up spacer

        G4VPhysicalVolume *anchorPosUp = new G4PVPlacement(0,
                        G4ThreeVector(0,0,0),
                        anchorVolUp,      // its logical volume
                        ancrVolName,      // its name
                        localMother,      // its mother  volume
                        false,            // no boolean operations
                        0,                // copy number
                        detailedCheck);
        anchorPosUp->GetCopyNo(); //just to remove the warning during compiling
        //------------- end finishing spacer and guard field wires -------------


        //delete [] spacerSpokesRot;
        for (int iSpSpoke=0; iSpSpoke<spacerCoreNumSpokes; ++iSpSpoke){
                delete spacerSpokesRot.at(iSpSpoke);
        }

}

void ITrackerBuilder::constructSignalCables(G4LogicalVolume* localMother, SimpleConfig const& config) {
        GeomHandle<ITracker> itracker;

        double cableFractionOnREP = config.getDouble("itracker.cableFractionOnREP",0.0);
        if (cableFractionOnREP>1.0) { cableFractionOnREP=1.0; }

        bool isUpStream = localMother->GetName().contains("_L");
        std::string sideName;
        bool locateCables(false);
        int NSlayerToSkip=0;
        if (isUpStream) {
                sideName="_L";
                if (cableFractionOnREP>0) {
                        locateCables = true;
                        NSlayerToSkip=(int) 1.0/cableFractionOnREP-1;
                }
        }
        else {
                sideName="_R";
                locateCables = true;
        }

        double minR, maxR;
        char tShapeName[50], tVolName[50];
        double motherDz = ((G4Tubs *)localMother->GetSolid())->GetZHalfLength();

        minR = ((G4Tubs *)localMother->GetSolid())->GetInnerRadius();
        maxR = ((G4Tubs *)localMother->GetSolid())->GetOuterRadius();
        G4Material* motherMat = localMother->GetMaterial();

        G4Material* cableDielMat =  findMaterialOrThrow( config.getString("itracker.cableDielMaterial") );
        G4Material* cableWireMat =  findMaterialOrThrow( config.getString("itracker.cableWireMaterial") );

        double cableDielThickness = config.getDouble("itracker.cableDielThickness");
        double cableDielWidth = config.getDouble("itracker.cableDielWidth");
        double cableWireDiameter = config.getDouble("itracker.cableWireDiameter");
        double cableWirePitch = config.getDouble("itracker.cableWirePitch");
        double cableWireRadius = 0.5*cableWireDiameter;

        double totalCableThickness = /*2.0**/cableDielThickness+cableWireDiameter;
        double ncellLayers = itracker->nSuperLayers()*itracker->nRing();

        if ( (totalCableThickness*(ncellLayers+1)*0.5)>motherDz ) { throw cet::exception("GEOM") <<"The Signal Cables exceed its mother volume dim\n"; }

        double cableZStep = (2.0*motherDz-totalCableThickness)/(ncellLayers-1);

        int nWirePerCable = (int)((cableDielWidth+0.001)/cableWirePitch);

        if (locateCables) {

                int countSkippedLayer=-1;

                sprintf(tShapeName,"cableCont%s",sideName.c_str());
                sprintf(tVolName,"cableContVol%s",sideName.c_str());
                double cableAngle = cableDielWidth/minR;
                G4VSolid* cableContShape = new G4Tubs(tShapeName, minR,maxR,motherDz,-0.5*cableAngle,cableAngle);
                G4LogicalVolume*   cableContVol = new G4LogicalVolume(cableContShape,motherMat,tVolName,0,0,0);

                cableDielThickness *= 0.5;
                cableDielWidth *= 0.5;
                double maxRCable = std::sqrt( maxR*maxR - cableDielWidth*cableDielWidth ); //max height of the cable if it is described like a box
                double calbeHalfHeight, cableZmin, calbeHalfOrzLeng;
                double iWireRelPos;

                double layerRmin, layerRmax, swR, tmpZincr;
                int superlayer, iring, icellLayer=0;
                for (int iSl = 0; iSl < itracker->nSuperLayers(); iSl++){

                        SuperLayer *SLayer = itracker->getSuperLayer(iSl);

                        for (int iLy=0; iLy < SLayer->nLayers(); iLy++ ){

                                VolumeInfo LayerInfo;
                                boost::shared_ptr<ITLayer> ily = SLayer->getLayer(iLy);
                                superlayer = ily->Id().getSuperLayer();
                                iring = ily->Id().getLayer();
                                if (ily->nCells()>0) {

                                        ++countSkippedLayer;
                                        if (isUpStream && countSkippedLayer>0&&countSkippedLayer<=NSlayerToSkip) {
                                                ++icellLayer;
                                                if (countSkippedLayer==NSlayerToSkip) { countSkippedLayer=-1; }
                                                continue;
                                        }

                                        layerRmin = ily->getDetail()->centerInnerRadiusRing();
                                        tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleInnerRing());
                                        layerRmin = std::sqrt( layerRmin*layerRmin + tmpZincr*tmpZincr);

                                        layerRmax = ily->getDetail()->centerOuterRadiusRing();
                                        tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleOuterRing());
                                        layerRmax = std::sqrt( layerRmax*layerRmax + tmpZincr*tmpZincr);

                                        swR = (layerRmax+layerRmin)*0.5;


                                        sprintf(tShapeName,"signCableInsulS%dR%d%s",superlayer,iring,sideName.c_str());
                                        sprintf(tVolName,"signCableInsulVolS%dR%d%s",superlayer,iring,sideName.c_str());
                                        calbeHalfHeight = 0.5*(maxRCable-swR);
                                        G4VSolid* signCableInsulShape = new G4Box(tShapeName, cableDielWidth, cableDielThickness, calbeHalfHeight);
                                        G4LogicalVolume*   signCableInsulVol = new G4LogicalVolume(signCableInsulShape,cableDielMat,tVolName,0,0,0);
                                        cableZmin = cableZStep*(ncellLayers-icellLayer-1)-motherDz;
                                        HepGeom::Transform3D baseCablePos ( HepGeom::TranslateX3D(swR+calbeHalfHeight) * HepGeom::TranslateZ3D(cableZmin+cableDielThickness) * HepGeom::RotateY3D( 90*CLHEP::degree ) * HepGeom::RotateZ3D( 90*CLHEP::degree ) );

                                        G4VPhysicalVolume *signCableInsulPos = new G4PVPlacement( baseCablePos,
                                                        signCableInsulVol,  // its logical volume
                                                        tVolName,           // its name
                                                        cableContVol,       // its mother  volume
                                                        false,              // no boolean operations
                                                        0,                  // copy number
                                                        detailedCheck);
                                        signCableInsulPos->GetCopyNo(); //just to remove the warning during compiling

                                        sprintf(tShapeName,"signCableWireS%dR%d%s",superlayer,iring,sideName.c_str());
                                        sprintf(tVolName,"signCableWireVolS%dR%d%s",superlayer,iring,sideName.c_str());
                                        G4VSolid* signCableWireShape = new G4Tubs(tShapeName, 0.0, cableWireRadius, calbeHalfHeight, 0.0, 360*CLHEP::degree);
                                        G4LogicalVolume*   signCableWireVol = new G4LogicalVolume(signCableWireShape,cableWireMat,tVolName,0,0,0);
                                        HepGeom::Transform3D baseCableWirePos( HepGeom::TranslateZ3D(cableDielThickness+cableWireRadius)*baseCablePos );

                                        iWireRelPos = -cableDielWidth+cableWirePitch*0.5;
                                        for (int iWirePerCable=0; iWirePerCable<nWirePerCable; ++iWirePerCable) {
                                                HepGeom::TranslateY3D wireRelPos(iWireRelPos);
                                                G4VPhysicalVolume *signCableWirePos = new G4PVPlacement( HepGeom::Transform3D( wireRelPos*baseCableWirePos ),
                                                                signCableWireVol,   // its logical volume
                                                                tVolName,           // its name
                                                                cableContVol,       // its mother  volume
                                                                false,              // no boolean operations
                                                                iWirePerCable,      // copy number
                                                                detailedCheck);
                                                signCableWirePos->GetCopyNo(); //just to remove the warning during compiling
                                                iWireRelPos+=cableWirePitch;
                                        }

                                        calbeHalfOrzLeng = 0.5*(motherDz+cableZmin);
                                        if (calbeHalfOrzLeng>0.0) {
                                                sprintf(tShapeName,"signCableInsulOrznS%dR%d%s",superlayer,iring,sideName.c_str());
                                                sprintf(tVolName,"signCableInsulOrznVolS%dR%d%s",superlayer,iring,sideName.c_str());
                                                G4VSolid* signCableInsulOrznShape = new G4Box(tShapeName, cableDielWidth, cableDielThickness, calbeHalfOrzLeng);
                                                G4LogicalVolume*   signCableInsulOrznVol = new G4LogicalVolume(signCableInsulOrznShape,cableDielMat,tVolName,0,0,0);
                                                HepGeom::Transform3D baseCablePosOrz ( HepGeom::TranslateX3D(swR+cableDielThickness) * HepGeom::TranslateZ3D(-motherDz+calbeHalfOrzLeng) * HepGeom::RotateZ3D( 90*CLHEP::degree ) );

                                                G4VPhysicalVolume *signCableInsulOrznPos = new G4PVPlacement( baseCablePosOrz,
                                                                signCableInsulOrznVol,  // its logical volume
                                                                tVolName,           // its name
                                                                cableContVol,       // its mother  volume
                                                                false,              // no boolean operations
                                                                0,                  // copy number
                                                                detailedCheck);
                                                signCableInsulOrznPos->GetCopyNo(); //just to remove the warning during compiling

                                                sprintf(tShapeName,"signCableWireOrznS%dR%d%s",superlayer,iring,sideName.c_str());
                                                sprintf(tVolName,"signCableWireOrznVolS%dR%d%s",superlayer,iring,sideName.c_str());
                                                G4VSolid* signCableWireOrznShape = new G4Tubs(tShapeName, 0.0, cableWireRadius, calbeHalfOrzLeng, 0.0, 360*CLHEP::degree);
                                                G4LogicalVolume*   signCableWireOrznVol = new G4LogicalVolume(signCableWireOrznShape,cableWireMat,tVolName,0,0,0);
                                                HepGeom::Transform3D baseCableWirePosOrz( HepGeom::TranslateX3D(-(cableDielThickness+cableWireRadius))*baseCablePosOrz );


                                                iWireRelPos = -cableDielWidth+cableWirePitch*0.5;
                                                for (int iWirePerCable=0; iWirePerCable<nWirePerCable; ++iWirePerCable) {
                                                        HepGeom::TranslateY3D wireRelPos(iWireRelPos);
                                                        G4VPhysicalVolume *signCableWireOrznPos = new G4PVPlacement( HepGeom::Transform3D( wireRelPos*baseCableWirePosOrz ),
                                                                        signCableWireOrznVol,  // its logical volume
                                                                        tVolName,              // its name
                                                                        cableContVol,          // its mother  volume
                                                                        false,                 // no boolean operations
                                                                        iWirePerCable,         // copy number
                                                                        detailedCheck);
                                                        signCableWireOrznPos->GetCopyNo(); //just to remove the warning during compiling
                                                        iWireRelPos+=cableWirePitch;
                                                }
                                        }

                                        ++icellLayer;
                                }

                        }

                }


                //---------- positioning cables pack for each sector ----------
                int spdWebSpokesNumber = config.getInt("itracker.spdWebSpokesNumber");
                //HepGeom::Transform3D baseSpokePos ();
                double spokeAngleStep = CLHEP::twopi/((double)spdWebSpokesNumber);
                for (int iSpokes=0; iSpokes<spdWebSpokesNumber; ++iSpokes) {
                        G4VPhysicalVolume *cableContPos = new G4PVPlacement(
                                        HepGeom::Transform3D( HepGeom::RotateZ3D( ((double)iSpokes)*spokeAngleStep ) ),
                                        cableContVol,         // its logical volume
                                        tVolName,             // its name
                                        localMother,          // its mother  volume
                                        false,                // no boolean operations
                                        iSpokes,              // copy number
                                        detailedCheck);
                        cableContPos->GetCopyNo(); //just to remove the warning during compiling
                }

        }

}

void ITrackerBuilder::constructHvCables(G4LogicalVolume* localMother, SimpleConfig const& config) {
        GeomHandle<ITracker> itracker;

        bool isDownStream = localMother->GetName().contains("_R");

        double minR, maxR;
        char tShapeName[50], tVolName[50];
        double motherDz = ((G4Tubs *)localMother->GetSolid())->GetZHalfLength();

        minR = ((G4Tubs *)localMother->GetSolid())->GetInnerRadius();
        maxR = ((G4Tubs *)localMother->GetSolid())->GetOuterRadius();
        G4Material* motherMat = localMother->GetMaterial();

        std::vector<double> cableRadii;
        std::vector<std::string> cableMatList;
        config.getVectorDouble("itracker.hvCableShellsThicknesses",cableRadii);
        config.getVectorString("itracker.hvCableMaterials",cableMatList);
        if (cableRadii.size()!=cableMatList.size()) { throw cet::exception("GEOM") <<"The number of Radii is different from that of materials for the HV cable\n"; }
        std::vector<G4Material *> cableMats;
        double totalCableThickness = 0.0;
        int nShells = cableRadii.size();
        for (int iShell=0; iShell<nShells; ++iShell) {
                totalCableThickness += cableRadii.at(iShell);
                cableMats.push_back( findMaterialOrThrow( cableMatList.at(iShell) ) );
        }
        totalCableThickness*=2.0;
        double cableDielWidth = totalCableThickness;

        double ncellLayers = itracker->nSuperLayers()*itracker->nRing();

        if ( (totalCableThickness*ncellLayers*0.5)>motherDz ) { throw cet::exception("GEOM") <<"The Hv Cables exceed its mother volume dim\n"; }

        double cableZStep = totalCableThickness;
        totalCableThickness*=0.5;

        if (isDownStream) {

                sprintf(tShapeName,"hvCableCont");
                sprintf(tVolName,"hvCableContVol");
                double cableAngle = cableDielWidth/minR;
                G4VSolid* hvCableContShape = new G4Tubs(tShapeName, minR,maxR,motherDz,-0.5*cableAngle,cableAngle);
                G4LogicalVolume*   hvCableContVol = new G4LogicalVolume(hvCableContShape,motherMat,tVolName,0,0,0);

                cableDielWidth *= 0.5;
                double maxRCable = std::sqrt( maxR*maxR - cableDielWidth*cableDielWidth ); //max height of the cable if it is described like a box
                double calbeHalfHeight, cableZmin, calbeHalfOrzLeng;
                double iCableShellRmin, iCableShellRmax;

                double layerRmin, layerRmax, swR, tmpZincr;
                int superlayer, iring, icellLayer=0;
                for (int iSl = 0; iSl < itracker->nSuperLayers(); iSl++){

                        SuperLayer *SLayer = itracker->getSuperLayer(iSl);

                        for (int iLy=0; iLy < SLayer->nLayers(); iLy++ ){

                                VolumeInfo LayerInfo;
                                boost::shared_ptr<ITLayer> ily = SLayer->getLayer(iLy);
                                superlayer = ily->Id().getSuperLayer();
                                iring = ily->Id().getLayer();
                                if (ily->nCells()>0) {
                                        layerRmin = ily->getDetail()->centerInnerRadiusRing();
                                        tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleInnerRing());
                                        layerRmin = std::sqrt( layerRmin*layerRmin + tmpZincr*tmpZincr);

                                        layerRmax = ily->getDetail()->centerOuterRadiusRing();
                                        tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleOuterRing());
                                        layerRmax = std::sqrt( layerRmax*layerRmax + tmpZincr*tmpZincr);

                                        swR = (layerRmax+layerRmin)*0.5;

                                        iCableShellRmin=0.0;
                                        iCableShellRmax=0.0;
                                        for (int iShell=0; iShell<nShells; ++iShell) {
                                                iCableShellRmax+=cableRadii.at(iShell);

                                                sprintf(tShapeName,"hvCableShellS%dR%d_%d",superlayer,iring,iShell);
                                                sprintf(tVolName,"hvCableShellVolS%dR%d_%d",superlayer,iring,iShell);
                                                calbeHalfHeight = 0.5*(maxRCable-swR);
                                                G4VSolid* hvCableShellShape = new G4Tubs(tShapeName, iCableShellRmin, iCableShellRmax, calbeHalfHeight, 0.0, 360*CLHEP::degree);
                                                G4LogicalVolume*   hvCableShellVol = new G4LogicalVolume(hvCableShellShape,cableMats.at(iShell),tVolName,0,0,0);
                                                cableZmin = cableZStep*(ncellLayers-icellLayer-1)-motherDz;
                                                HepGeom::Transform3D baseCablePos ( HepGeom::TranslateX3D(swR+calbeHalfHeight) * HepGeom::TranslateZ3D(cableZmin+totalCableThickness) * HepGeom::RotateY3D( 90*CLHEP::degree ) );

                                                G4VPhysicalVolume *hvCableShellPos = new G4PVPlacement( baseCablePos,
                                                                hvCableShellVol,    // its logical volume
                                                                tVolName,           // its name
                                                                hvCableContVol,     // its mother  volume
                                                                false,              // no boolean operations
                                                                0,                  // copy number
                                                                detailedCheck);
                                                hvCableShellPos->GetCopyNo(); //just to remove the warning during compiling

                                                calbeHalfOrzLeng = 0.5*(motherDz+cableZmin);
                                                if (calbeHalfOrzLeng>0.0) {
                                                        sprintf(tShapeName,"hvCableShellOrznS%dR%d_%d",superlayer,iring,iShell);
                                                        sprintf(tVolName,"hvCableShellOrznVolS%dR%d_%d",superlayer,iring,iShell);
                                                        G4VSolid* hvCableShellOrznShape = new G4Tubs(tShapeName, iCableShellRmin, iCableShellRmax, calbeHalfOrzLeng, 0.0, 360*CLHEP::degree);
                                                        G4LogicalVolume*   hvCableShellOrznVol = new G4LogicalVolume(hvCableShellOrznShape,cableMats.at(iShell),tVolName,0,0,0);
                                                        HepGeom::Transform3D baseCablePosOrz ( HepGeom::TranslateX3D(swR) * HepGeom::TranslateZ3D(-motherDz+calbeHalfOrzLeng) );

                                                        G4VPhysicalVolume *hvCableShellOrznPos = new G4PVPlacement( baseCablePosOrz,
                                                                        hvCableShellOrznVol,  // its logical volume
                                                                        tVolName,             // its name
                                                                        hvCableContVol,       // its mother  volume
                                                                        false,                // no boolean operations
                                                                        0,                    // copy number
                                                                        detailedCheck);
                                                        hvCableShellOrznPos->GetCopyNo(); //just to remove the warning during compiling
                                                }

                                                iCableShellRmin = iCableShellRmax;
                                        }

                                        ++icellLayer;
                                }

                        }

                }


                //---------- positioning cables pack for each sector ----------
                int spdWebSpokesNumber = config.getInt("itracker.spdWebSpokesNumber");
                double spokeAngleStep = CLHEP::twopi/((double)spdWebSpokesNumber);
                HepGeom::Transform3D baseSpokePos (HepGeom::RotateZ3D( 0.5*spokeAngleStep ));
                for (int iSpokes=0; iSpokes<spdWebSpokesNumber; ++iSpokes) {
                        G4VPhysicalVolume *hvCableContPos = new G4PVPlacement(
                                        HepGeom::Transform3D( HepGeom::RotateZ3D( ((double)iSpokes)*spokeAngleStep )*baseSpokePos ),
                                        hvCableContVol,        // its logical volume
                                        tVolName,              // its name
                                        localMother,           // its mother  volume
                                        false,                 // no boolean operations
                                        iSpokes,               // copy number
                                        detailedCheck);
                        hvCableContPos->GetCopyNo(); //just to remove the warning during compiling
                }

        }

}

} // end namespace mu2e
