/*
 * ITrackerBuilder.cpp
 *
 *  Created on: Feb 11, 2010
 *      Author: tassiell
 */

// C++ includes
#include <iostream>
#include <sstream>

#include <boost/regex.hpp>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/ITrackerBuilder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/ITGasLayerSD_ExtWireData.hh"
#include "Mu2eG4/inc/ITGasLayerSD_v2.hh"
#include "Mu2eG4/inc/ITGasLayerSD_v3.hh"

// G4 includes
#include "G4Tubs.hh"
#include "G4Hype.hh"
#include "G4Sphere.hh"
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "globals.hh"

using namespace std;

namespace mu2e {

VolumeInfo ITrackerBuilder::constructTracker( G4LogicalVolume* mother, double zOff ){

        // Master geometry for the ITracker.
        GeomHandle<ITracker> itracker;

        VolumeInfo trackerInfo;

        // Make the mother volume for the ITracker.
        string trackerName("TrackerMother");

        double z0    = CLHEP::mm * itracker->z0();
        G4ThreeVector trackerOffset(0.,0.,z0-zOff);

        if ( itracker->rOut() >= ((G4Tubs *)mother->GetSolid())->GetOuterRadius() )
                throw cms::Exception("GEOM") <<"The ITracker doesn't fit inside the DS, please check the external radius of the ITracker\n";

        if (itracker->isExternal()) {
                throw cms::Exception("GEOM") <<"This GDML file option is temporarily disabled\n";
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

                G4Material* Vacuum = findMaterialOrThrow( "WAGVacuum" );
                trackerInfo.solid = new G4Tubs(itracker->name(),0.0,itracker->rOut(),itracker->maxEndCapDim(),0.0,360.0*CLHEP::degree);
                trackerInfo.logical = new G4LogicalVolume(trackerInfo.solid , Vacuum, trackerName,0,0,0);
                trackerInfo.logical->SetVisAttributes(visAtt);

                char shape[30], vol[30], shape_name_FD[30], shape_name_SD[30], vol_name_FD[30], vol_name_SD[30], wire_name[40];
                sprintf(shape_name_FD,"tube_Field");
                sprintf(shape_name_SD,"tube_Sense");

                int superlayer,iring;
                G4SDManager* SDman   = G4SDManager::GetSDMpointer();

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

                double outerWallInnerRadius;

                multimap<Wall::Walltype,boost::shared_ptr<Wall> >::iterator walls_it;
                for ( walls_it=itracker->getWalls()->begin() ; walls_it != itracker->getWalls()->end(); walls_it++ ) {
                        VolumeInfo WallInfo;
                        boost::shared_ptr<Wall> iwall = walls_it->second;
                        WallInfo = buildWall(iwall.get(),itracker->endcapType());
                        WallInfo.logical->SetVisAttributes(visAttWall);
                        WallInfo.physical = new G4PVPlacement(iwall->getPos(),
                                        WallInfo.logical,                // its logical volume
                                        iwall->getName(),                // its name
                                        trackerInfo.logical,                // its mother  volume
                                        false,                                        // no boolean operations
                                        0);                                                // copy number
                        if (iwall->getType() == Wall::outer) outerWallInnerRadius = iwall->getRmin();
                }

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

                                //cout<<ily->getDetail()->centerOuterRadiusRing()<<" "<<sqrt( pow(ily->getDetail()->centerOuterRadiusRing(),2) +
                                //                                                pow(ily->getDetail()->halfLength()*tan(ily->getDetail()->stereoAngleOuterRing()),2) ) <<" "<< outerWallInnerRadius<<endl;
                                //if ( sqrt( pow(ily->getDetail()->centerOuterRadiusRing(),2) +
                                //                pow(ily->getDetail()->halfLength()*tan(ily->getDetail()->stereoAngleOuterRing()),2) ) > outerWallInnerRadius )
                                //        throw cms::Exception("GEOM") <<"The ITracker layer "<<ily->Id()<<" doesn't fit inside the ITracker outer wall\n";

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
                                                0);                      // copy number

                                if (ily->getLayerType() != ITLayer::undefined) {
                                        ITGasLayerSD* glSD;
                                        if ( itracker->geomType()==ITracker::Hexagonal ) glSD = new ITGasLayerSD_v2( vol );
                                        else if ( itracker->geomType()==ITracker::Square ) glSD = new ITGasLayerSD_v3( vol );
                                        SDman->AddNewDetector( glSD );
                                        LayerInfo.logical->SetSensitiveDetector( glSD );
                                }

                                for ( int iFw=0; iFw < ily->nFieldWires(); iFw++){
                                        sprintf(vol_name_FD,"tubeFD_%d_%d",superlayer,iring);
                                        VolumeInfo FieldWireInfo;
                                        boost::shared_ptr<Wire> iwire = ily->getFWire(iFw);
                                        boost::shared_ptr<WireDetail> wdet = iwire->getDetail();
                                        FieldWireInfo = buildWire(wdet->outerRadius(),wdet->halfLength(),shape_name_FD,vol_name_FD,wdet->materialNames(),wdet->shellsThicknesses());
                                        sprintf(wire_name,"%s_%i",vol_name_FD,iwire->Id().getWire());
                                        FieldWireInfo.logical->SetVisAttributes(visAttFw);
                                        FieldWireInfo.physical = new G4PVPlacement(iwire->get3DTransfrom(),
                                                        FieldWireInfo.logical,         // its logical volume
                                                        wire_name,                     // its name
                                                        LayerInfo.logical,             // its mother  volume
                                                        false,                         // no boolean operations
                                                        iwire->Id().getWire());  // copy number
                                }

                                for ( int iSw=0; iSw < ily->nCells(); iSw++){
                                        sprintf(vol_name_SD,"tubeSD_%d_%d",superlayer,iring);
                                        VolumeInfo SenseWireInfo;
                                        boost::shared_ptr<Wire> iwire = ily->getCell(iSw)->getWire();
                                        boost::shared_ptr<WireDetail> wdet = iwire->getDetail();
                                        SenseWireInfo = buildWire(wdet->outerRadius(),wdet->halfLength(),shape_name_SD,vol_name_SD,wdet->materialNames(),wdet->shellsThicknesses());
                                        sprintf(wire_name,"%s_%i",vol_name_SD,iwire->Id().getWire());
                                        SenseWireInfo.logical->SetVisAttributes(visAttSw);
                                        SenseWireInfo.physical = new G4PVPlacement(iwire->get3DTransfrom(),
                                                        SenseWireInfo.logical,         // its logical volume
                                                        wire_name,                     // its name
                                                        LayerInfo.logical,             // its mother  volume
                                                        false,                         // no boolean operations
                                                        iwire->Id().getWire());  // copy number
                                }

                        }

                }




                trackerInfo.physical =  new G4PVPlacement( 0,
                                trackerOffset,
                                trackerInfo.logical,
                                trackerName,
                                mother,
                                0,
                                0);

        }

        return trackerInfo;

}

VolumeInfo ITrackerBuilder::buildWire(float radius, float length, char *shapeName, char *volName, const std::vector<std::string> &materialName, const std::vector<double> &thicknesses){

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
                        G4VPhysicalVolume *tphysWire = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,0),
                                        tlogicWire,      // its logical volume
                                        tVolName,             // its name
                                        wire.logical,    // its mother  volume
                                        false,                 // no boolean operations
                                        0);                  // copy number
                }

        }
        return wire;
}

VolumeInfo ITrackerBuilder::buildWall(Wall *wall, ITracker::EnCapType endcapType){

        VolumeInfo wallInfo;
        char volName[50], shapeName[50];
        sprintf(volName,"vol_%s",wall->getName().c_str());
        sprintf(shapeName,"shape_%s",wall->getName().c_str());
        if (wall->getType()==Wall::endcap && endcapType == ITracker::Spherical ) {
                wallInfo.solid = new G4Sphere( shapeName,wall->getRmin(),wall->getRmax(),wall->getSPhi(),wall->getDPhi(),wall->getSTheta(),wall->getDTheta() );
        } else {
                wallInfo.solid = new G4Tubs( shapeName,wall->getRmin(),wall->getRmax(),wall->getDz(),wall->getSPhi(),wall->getDPhi() );
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

                for (int ishell=0; ishell<nSub; ishell++){
                        sprintf(tShapeName,"%s_sub%i",shapeName,ishell);
                        sprintf(tVolName,"%s_sub%i",volName,ishell);
                        iRadius+=wall->getThicknesses()->at(ishell);
                        G4VSolid *tswall;
                        if (wall->getType()==Wall::endcap && endcapType == ITracker::Spherical ) {
                                tswall = new G4Sphere( tShapeName,oldRadius,iRadius,wall->getSPhi(),wall->getDPhi(),wall->getSTheta(),wall->getDTheta() );
                        } else {
                                tswall = new G4Tubs(tShapeName,oldRadius,iRadius,wall->getDz(),wall->getSPhi(),wall->getDPhi() );
                        }
                        //          cout<<tShapeName<<" "<<oldRadius<<" "<<iRadius<<" "<<length<<" "<<wall->getMaterialsName()->at(ishell)<<endl;
                        oldRadius=iRadius;

                        G4LogicalVolume *tlogicwall = new G4LogicalVolume(tswall,findMaterialOrThrow(wall->getMaterialsName()->at(ishell).c_str()),tVolName,0,0,0);
                        G4VPhysicalVolume *tphyswall = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,0),
                                        tlogicwall,      // its logical volume
                                        tVolName,             // its name
                                        wallInfo.logical,    // its mother  volume
                                        false,                 // no boolean operations
                                        0);                  // copy number
                }

        }
        return wallInfo;
}

} // end namespace mu2e
