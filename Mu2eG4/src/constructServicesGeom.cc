//
// David Norvil Brown, University of Louisville, March 2016
//
// Build Electronic Racks, Pipes, etc.  (AKA Services)

#include "Offline/Mu2eG4/inc/constructServicesGeom.hh"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "cetlib_except/exception.h"

// Include each ServicesGeom element here...
#include "Offline/ServicesGeom/inc/Pipe.hh"
#include "Offline/ServicesGeom/inc/ElectronicRack.hh"

// etc...
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/GeomPrimitives/inc/Torus.hh"
#include "Offline/GeomPrimitives/inc/TorusParams.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestTorus.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/GeneralUtilities/inc/OrientationResolver.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Orb.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4TwoVector.hh"
#include "Geant4/G4NistManager.hh"

#include "Geant4/G4LogicalVolume.hh"

#include <vector>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <string>
#include <array>

namespace mu2e {


  void constructServicesGeom(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<Pipe>   pipeSet;
    GeomHandle<ElectronicRack>   rackSet;

    // Utility for converting orientations to rotations
    OrientationResolver* OR = new OrientationResolver();

    // Get config info
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "services",       "services");
    geomOptions->loadEntry( config, "electronicRack", "electronicRack");
    geomOptions->loadEntry( config, "pipe",           "pipe");


    const bool electronicRackIsVisible = geomOptions->isVisible("electronicRack");
    const bool electronicRackIsSolid   = geomOptions->isSolid("electronicRack");
    const bool pipeIsVisible           = geomOptions->isVisible("pipe");
    const bool pipeIsSolid             = geomOptions->isSolid("pipe");
    const bool doSurfaceCheck          = geomOptions->doSurfaceCheck("services");
    const bool forceAuxEdgeVisible     = geomOptions->forceAuxEdgeVisible("services");
    const bool placePV                 = geomOptions->placePV("services");


    // *******************************************************
    // ==> Make Pipes <========
    // *******************************************************

    int                        version = pipeSet->getVersion();

    // Load up the vectors needed for building.
    std::vector<int>            nPipes = pipeSet->getNPipes();
    std::vector<int>            nComps = pipeSet->getNComponentsInPipe();
    std::vector<double>         pLeng  = pipeSet->getLengths();
    std::vector<std::string>    pFlav  = pipeSet->getFlavor();
    std::vector<std::string>    pFill  = pipeSet->getFillMaterialNames();

    std::vector<std::vector<CLHEP::Hep3Vector> > pCent = pipeSet->getCentersOfPipes();
    std::vector<std::vector<std::string> > pOrient = pipeSet->getOrientations();

    std::vector<std::vector<double> > cInRad = pipeSet->getInnerRads();
    std::vector<std::vector<double> > cOutRad = pipeSet->getOuterRads();
    std::vector<std::vector<std::string> > cMats = pipeSet->getMaterialNames();
    std::vector<std::vector<double> > cUOff = pipeSet->getUOffsets();
    std::vector<std::vector<double> > cVOff = pipeSet->getVOffsets();

    // ***************************************************
    // *** Loop over the types and construct all *********
    // ***************************************************

    int whichPipe = 0;
    for ( unsigned int it = 0; it < nPipes.size(); it++ ) {

      int nPipe = nPipes[it];
      int nComp = nComps[it];
      double len = 0.0;
      if (version == 1) len = pLeng[it];
      std::string flav = pFlav[it];
      bool isBend = false;
      if ( flav != "straight" ) {
        isBend = true;
      } // end of if not straight

      std::string fillMat = pFill[it];

      // *** Now loop over the individual pipes ***
      for ( int ip = 0; ip < nPipe; ip++ ) {

        if ( version > 1 ) len = pLeng[whichPipe];
        whichPipe++;
        // Make the container ("mother") volume for the type.
        TubsParams  pipeParams(0.0,cOutRad[it][0],len/2.0);
        TorusParams bendParams(0.0,cOutRad[it][0],len);

        std::ostringstream motherTubeName;
        if ( isBend ) {
          motherTubeName << "pBendType" << it+1 << "PBend" << ip+1;
        } else {
          motherTubeName << "pipeType" << it+1 << "Pipe" << ip+1;
        }

        std::ostringstream motherName;
        if ( isBend ) {
          motherName << "pBendLogicVolType" << it+1 << "PBend" << ip+1;
        } else {
          motherName << "pipeLogicVolType" << it+1 << "Pipe" << ip+1;
        }

        CLHEP::HepRotation* pipeRotat = new
          CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
        OR->getRotationFromOrientation( *pipeRotat, pOrient[it][ip] );

        CLHEP::Hep3Vector pipePosInMu2e = pCent[it][ip]
          - parent.centerInMu2e();


        VolumeInfo motherLogVol;
        if ( isBend ) {
          motherLogVol = nestTorus(motherName.str(),
                                   bendParams,
                                   findMaterialOrThrow(fillMat),
                                   pipeRotat,
                                   pipePosInMu2e,
                                   parent,
                                   ip,
                                   pipeIsVisible,
                                   G4Colour::Magenta(),
                                   pipeIsSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck);
        } else {
          motherLogVol = nestTubs(motherName.str(),
                                  pipeParams,
                                  findMaterialOrThrow(fillMat),
                                  pipeRotat,
                                  pipePosInMu2e,
                                  parent.logical,
                                  ip,
                                  pipeIsVisible,
                                  G4Colour::Magenta(),
                                  pipeIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck);
        } // end of making mother logical if isBend, etc

        for ( int ic = 0; ic < nComp; ic++ ) {
          // Create a tube for each "component" pipe and embed in the mother
          // Logical volume
          TubsParams componentParams(cInRad[it][ic],cOutRad[it][ic],len/2.0);
          TorusParams compTorParams(cInRad[it][ic],cOutRad[it][ic],len-cUOff[it][ic]);

          std::ostringstream compName;
          if ( isBend ) {
            compName << "pBendType" << it+1 << "PBend" << ip+1 <<
            "Component" << ic+1;
          } else {
            compName << "pipeType" << it+1 << "Pipe" << ip+1 <<
            "Component" << ic+1;
          }

          // Until we place it, the mother should be centered at 0,0,0
          // So place relative to that using u- and v-components and w=0.
          CLHEP::Hep3Vector posComp(cUOff[it][ic],cVOff[it][ic],0.0);
          CLHEP::Hep3Vector tPosComp(0.0,0.0,cVOff[it][ic]);
          if ( isBend ) {
            nestTorus ( compName.str(), compTorParams,
                        findMaterialOrThrow(cMats[it][ic]),
                        0, tPosComp, motherLogVol,
                        ic,
                        pipeIsVisible,
                        G4Colour::Magenta(),
                        pipeIsSolid,
                        forceAuxEdgeVisible,
                        placePV,
                        doSurfaceCheck);
          } else {
            nestTubs( compName.str(),  componentParams,
                      findMaterialOrThrow(cMats[it][ic]),
                      0, posComp, motherLogVol,
                      ic,
                      pipeIsVisible,
                      G4Colour::Magenta(),
                      pipeIsSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
          } // end of if isBend
        } // end of loop over components

      } // end of loop over pipes of this type.

    } // end of loop over pipe types

    //----------------------------------------------------------------
    // Electronics Rack Boxes <=====
    //----------------------------------------------------------------
    // Fill the vectors of information about the boxes that are
    // the racks

    std::vector<std::vector<double> > dimsER = rackSet->getBoxDimensions();
    std::vector<std::string> matsER = rackSet->getMaterialNames();
    std::vector<CLHEP::Hep3Vector> sitesER = rackSet->getCentersOfBoxes();
    std::vector<std::string> orientsER = rackSet->getOrientations();

    int nBoxER = dimsER.size();

    for(int i = 0; i < nBoxER; i++)
   {

        // Dimensions for this rack
        std::vector<double> lwhsER = dimsER[i];

        //  Make the name of the box
        std::ostringstream nameER;
        nameER << "ElectronicRackBox_" << i+1 ;

        // Make the needed rotation by parsing orientation
        CLHEP::HepRotation* itsRotatER= new
          CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
        std::string orientInitER = orientsER[i];

        OR->getRotationFromOrientation(*itsRotatER, orientInitER);

        // Build each box here

        nestBox( nameER.str(), lwhsER, findMaterialOrThrow(matsER[i]),
                 itsRotatER, sitesER[i]-parent.centerInMu2e(),
                 parent.logical,
                 0,
                 electronicRackIsVisible,
                 G4Colour::Magenta(),
                 electronicRackIsSolid,
                 forceAuxEdgeVisible,
                 placePV,
                 doSurfaceCheck);
   } // end loop over ElectronicRack boxes


  } // end of constructServicesGeom fn

} // namespace mu2e

