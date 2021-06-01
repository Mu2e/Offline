// Function to build G4 rep of saddles for cryostats.
//
//
// David Norvil Brown, University of Louisville, 02 June 2017, based
// on constructExternalShielding, written November 2014
//
//

#include "Mu2eG4/inc/constructSaddles.hh"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "cetlib_except/exception.h"


#include "ExternalShieldingGeom/inc/Saddle.hh"

// etc...
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "GeneralUtilities/inc/OrientationResolver.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4Orb.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4TwoVector.hh"
#include "Geant4/G4NistManager.hh"

#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4LogicalVolume.hh"

#include <vector>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

namespace mu2e {


  void constructSaddles(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<Saddle> saddleSet;

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    // Utility for converting orientations to rotations
    OrientationResolver* OR = new OrientationResolver();

    // Get config info
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "saddle", "saddle");

    const bool saddleIsVisible     = geomOptions->isVisible("saddle");
    const bool saddleIsSolid       = geomOptions->isSolid("saddle");
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("saddle");
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck("saddle");
    const bool placePV             = geomOptions->placePV("saddle");


    //================================================================
    // OK, so in version 1, Saddles were just extruded boxes that,
    // when taken as a whole, looked like a bunch of saddles and
    // stands for (mainly) cryostats.
    // In version 2,3, Saddles are logical structures, each of which looks
    // like an individual saddle or stand.  Each saddle contains
    // an extruded box with holes and/or notches.
    //================================================================

    if ( saddleSet->getVersion() <= 3 ) {
      //----------------------------------------------------------------
      // Saddle Boxes <=====
      //----------------------------------------------------------------
      // Fill the vectors of information about the boxes making up the saddles.

      std::vector<std::vector<std::vector<double> > > outlSA = saddleSet->getOutlines();
      std::vector<double> lengSA               = saddleSet->getLengths();
      std::vector<std::vector<double> > tolsSA = saddleSet->getTolerances();
      std::vector<std::string> matsSA          = saddleSet->getMaterialNames();
      std::vector<CLHEP::Hep3Vector> sitesSA   = saddleSet->getCentersOfBoxes();
      std::vector<std::string> orientsSA       = saddleSet->getOrientations();
      std::vector<int> nHolesSA                = saddleSet->getNHoles();
      std::vector<int> nNotchesSA              = saddleSet->getNNotches();
      std::vector<int> holeIDSA                = saddleSet->getHoleIndices();
      std::vector<CLHEP::Hep3Vector> holeLocSA = saddleSet->getHoleLocations();
      std::vector<double> holeRadSA            = saddleSet->getHoleRadii();
      std::vector<double> holeLenSA            = saddleSet->getHoleLengths();
      std::vector<std::string> holeOrientsSA   = saddleSet->getHoleOrientations();
      std::vector<int> notchIDSA               = saddleSet->getNotchIndices();
      std::vector<CLHEP::Hep3Vector> notchLocSA= saddleSet->getNotchLocations();
      std::vector<std::vector<double> > notchDimSA = saddleSet->getNotchDimensions();

      int nBox = outlSA.size();

      for(int i = 0; i < nBox; i++)
        {
          // combine the tolerances with the outline dimensions.
          std::vector<G4TwoVector> itsOutline;
          std::vector<std::vector<double> > vertices = outlSA[i];
          double du = tolsSA[i][0];
          double dv = tolsSA[i][1];
          double dw = tolsSA[i][2];
          double hlen = lengSA[i];

          for ( unsigned int idim = 0; idim < vertices.size(); idim++ ) {
            if ( vertices[idim][0] > -10.0 * CLHEP::mm) {
              vertices[idim][0] += du;
            } else {
              vertices[idim][0] -= du;
            }
            if (vertices[idim][1] > -10.0 * CLHEP::mm ) {
              vertices[idim][1] += dv;
            } else {
              vertices[idim][1] -= dv;
            }
            G4TwoVector vertex( vertices[idim][0], vertices[idim][1] );
            itsOutline.push_back(vertex);
          }
          hlen += dw;

          //  Make the name of the box
          std::ostringstream name;
          name << "SaddleBox_" << i+1 ;

          // Make the needed rotation by parsing orientation
          std::string orientSAInit = orientsSA[i];

          CLHEP::HepRotation* itsSARotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          OR->getRotationFromOrientation( *itsSARotat, orientSAInit );

          // ****************************************
          // Now add holes and notches
          // ***************************************
          if ( nHolesSA[i] + nNotchesSA[i] > 0 ) {

            // This box has at least one window or notch.
            // Each is implemented as a
            // G4SubtractionSolid to allow for another volume placement
            // through it.  First go ahead and make the extrusion for the block.
            std::ostringstream name1;
            name1 << "SaddleBox_sub_" << i+1;

            G4ExtrudedSolid* block = new G4ExtrudedSolid( name1.str(),
                                                          itsOutline,
                                                          hlen,
                                                          G4TwoVector(0,0), 1.,
                                                          G4TwoVector(0,0), 1.);

            // Now create the volume for the block
            VolumeInfo extShieldVol(name.str(),
                                    sitesSA[i]-parent.centerInMu2e(),
                                    parent.centerInWorld);


            // Create a subtraction solid that will become the solid for the
            // volume.

            G4SubtractionSolid* aSolid = 0;

            // Find the index of the first hole, for accessing info lists...
            int hID = holeIDSA[i];

            // Now loop over the holes and make them.
            for ( int jHole = 0; jHole < nHolesSA[i]; jHole++ ) {
              // Now make the window (AKA "Hole")
              name << "h" << jHole+1;
              //            std::cout << __func__ << " making " << name.str() << std::endl;

              int thisHID = hID + jHole; // get pointer for right hole

              const TubsParams windparams(0.0,  //inner radius
                                          holeRadSA[thisHID], // outer
                                          holeLenSA[thisHID]  //obvious?
                                          );

              std::ostringstream name2;
              name2 << "SaddleBox" << i+1 << "hole" << jHole+1;

              G4Tubs* aWindTub = new G4Tubs( name2.str(),
                                             windparams.data()[0],
                                             windparams.data()[1],
                                             windparams.data()[2]+2.,// to satisfy a G4SubtractionSolid feature
                                             windparams.data()[3],
                                             windparams.data()[4]);

              CLHEP::HepRotation* windRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
              OR->getRotationFromOrientation( *windRotat, holeOrientsSA[hID] );


              if ( 0 == aSolid ) {
                aSolid = new G4SubtractionSolid( extShieldVol.name,
                                                 block,
                                                 aWindTub,
                                                 windRotat,
                                                 holeLocSA[thisHID]);
              } else {
                G4SubtractionSolid * bSolid = new G4SubtractionSolid
                  ( extShieldVol.name,
                    aSolid,
                    aWindTub,
                    windRotat,
                    holeLocSA[thisHID]);
                aSolid = bSolid;
              }
            } // End of loop over holes. Now loop over notches.


            // Find the index of the first notch, for accessing info lists...
            int notchID = notchIDSA[i];

            for ( int jNotch = 0; jNotch < nNotchesSA[i]; jNotch++ ) {
              int thisNID = notchID + jNotch; //index for this notch

              // Put notch(es) into box now

              // Each notch is implemented as a
              // G4SubtractionSolid to allow for another volume placement
              // through it

              // Now make the notch
              name << "n" << jNotch+1;


              // Get dimensions of this box
              std::vector<double>tempDims = notchDimSA[thisNID];

              std::ostringstream name2;
              name2 << "SaddleBox" << i+1 << "Notch" << jNotch+1;

              G4Box* aNotchBox = new G4Box(  name2.str(),
                                             tempDims[0],
                                             tempDims[1],
                                             tempDims[2]);

              CLHEP::HepRotation* notchRotat(nullptr);


              if ( 0 == aSolid ) {
                aSolid = new G4SubtractionSolid( extShieldVol.name,
                                                 block,
                                                 aNotchBox,
                                                 notchRotat,
                                                 notchLocSA[thisNID]);
              } else {
                G4SubtractionSolid * bSolid = new G4SubtractionSolid
                  ( extShieldVol.name,
                    aSolid,
                    aNotchBox,
                    notchRotat,
                    notchLocSA[thisNID]);
                aSolid = bSolid;
              }
            } // End of loop over notches

            extShieldVol.solid = aSolid;

            finishNesting(extShieldVol,
                          findMaterialOrThrow(matsSA[i]),
                          itsSARotat,
                          extShieldVol.centerInParent,
                          parent.logical,
                          0,
                          saddleIsVisible,
                          G4Colour::Magenta(),
                          saddleIsSolid,
                          forceAuxEdgeVisible,
                          placePV,
                          doSurfaceCheck);
          } else {

            // Build each normal box here.  Normal means no holes or notches.

            VolumeInfo extShieldVol(name.str(),
                                    sitesSA[i]-parent.centerInMu2e(),
                                    parent.centerInWorld);

            extShieldVol.solid = new G4ExtrudedSolid( extShieldVol.name,
                                                      itsOutline,
                                                      hlen,
                                                      G4TwoVector(0,0), 1.,
                                                      G4TwoVector(0,0), 1.);


            finishNesting(extShieldVol,
                          findMaterialOrThrow(matsSA[i]),
                          itsSARotat,
                          extShieldVol.centerInParent,
                          parent.logical,
                          0,
                          saddleIsVisible,
                          G4Colour::Magenta(),
                          saddleIsSolid,
                          forceAuxEdgeVisible,
                          placePV,
                          doSurfaceCheck);

          } // end of if...else...

        } // end of for loop over saddles

    } else { // end of version 1, 2, 3.  At this time, other versions are not
      // implemented.
      std::cout << "Requested Saddle version " << saddleSet->getVersion() << ".  Only versions 1 and 2 currently supported." << std::endl;
    }
  } // end of constructSaddles fn

} // namespace mu2e
