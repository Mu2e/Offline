// David Norvil Brown, University of Louisville, November 2014
//
//

#include "Mu2eG4/inc/constructExternalShielding.hh"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "cetlib_except/exception.h"

// Include each shield here...
#include "ExternalShieldingGeom/inc/ExtShieldUpstream.hh"
#include "ExternalShieldingGeom/inc/ExtShieldDownstream.hh"

// etc...
#include "GeometryService/inc/GeomHandle.hh"
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


  void constructExternalShielding(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<ExtShieldUpstream> extshldUp;
    GeomHandle<ExtShieldDownstream> extshldDn;

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    // Utility for converting orientations to rotations
    OrientationResolver* OR = new OrientationResolver();

    // Get config info
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "ExtShield",           "ExtShield");
    geomOptions->loadEntry( config, "ExtShielUpstream",    "ExtShieldUpstream");
    geomOptions->loadEntry( config, "ExtShieldDownstream", "ExtShieldDownstream");

    const bool isExtShieldUpstreamVisible   = geomOptions->isVisible("ExtShieldUpstream");
    const bool isExtShieldUpstreamSolid     = geomOptions->isSolid("ExtShieldUpstream");
    const bool isExtShieldDownstreamVisible = geomOptions->isVisible("ExtShieldDownstream");
    const bool isExtShieldDownstreamSolid   = geomOptions->isSolid("ExtShieldDownstream");
    const bool forceAuxEdgeVisible          = geomOptions->forceAuxEdgeVisible("ExtShield");
    const bool doSurfaceCheck               = geomOptions->doSurfaceCheck("ExtShield");
    const bool placePV                      = geomOptions->placePV("ExtShield");


    //----------------------------------------------------------------
    // Upstream External Shielding Boxes <=====
    //----------------------------------------------------------------
    // Fill the vectors of information about the boxes making up the Upstream
    // shields
    std::vector<std::vector<double> > dims = extshldUp->getBoxDimensions();
    std::vector<std::vector<double> > tols = extshldUp->getBoxTolerances();
    std::vector<std::string> mats = extshldUp->getMaterialNames();
    std::vector<CLHEP::Hep3Vector> sites = extshldUp->getCentersOfBoxes();
    std::vector<std::string> orients = extshldUp->getOrientations();

    int nBox = dims.size();

    for(int i = 0; i < nBox; i++)
    {

        // combine the tolerances with the dimensions.
        std::vector<double> lwhs  = dims[i];
        std::vector<double> dlwhs = tols[i];
        for ( unsigned int idim = 0; idim < lwhs.size(); idim++ ) {
          lwhs[idim] += dlwhs[idim]/2.0;
        }

        //  Make the name of the box
        std::ostringstream name;
        name << "ExtShieldUpstreamBox_" << i+1 ;

        // Make the needed rotation by parsing orientation
        CLHEP::HepRotation* itsRotat= reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
        std::string orientInit = orients[i];

        OR->getRotationFromOrientation(*itsRotat, orientInit);

        // Build each box here
        nestBox( name.str(), lwhs, findMaterialOrThrow(mats[i]),
                 itsRotat, sites[i]-parent.centerInMu2e(),
                 parent.logical,
                 0,
                 isExtShieldUpstreamVisible,
                 G4Colour::Magenta(),
                 isExtShieldUpstreamSolid,
                 forceAuxEdgeVisible,
                 placePV,
                 doSurfaceCheck);
    }


    //----------------------------------------------------------------
    // Downstream External Extrusion Shielding Boxes <=====
    //----------------------------------------------------------------
    // Fill the vectors of information about the boxes making up the Downstream
    // shields
    std::vector<std::vector<std::vector<double> > > outlDS = extshldDn->getOutlines();
    std::vector<double> lengDS               = extshldDn->getLengths();
    std::vector<std::vector<double> > tolsDS = extshldDn->getTolerances();
    std::vector<std::string> matsDS          = extshldDn->getMaterialNames();
    std::vector<CLHEP::Hep3Vector> sitesDS   = extshldDn->getCentersOfBoxes();
    std::vector<std::string> orientsDS       = extshldDn->getOrientations();
    std::vector<int> nHolesDS                = extshldDn->getNHoles();
    std::vector<int> nNotchesDS              = extshldDn->getNNotches();
    std::vector<int> holeIDDS                = extshldDn->getHoleIndices();
    std::vector<CLHEP::Hep3Vector> holeLocDS = extshldDn->getHoleLocations();
    std::vector<double> holeRadDS            = extshldDn->getHoleRadii();
    std::vector<double> holeLenDS            = extshldDn->getHoleLengths();
    std::vector<std::string> holeOrientsDS   = extshldDn->getHoleOrientations();
    std::vector<int> notchIDDS               = extshldDn->getNotchIndices();
    std::vector<CLHEP::Hep3Vector> notchLocDS= extshldDn->getNotchLocations();
    std::vector<std::vector<double> > notchDimDS = extshldDn->getNotchDimensions();

    nBox = outlDS.size();

    for(int i = 0; i < nBox; i++)
      {
        // combine the tolerances with the outline dimensions.
        std::vector<G4TwoVector> itsOutline;
        std::vector<std::vector<double> > vertices = outlDS[i];
        double du = tolsDS[i][0];
        double dv = tolsDS[i][1];
        double dw = tolsDS[i][2];
        double hlen = lengDS[i];
        double xCnt[vertices.size()];
        double yCnt[vertices.size()];
        double epsilon = 1.0e-4;

        for ( unsigned int idim = 0; idim < vertices.size(); idim++ ) {
          // Decide how to apply tolerances by looking where other vertices
          // lie relative to this vertex.
          xCnt[idim] = 0.0;
          yCnt[idim] = 0.0;
          bool Q1, Q2, Q3, Q4;  // The four quadrants on the Cartesian plane
          Q1 = Q2 = Q3 = Q4 = false;
          bool N, S, E, W;      // and the cardinal directions
          N = S = E = W = false;
          for ( unsigned int j = 0; j < vertices.size(); j++ ) {
            if ( j != idim ) {
              // Here we figure out where other vertices are relative to
              // This vertex.
              double xDiff = vertices[j][0] - vertices[idim][0];
              double yDiff = vertices[j][1] - vertices[idim][1];
              if ( xDiff > epsilon ) { // Must be Q1 or Q4 or E
                if ( yDiff > epsilon ) {
                  Q1 = true;
                } else if ( yDiff < -epsilon ) {
                  Q4 = true;
                } else {
                  E = true;
                }
              } else if ( xDiff < -epsilon ) { // Must be Q2 or Q3 or W
                if ( yDiff > epsilon ) {
                  Q2 = true;
                } else if ( yDiff < -epsilon ) {
                  Q3 = true;
                } else {
                  W = true;
                }
              } else { // Must be North or South
                if ( yDiff > epsilon ) {
                  N = true;
                } else if ( yDiff < epsilon ) {
                  S = true;
                } // There should never be an else case since I test for j!=idim
              }
            } // endif j != idim
          } // end for loop over all other vertices
          // Now for idim, decide how to apply tols.

          if ( ( !Q1 && !Q4 && !E ) || ( Q2 && Q3 && W && !(Q1 && Q4 && E))
               || ( W && !E ) || ( N && E && !Q1 ) || (S && E && !Q4) ) {
            xCnt[idim] = 1.0;
          } else if ( ( !Q2 && !Q3 && !W )
                      || ( Q1 && Q4 && E && !(Q2 && Q3 && W)) || ( E && !W )
                      || ( N && W && !Q2 ) || ( S && W && !Q3 ) ) {
            xCnt[idim] = -1.0;
          }
          if ( ( !Q1 && !Q2 && !N ) || ( ((S && W) || (N && E)) && Q3 && !Q1 )
               || ( ((S && E) || (N && W)) && Q4 && !Q2) || ( S && !N ) ) {
            yCnt[idim] = 1.0;
          } else if ( ( !Q3 && !Q4 && !S ) // Along bottom edge
                      || ( ((S && W) || (N && E)) && Q1 && !Q3 ) // mid corner
                      || ( ((S && E) || (N && W)) && Q2 && !Q4 ) // int corner
                      || ( N && !S )) {
            yCnt[idim] = -1.0;
          }
          if ( fabs(xCnt[idim]) < epsilon ) { // Special vertex - check
            unsigned int neighbor1 = idim + 1;
            if ( idim == vertices.size() ) neighbor1 = 0;
            unsigned int neighbor2 = idim - 1;
            if ( idim == 0 ) neighbor2 = vertices.size() - 1;
            double xC = (vertices[neighbor1][0] + vertices[neighbor2][0])/2.0;
            if ( xC < vertices[idim][0] ) {
              xCnt[idim] = 1.0;
            } else if ( xC > vertices[idim][0] ) {
              xCnt[idim] = -1.0;
            }
          }
          if ( fabs(yCnt[idim]) < epsilon ) { // Special vertex - check
            unsigned int neighbor1 = idim + 1;
            if ( idim == vertices.size() ) neighbor1 = 0;
            unsigned int neighbor2 = idim - 1;
            if ( idim == 0 ) neighbor2 = vertices.size();
            double yC = (vertices[neighbor1][1] + vertices[neighbor2][1])/2.0;
            if ( yC < vertices[idim][1] ) {
              yCnt[idim] = 1.0;
            } else if ( yC > vertices[idim][1] ) {
              yCnt[idim] = -1.0;
            }
          }

        } // loop over all vertices to find xCnt and yCnt

        for ( unsigned int idim = 0; idim < vertices.size(); idim++ ) {
          vertices[idim][0] += xCnt[idim] * du;
          vertices[idim][1] += yCnt[idim] * dv;

          G4TwoVector vertex( vertices[idim][0], vertices[idim][1] );
          itsOutline.push_back(vertex);
        } // modified all the vertices appropriately
        hlen += dw;

        //  Make the name of the box
        std::ostringstream name;
        name << "ExtShieldDownstreamBox_" << i+1 ;

        // Make the needed rotation by parsing orientation
        std::string orientDSInit = orientsDS[i];

        CLHEP::HepRotation* itsDSRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
        OR->getRotationFromOrientation( *itsDSRotat, orientDSInit );

        // ****************************************
        // Now add holes and notches
        // ***************************************
        if ( nHolesDS[i] + nNotchesDS[i] > 0 ) {

          // This box has at least one window or notch.
          // Each is implemented as a
          // G4SubtractionSolid to allow for another volume placement
          // through it.  First go ahead and make the extrusion for the block.
          std::ostringstream name1;
          name1 << "ExtShieldDownstreamBox_sub_" << i+1;

          G4ExtrudedSolid* block = new G4ExtrudedSolid( name1.str(),
                                                        itsOutline,
                                                        hlen,
                                                        G4TwoVector(0,0), 1.,
                                                        G4TwoVector(0,0), 1.);

          // Now create the volume for the block
          VolumeInfo extShieldVol(name.str(),
                                  sitesDS[i]-parent.centerInMu2e(),
                                  parent.centerInWorld);


          // Create a subtraction solid that will become the solid for the
          // volume.

          G4SubtractionSolid* aSolid = 0;

          // Find the index of the first hole, for accessing info lists...
          int hID = holeIDDS[i];

          // Now loop over the holes and make them.
          for ( int jHole = 0; jHole < nHolesDS[i]; jHole++ ) {
            // Now make the window (AKA "Hole")
            name << "h" << jHole+1;
            //      std::cout << __func__ << " making " << name.str() << std::endl;

            int thisHID = hID + jHole; // get pointer for right hole

            const TubsParams windparams(0.0,  //inner radius
                                        holeRadDS[thisHID], // outer
                                        holeLenDS[thisHID]  //obvious?
                                        );

            std::ostringstream name2;
            name2 << "ExtShieldDownstreamBox" << i+1 << "hole" << jHole+1;

            G4Tubs* aWindTub = new G4Tubs( name2.str(),
                                           windparams.data()[0],
                                           windparams.data()[1],
                                           windparams.data()[2]+2.,// to satisfy a G4SubtractionSolid feature
                                           windparams.data()[3],
                                           windparams.data()[4]);

            CLHEP::HepRotation* windRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
            OR->getRotationFromOrientation( *windRotat, holeOrientsDS[hID] );


            if ( 0 == aSolid ) {
              aSolid = new G4SubtractionSolid( extShieldVol.name,
                                               block,
                                               aWindTub,
                                               windRotat,
                                               holeLocDS[thisHID]);
            } else {
              G4SubtractionSolid * bSolid = new G4SubtractionSolid
                ( extShieldVol.name,
                  aSolid,
                  aWindTub,
                  windRotat,
                  holeLocDS[thisHID]);
              aSolid = bSolid;
            }
          } // End of loop over holes. Now loop over notches.


          // Find the index of the first notch, for accessing info lists...
          int notchID = notchIDDS[i];

          for ( int jNotch = 0; jNotch < nNotchesDS[i]; jNotch++ ) {
            int thisNID = notchID + jNotch; //index for this notch

            // Put notch(es) into box now

            // Each notch is implemented as a
            // G4SubtractionSolid to allow for another volume placement
            // through it

            // Now make the notch
            name << "n" << jNotch+1;

            ///     std::cout << __func__ << " making " << name.str() << std::endl;


            // Get dimensions of this box
            std::vector<double>tempDims = notchDimDS[thisNID];

            std::ostringstream name2;
            name2 << "ExtShieldDownstreamBox" << i+1 << "Notch" << jNotch+1;

            G4Box* aNotchBox = new G4Box(  name2.str(),
                                           tempDims[0],
                                           tempDims[1],
                                           tempDims[2]);

            CLHEP::HepRotation* notchRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));


            if ( 0 == aSolid ) {
              aSolid = new G4SubtractionSolid( extShieldVol.name,
                                               block,
                                               aNotchBox,
                                               notchRotat,
                                               notchLocDS[thisNID]);
            } else {
              G4SubtractionSolid * bSolid = new G4SubtractionSolid
                ( extShieldVol.name,
                  aSolid,
                  aNotchBox,
                  notchRotat,
                  notchLocDS[thisNID]);
              aSolid = bSolid;
            }
          } // End of loop over notches

          extShieldVol.solid = aSolid;

          finishNesting(extShieldVol,
                        findMaterialOrThrow(matsDS[i]),
                        itsDSRotat,
                        extShieldVol.centerInParent,
                        parent.logical,
                        0,
                        isExtShieldDownstreamVisible,
                        G4Colour::Magenta(),
                        isExtShieldDownstreamSolid,
                        forceAuxEdgeVisible,
                        placePV,
                        doSurfaceCheck);
        } else {

          // Build each normal box here.  Normal means no holes or notches.

          VolumeInfo extShieldVol(name.str(),
                                  sitesDS[i]-parent.centerInMu2e(),
                                  parent.centerInWorld);

          extShieldVol.solid = new G4ExtrudedSolid( extShieldVol.name,
                                                    itsOutline,
                                                    hlen,
                                                    G4TwoVector(0,0), 1.,
                                                    G4TwoVector(0,0), 1.);


          finishNesting(extShieldVol,
                        findMaterialOrThrow(matsDS[i]),
                        itsDSRotat,
                        extShieldVol.centerInParent,
                        parent.logical,
                        0,
                        isExtShieldDownstreamVisible,
                        G4Colour::Magenta(),
                        isExtShieldDownstreamSolid,
                        forceAuxEdgeVisible,
                        placePV,
                        doSurfaceCheck);

        } // end of if...else...

      } // end of for loop over boxes
    // Do a bit of cleanup
    if ( 0 != OR ) {
      delete OR;
      OR = 0;
    }
  } // end of constructExternalShielding fn

} // namespace mu2e
