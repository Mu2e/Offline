// $Id: constructExternalShielding.cc,v 1.6 2014/09/19 19:15:02 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/19 19:15:02 $
// David Norvil Brown, University of Louisville, November 2014
//
// 

#include "Mu2eG4/inc/constructExternalShielding.hh"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

// Include each shield here...
#include "ExternalShieldingGeom/inc/ExtShieldUpstream.hh"
#include "ExternalShieldingGeom/inc/ExtShieldDownstream.hh"
// etc...
#include "GeometryService/inc/GeomHandle.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
//#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4TwoVector.hh"
#include "CLHEP/Vector/Rotation.h"
#include "G4NistManager.hh"

#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"

#include <vector>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <string>

namespace mu2e {

  void getRotationFromOrientation ( CLHEP::HepRotation& aRotation, 
				    std::string orient ){
    // This is a helper function for determining the correct HepRotation
    // starting from an orientation number, described in docdb #xxxx
    // Note:  the length of orient should always be 3, enforced in the
    // Maker functions.
    if ( orient == "000" ) return;  // Don't need to reorient
    int uRotNum=0,vRotNum=0,wRotNum=0;
    for ( int iChar = 0; iChar < 3; iChar++ ) {
      // unpack the orientation number
      int genRotNum = 0; // Default to no rotation if mis-specified.  Throw instead, maybe?
      char digit = orient[iChar];
      if ( digit == '1' ) genRotNum = 1;
      if ( digit == '2' ) genRotNum = 2;
      if ( digit == '3' ) genRotNum = 3;

      if ( iChar == 0 ) uRotNum = genRotNum;
      if ( iChar == 1 ) vRotNum = genRotNum;
      if ( iChar == 2 ) wRotNum = genRotNum;
    }
//     int wRotNum  = orient%4;
//     orient  -= wRotNum;
//     int vRotNum  = orient%16;
//     orient  -= vRotNum;
//     int uRotNum  = orient/16;
//     vRotNum      = vRotNum/4;
//    std::cout << "DNB:  uRotNum, vRotNum, wRotNum = " << uRotNum << ", " << vRotNum << ", " << wRotNum << std::endl;
    double phi   = (double)uRotNum*90*CLHEP::degree;
    double theta = (double)vRotNum*90.0*CLHEP::degree;
    double psi   = (double)wRotNum*90.0*CLHEP::degree;
    //    std::cout << "DNB:  phi, theta, psi = " << phi << ", " << theta << ", " << psi << std::endl;  
    if ( wRotNum != 0 ) aRotation.rotateZ(-psi);
    // Now have to accommodate the change of axes based on Z rotation
    // If 90 rot in Z, rotation around y becomes rot around neg x
    // If 180 in Z, rotation around y becomes rot around neg y
    // If 270 in Z, rot around y becomes rot around x
    if ( wRotNum == 0 ) {
      if ( vRotNum != 0 ) aRotation.rotateY(-theta);
      if ( vRotNum == 0 ) {
	if ( uRotNum != 0 ) aRotation.rotateX(-phi);
      } else if ( vRotNum == 1 ) {
	if ( uRotNum != 0 ) aRotation.rotateZ(-phi);
      } else if ( vRotNum == 2 ) {
	if ( uRotNum != 0 ) aRotation.rotateX(phi);
      } else if ( vRotNum == 3 ) {
	if ( uRotNum != 0 ) aRotation.rotateZ(phi);
      }
    } else if ( wRotNum == 1 ) {
      // y-axis has moved to x and x to -y 
      if ( vRotNum != 0 ) aRotation.rotateX(theta);
      if ( vRotNum == 0 ) {
	if ( uRotNum != 0 ) aRotation.rotateY(-phi);
      } else if ( vRotNum == 1 ) {
	if ( uRotNum != 0 ) aRotation.rotateZ(-phi);
      } else if ( vRotNum == 2 ) {
	if ( uRotNum != 0 ) aRotation.rotateY(phi);
      } else if ( vRotNum == 3 ) {
	if ( uRotNum != 0 ) aRotation.rotateZ(phi);
      }
    } else if ( wRotNum == 2 ) {
      // y axis now -y and x axis now -x
      if ( vRotNum != 0 ) aRotation.rotateY(theta);
      if ( vRotNum == 0 ) {
	if ( uRotNum != 0 ) aRotation.rotateX(phi);
      } else if ( vRotNum == 1 ) {
	if ( uRotNum != 0 ) aRotation.rotateZ(-phi);
      } else if ( vRotNum == 2 ) {
	if ( uRotNum != 0 ) aRotation.rotateX(-phi);
      } else if ( vRotNum == 3 ) {
	if ( uRotNum != 0 ) aRotation.rotateZ(phi);
      }
    } else if ( wRotNum == 3 ) {
      // y axis is now -x and x axis is now +y
      if ( vRotNum != 0 ) aRotation.rotateX(-theta);
      if ( vRotNum == 0 ) {
	if ( uRotNum != 0 ) aRotation.rotateY(phi);
      } else if ( vRotNum == 1 ) {
	if ( uRotNum != 0 ) aRotation.rotateZ(-phi);
      } else if ( vRotNum == 2 ) {
	if ( uRotNum != 0 ) aRotation.rotateY(-phi);
      } else if ( vRotNum == 3 ) {
	if ( uRotNum != 0 ) aRotation.rotateZ(phi);
      }
    }
    return;
  }
  // End of utility function for converting orientation to rotation
 

  //================================================================

  void constructExternalShielding(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<ExtShieldUpstream> extshldUp;
    GeomHandle<ExtShieldDownstream> extshldDn;


    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

    
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
	std::vector<double> lwhs = dims[i];
	std::vector<double> dlwhs = tols[i];
	for ( unsigned int idim = 0; idim < lwhs.size(); idim++ ) {
	  lwhs[idim] += dlwhs[idim]/2.0;
	}

	//  Make the name of the box
	std::ostringstream name;
	name << "ExtShieldUpstreamBox_" << i+1 ;

	// Make the needed rotation by parsing orientation
        CLHEP::HepRotation* itsRotat= new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	std::string orientInit = orients[i];
	//	std::cout << "DNB:  orientInit is: " << orientInit << std::endl;
	getRotationFromOrientation(*itsRotat, orientInit);

	// Build each box here

	nestBox( name.str(), lwhs, findMaterialOrThrow(mats[i]),
		 itsRotat, sites[i]-parent.centerInMu2e(),
		 parent.logical,
		 0,
		 config.getBool("ExtShieldUpstream.visible"),
		 G4Colour::Magenta(),
		 config.getBool("ExtShieldUpstream.solid"),
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
    std::vector<bool> hasHoleDS              = extshldDn->getHasHole();
    std::vector<bool> hasNotchDS             = extshldDn->getHasNotch();
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

	for ( unsigned int idim = 0; idim < vertices.size(); idim++ ) {
	  vertices[idim][0] += du;
	  vertices[idim][1] += dv;
	  G4TwoVector vertex( vertices[idim][0], vertices[idim][1] );
	  itsOutline.push_back(vertex);
	}
	hlen += dw;

	//  Make the name of the box
	std::ostringstream name;
	name << "ExtShieldDownstreamBox_" << i+1 ;

	// Make the needed rotation by parsing orientation
	std::string orientDSInit = orientsDS[i];
	//	std::cout << "DNB:  orientDSInit is: " << orientDSInit << std::endl;
	CLHEP::HepRotation* itsDSRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	getRotationFromOrientation( *itsDSRotat, orientDSInit );

	if ( hasHoleDS[i] ) {

	  // Build each box here, with window.

 	  // This box has a window.  implemented as a
	  // G4SubtractionSolid to allow for another volume placement
	  // through it

  	  // Now make the window (AKA "Hole")
  	  name << "window";

	  std::cout << __func__ << " making " << name.str() << std::endl;

	  std::ostringstream name1;
	  name1 << "ExtShieldDownstreamBox_sub_" << i+1;


	  G4ExtrudedSolid* block = new G4ExtrudedSolid( name1.str(),
							itsOutline,
							hlen,
							G4TwoVector(0,0), 1.,
							G4TwoVector(0,0), 1.);


  	  // Find the index of this hole, for accessing info lists...
  	  int hID = holeIDDS[i];

	  const TubsParams windparams(0.0,  //inner radius
				      holeRadDS[hID], // outer
				      holeLenDS[hID]  //obvious?
				      );

	  std::ostringstream name2;
	  name2 << "ExtShieldDownstreamBox" << i+1 << "window";

          G4Tubs* aWindTub = new G4Tubs( name2.str(), 
                                         windparams.data()[0], 
                                         windparams.data()[1], 
                                         windparams.data()[2]+2.,// to satisfy a G4SubtractionSolid feature
                                         windparams.data()[3], 
                                         windparams.data()[4]);

	CLHEP::HepRotation* windRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	getRotationFromOrientation( *windRotat, holeOrientsDS[hID] );

	//	std::cout << "DNB** Making extShieldVol" << std::endl;

	  VolumeInfo extShieldVol(name.str(),
				  sitesDS[i]-parent.centerInMu2e(),
				  parent.centerInWorld);

	  //	  std::cout << "DNB** About to make SubtractionSolid" << std::endl;

	  extShieldVol.solid = new G4SubtractionSolid( extShieldVol.name,
						       block,
						       aWindTub,
						       windRotat,
						       holeLocDS[hID]);

	  //	  std::cout << "DNB** About to finish nesting" << std::endl;

	  finishNesting(extShieldVol,
			findMaterialOrThrow(matsDS[i]),
			itsDSRotat,
			extShieldVol.centerInParent,
			parent.logical,
			0,
			config.getBool("ExtShieldDownstream.visible"),
			G4Colour::Magenta(),
			config.getBool("ExtShieldDownstream.solid"),
			forceAuxEdgeVisible,
			placePV,
			doSurfaceCheck);

	} else if ( hasNotchDS[i] ) {

	  // Build each box here, with Notch.

 	  // This box has a notch.  implemented as a
	  // G4SubtractionSolid to allow for another volume placement
	  // through it

  	  // Now make the notch 
  	  name << "notch";

	  std::cout << __func__ << " making " << name.str() << std::endl;

	  std::ostringstream name1;
	  name1 << "ExtShieldDownstreamBox_sub_" << i+1;


	  G4ExtrudedSolid* block = new G4ExtrudedSolid( name1.str(),
							itsOutline,
							hlen,
							G4TwoVector(0,0), 1.,
							G4TwoVector(0,0), 1.);


  	  // Find the index of this hole, for accessing info lists...
  	  int hID = notchIDDS[i];

	  // Get dimensions of this box
	  std::vector<double>tempDims = notchDimDS[hID];

	  std::ostringstream name2;
	  name2 << "ExtShieldDownstreamBox" << i+1 << "Notch";

          G4Box* aNotchBox = new G4Box(  name2.str(), 
                                         tempDims[0],
					 tempDims[1],
                                         tempDims[2]);

	  CLHEP::HepRotation* notchRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);

	  VolumeInfo extShieldVol(name.str(),
				  sitesDS[i]-parent.centerInMu2e(),
				  parent.centerInWorld);

	  //	  std::cout << "DNB** About to make SubtractionSolid" << std::endl;

	  extShieldVol.solid = new G4SubtractionSolid( extShieldVol.name,
						       block,
						       aNotchBox,
						       notchRotat,
						       notchLocDS[hID]);

	  //	  std::cout << "DNB** About to finish nesting" << std::endl;

	  finishNesting(extShieldVol,
			findMaterialOrThrow(matsDS[i]),
			itsDSRotat,
			extShieldVol.centerInParent,
			parent.logical,
			0,
			config.getBool("ExtShieldDownstream.visible"),
			G4Colour::Magenta(),
			config.getBool("ExtShieldDownstream.solid"),
			forceAuxEdgeVisible,
			placePV,
			doSurfaceCheck);

	} else {

	  // Build each box here

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
			config.getBool("ExtShieldDownstream.visible"),
			G4Colour::Magenta(),
			config.getBool("ExtShieldDownstream.solid"),
			forceAuxEdgeVisible,
			placePV,
			doSurfaceCheck);
	  
	}

      }

//     //----------------------------------------------------------------
//     // DS-Cend boxes <=====
//     //----------------------------------------------------------------

//     std::vector<std::vector<double>> dimes = enscendb->getDimensions();
//     std::vector<std::string> mates = enscendb->materialNames();
//     std::vector<CLHEP::Hep3Vector> sitees = enscendb->centersOfBoxes();
    
//     nBox = dimes.size();

//     for(int i = 0; i < nBox; i++)
//       {

// 	std::vector<double> lwhs = dimes[i];

// 	//  Make the name of the box
// 	std::ostringstream name;
// 	name << "ExtShieldCendBox_" << i+1 ;

// 	// Make a non-rotating rotation
// 	static CLHEP::HepRotation fakeRotat(CLHEP::HepRotation::IDENTITY);

// 	// Build each box here


//  	  const TubsParams windparams(0.0,  //inner radius
//  				      enscendb->holeRadius(hID), // outer
//  				      enscendb->holeHalfLength(hID)  //obvious?
//  				      );

//           std::ostringstream name2;
//           name2 << "ExtShieldCendBox_sub_" << i+1 << "window";

//           G4Tubs* awindTub = new G4Tubs( name2.str(), 
//                                          windparams.data()[0], 
//                                          windparams.data()[1], 
//                                          windparams.data()[2]+2.,// to satisfy a G4SubtractionSolid feature
//                                          windparams.data()[3], 
//                                          windparams.data()[4]);

// 	  VolumeInfo awindBox;
//           awindBox.name = name.str();
          
//           // we need to put the window on the DS axis

//           // fixme, this is specific to the Box #4 (but windows in any
//           // other boxes will probably never be needed anyway)

//           GeomHandle<DetectorSolenoid> ds;
//           G4ThreeVector const & dsP ( ds->position() );

//           G4ThreeVector offsetWRTDS(sitees[i].x()-dsP.x(), sitees[i].y()-dsP.y(), 0.0);

//           awindBox.solid = new G4SubtractionSolid(awindBox.name,awindBox1,awindTub,0,-offsetWRTDS);

//           finishNesting(awindBox,
//                         findMaterialOrThrow(mates[i]),
//                         0,
//                         sitees[i]-parent.centerInMu2e(),
//                         parent.logical,
//                         0,
//                         config.getBool("ExtShieldCendBoxes.visible"),
//                         G4Colour::Magenta(),
//                         config.getBool("ExtShieldCendBoxes.solid"),
//                         forceAuxEdgeVisible,
//                         placePV,
//                         doSurfaceCheck);

// 	} else {

// 	  // Just put the box in the world, no window
// 	  nestBox( name.str(), lwhs, findMaterialOrThrow(mates[i]),
// 		   &fakeRotat, sitees[i]-parent.centerInMu2e(),
// 		   parent.logical,
// 		   0,
// 		   config.getBool("ExtShieldCendBoxes.visible"),
// 		   G4Colour::Magenta(),
// 		   config.getBool("ExtShieldCendBoxes.solid"),
// 		   forceAuxEdgeVisible,
// 		   placePV,
// 		   doSurfaceCheck);

// 	}

//       }
    
  }

}

