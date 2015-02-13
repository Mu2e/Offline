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
	    std::cout << __func__ << " making " << name.str() << std::endl;

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
	    
	    CLHEP::HepRotation* windRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	    getRotationFromOrientation( *windRotat, holeOrientsDS[hID] );

	    //	std::cout << "DNB** Making extShieldVol" << std::endl;

	    //	  std::cout << "DNB** About to make SubtractionSolid" << std::endl;

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

	    std::cout << __func__ << " making " << name.str() << std::endl;


	    // Get dimensions of this box
	    std::vector<double>tempDims = notchDimDS[thisNID];

	    std::ostringstream name2;
	    name2 << "ExtShieldDownstreamBox" << i+1 << "Notch" << jNotch+1;

	    G4Box* aNotchBox = new G4Box(  name2.str(), 
					   tempDims[0],
					   tempDims[1],
					   tempDims[2]);

	    CLHEP::HepRotation* notchRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);

	    //	  std::cout << "DNB** About to make SubtractionSolid" << std::endl;

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
			config.getBool("ExtShieldDownstream.visible"),
			G4Colour::Magenta(),
			config.getBool("ExtShieldDownstream.solid"),
			forceAuxEdgeVisible,
			placePV,
			doSurfaceCheck);
	  
	}

      }

  }

}

