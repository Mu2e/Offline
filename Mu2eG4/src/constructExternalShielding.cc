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
#include "CLHEP/Units/SystemOfUnits.h"
#include "cetlib/exception.h"

// Include each shield here...
#include "ExternalShieldingGeom/inc/ExtShieldUpstream.hh"
#include "ExternalShieldingGeom/inc/ExtShieldDownstream.hh"
#include "ExternalShieldingGeom/inc/Saddle.hh"
#include "ServicesGeom/inc/Pipe.hh"
#include "ServicesGeom/inc/ElectronicRack.hh"

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
    if ( orient == "550" ) {  // Special for proton beamline elements
      aRotation.rotateY(-13.62*CLHEP::degree);
      aRotation.rotateX(2.71*CLHEP::degree);
      return;
    }
    if ( orient == "040" ) { // Special 45 degree turn around y
      aRotation.rotateY(45.0*CLHEP::degree);
      return;
    }

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

    double phi   = (double)uRotNum*90*CLHEP::degree;
    double theta = (double)vRotNum*90.0*CLHEP::degree;
    double psi   = (double)wRotNum*90.0*CLHEP::degree;

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
    GeomHandle<Saddle> saddleSet;
    GeomHandle<Pipe>   pipeSet;
    GeomHandle<ElectronicRack>   rackSet;

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
	name << "ExtShieldDownstreamBox_" << i+1 ;

	// Make the needed rotation by parsing orientation
	std::string orientDSInit = orientsDS[i];

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
	    //	    std::cout << __func__ << " making " << name.str() << std::endl;

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

	    ///	    std::cout << __func__ << " making " << name.str() << std::endl;


	    // Get dimensions of this box
	    std::vector<double>tempDims = notchDimDS[thisNID];

	    std::ostringstream name2;
	    name2 << "ExtShieldDownstreamBox" << i+1 << "Notch" << jNotch+1;

	    G4Box* aNotchBox = new G4Box(  name2.str(), 
					   tempDims[0],
					   tempDims[1],
					   tempDims[2]);

	    CLHEP::HepRotation* notchRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);


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
	  
	} // end of if...else...

      } // end of for loop over boxes

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

    nBox = outlSA.size();

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

	CLHEP::HepRotation* itsSARotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	getRotationFromOrientation( *itsSARotat, orientSAInit );

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
	    //	    std::cout << __func__ << " making " << name.str() << std::endl;

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
	    
	    CLHEP::HepRotation* windRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	    getRotationFromOrientation( *windRotat, holeOrientsSA[hID] );


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

	    //	    std::cout << __func__ << " making " << name.str() << std::endl;


	    // Get dimensions of this box
	    std::vector<double>tempDims = notchDimSA[thisNID];

	    std::ostringstream name2;
	    name2 << "SaddleBox" << i+1 << "Notch" << jNotch+1;

	    G4Box* aNotchBox = new G4Box(  name2.str(), 
					   tempDims[0],
					   tempDims[1],
					   tempDims[2]);

	    CLHEP::HepRotation* notchRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);


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
			config.getBool("Saddle.visible"),
			G4Colour::Magenta(),
			config.getBool("Saddle.solid"),
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
			config.getBool("Saddle.visible"),
			G4Colour::Magenta(),
			config.getBool("Saddle.solid"),
			forceAuxEdgeVisible,
			placePV,
			doSurfaceCheck);
	  
	} // end of if...else...

      } // end of for loop over saddles

    // *******************************************************
    // ==> Make Pipes <========
    // *******************************************************

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

    for ( unsigned int it = 0; it < nPipes.size(); it++ ) {

      int nPipe = nPipes[it];
      int nComp = nComps[it];
      double len = pLeng[it];
      std::string flav = pFlav[it];
      if ( flav != "straight" ) {
	// Bow out not-so-gracefully
	throw cet::exception("GEOM") << " in constructExternalShielding, have not yet implemented non-straight segments of pipe."<<"\n";
      } // end of if not straight
      std::string fillMat = pFill[it];

      for ( int ip = 0; ip < nPipe; ip++ ) {

	// Make the container ("mother") volume for the type.
	TubsParams  pipeParams(0.0,cOutRad[it][0],len/2.0);
	std::ostringstream motherTubeName;
	motherTubeName << "pipeType" << it+1 << "Pipe" << ip+1;
	//	G4Tubs* motherSolidVol = new G4Tubs(motherTubeName.str(),0.0,
	//					    cOutRad[it][0],len/2.0,
	//					    0.0,CLHEP::twopi );
	std::ostringstream motherName;
	motherName << "pipeLogicVolType" << it+1 << "Pipe" << ip+1;

	CLHEP::HepRotation* pipeRotat = new 
	  CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	getRotationFromOrientation( *pipeRotat, pOrient[it][ip] );

	CLHEP::Hep3Vector pipePosInMu2e = pCent[it][ip] 
	  - parent.centerInMu2e();

	VolumeInfo motherLogVol = nestTubs(motherName.str(),
					   pipeParams,
					   findMaterialOrThrow(fillMat),
					   pipeRotat,
					   pipePosInMu2e,
					   parent.logical,
					   ip,
					   config.getBool("Pipe.visible"),
					   G4Colour::Magenta(),
					   config.getBool("Pipe.solid"),
					   forceAuxEdgeVisible,
					   placePV,
					   doSurfaceCheck);

	for ( int ic = 0; ic < nComp; ic++ ) {
	  // Create a tube for each "component" pipe and embed in the mother
	  // Logical volume
	  TubsParams componentParams(cInRad[it][ic],cOutRad[it][ic],len/2.0);
	  std::ostringstream compName;
	  compName << "pipeType" << it+1 << "Pipe" << ip+1 << 
	    "Component" << ic+1;
	  // Until we place it, the mother should be centered at 0,0,0
	  // So place relative to that using u- and v-components and w=0.
	  CLHEP::Hep3Vector posComp(cUOff[it][ic],cVOff[it][ic],0.0); 
	  nestTubs( compName.str(),  componentParams, 
		    findMaterialOrThrow(cMats[it][ic]),
		    0, posComp, motherLogVol,
		    ic,
		    config.getBool("Pipe.visible"),
		    G4Colour::Magenta(),
		    config.getBool("Pipe.solid"),
		    forceAuxEdgeVisible,
		    placePV,
		    doSurfaceCheck);

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

	getRotationFromOrientation(*itsRotatER, orientInitER);

	// Build each box here

	nestBox( nameER.str(), lwhsER, findMaterialOrThrow(matsER[i]),
		 itsRotatER, sitesER[i]-parent.centerInMu2e(),
		 parent.logical,
		 0,
		 config.getBool("ElectronicRack.visible"),
		 G4Colour::Magenta(),
		 config.getBool("ElectronicRack.solid"),
		 forceAuxEdgeVisible,
		 placePV,
		 doSurfaceCheck);

      } // end loop over ElectronicRack boxes


  } // end of constructExternalShielding fn

} // namespace mu2e

