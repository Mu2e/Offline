// ***********************************************************
// This is the definitions file for the OrientationResolver 
// class, which is a helper to turn geometry orientation 
// strings into CLHEP::HepRotations.  This is described in 
// more detail in the paper in docdb #4678 and the slides
// in docdb #4999.
// Original author:  David Norvil Brown
// University of Louisville
// Code written in Nov. 2014, moved to this file in May 2016.
// ***********************************************************


#include "GeneralUtilities/inc/OrientationResolver.hh"
#include <string>

namespace mu2e {

  void OrientationResolver::getRotationFromOrientation( CLHEP::HepRotation& aRotation,
							std::string orient ) {
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
    if ( orient == "400" ) { // Special 45 degree turn around x
      aRotation.rotateX(45.0*CLHEP::degree);
      return;
    }
    if ( orient == "500" ) { // Special 45 degree turn around x other way
      aRotation.rotateX(-45.0*CLHEP::degree);
      return;
    }
    if ( orient == "060" ) { // Special -32 degree turn around y for pipe in TS
      aRotation.rotateY(-32.0*CLHEP::degree);
      return;
    }
    if ( orient == "0b0" ) { // Special -(90+32) degree turn for bend in TS
      aRotation.rotateY(-122.0*CLHEP::degree);
      return;
    }
    if ( orient == "ll0" ) { // Special 78 and 2 degree turn for PS pump out pipe
      aRotation.rotateY(-78.0*CLHEP::degree);
      aRotation.rotateX(-2.0*CLHEP::degree);
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
} // End namespace mu2e
