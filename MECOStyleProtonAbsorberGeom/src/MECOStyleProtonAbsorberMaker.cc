//
// Construct and return MECOStyleProtonAbsorber
//
//
// $Id: MECOStyleProtonAbsorberMaker.cc,v 1.3 2012/08/09 22:20:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/08/09 22:20:16 $
//
// Original author MyeongJae Lee
//
// All in mu2e coordinate system, except target.z0 which is given by old Detector coordinate system. (-3904, 0, 12000)
// Geometry is based on old construcProtonAbsorber codes.
// distFromTargetEnd value is prepared for forward compatibility. Do not use now.

#include <iostream>
#include <iomanip>
#include <cmath>


// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorberMaker.hh"
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  MECOStyleProtonAbsorberMaker::MECOStyleProtonAbsorberMaker( SimpleConfig const& _config)
  {

    BuildIt (_config);

    int verbosity(_config.getInt("protonabsorber.verbosity",0));
    if ( verbosity > 0 ) PrintConfig();
  }

  void MECOStyleProtonAbsorberMaker::BuildIt ( SimpleConfig const& _config)
  {

    //////////////////////////////////////
    // All variables read from geom_01.txt
    //////////////////////////////////////

    // geometry and material
    double r1out0        = _config.getDouble("protonabsorber.OutRadius0", 335.2);
    double r2out1        = _config.getDouble("protonabsorber.OutRadius1", 380.0);
    double pabsZHalfLen      = _config.getDouble("protonabsorber.halfLength", 1250.0);
    double thick             = _config.getDouble("protonabsorber.thickness", 0.5);
    double distFromTargetEnd = 0.;
    // double distFromTargetEnd = _config.getDouble("protonabsorber.distFromTargetEnd", 0);
    std::string materialName = _config.getString("protonabsorber.materialName", "Polyethylene092");
    if (distFromTargetEnd > pabsZHalfLen *2. || distFromTargetEnd <0) {
      std::cerr <<"MECOStyleProtonAbsorberMaker: Wrong value for distFromTargetEnd. Set to 0" << std::endl;
      distFromTargetEnd = 0;
    }

    // TS and DS geometry for locating the center of Proton Absorber
    double rTorus = _config.getDouble("toyTS.rTorus");
    double ts5HalfLength = _config.getDouble("toyTS5.halfLength");
    double ds2HalfLength = _config.getDouble("toyDS2.halfLength");
    double solenoidOffset = _config.getDouble("mu2e.solenoidOffset");

    // Info on targets to locate the center of Proton Absorber
    vector<double> targetRadius;  _config.getVectorDouble("target.radii", targetRadius);
    double foilwid=_config.getDouble("target.deltaZ");
    double z0valt =_config.getDouble("target.z0");

    // adding virtual detector before and after target
    double vdHL = 0.;
    if (_config.getBool("hasVirtualDetector", true) ) {
      vdHL = _config.getDouble("vd.halfLength");
    }

    // subtract virtual detector thickness from the larger outer
    // radius of the proton absorber
    r2out1 -= 2.*vdHL;



    /////////////////////////////
    // Geometry related variables
    /////////////////////////////

    // DS2 offset
    double z0DSup = rTorus + 2.*ts5HalfLength + ds2HalfLength;

    // number of target foil
    double numoftf = (targetRadius.size()-1.0)*0.5;

    // we add space for the virtual detector here, and get the half length of target
    double taghalflen =(foilwid*numoftf) + 5.0 + 2.*vdHL; // what/why is 5.0 hardcoded here?

    // target offset in mu2e coordinate
    double tagoff =z0valt + 12000.0;

    // z of target end in mu2e coordinate
    double targetEnd = tagoff + taghalflen;

    // distance from target end to ds2-ds3 boundary
    double targetEndToDS2End = ds2HalfLength + (z0DSup - targetEnd);



    //////////////////////////////////////
    // Decide which pabs will be turned on
    //////////////////////////////////////

    bool pabs1 = true;
    bool pabs2 = true;

    // if pabs starts from DS3 region
    if (distFromTargetEnd > targetEndToDS2End) {
      std::cerr <<"MECOStyleProtonAbsorberMaker: pabs1 turned off." << std::endl;
      pabs1 = false;
    }

    // if pabs is short enouhg to locate at DS2 region only
    if (distFromTargetEnd + pabsZHalfLen*2.< targetEndToDS2End) {
      std::cerr <<"MECOStyleProtonAbsorberMaker: pabs2 turned off." << std::endl;
      pabs2 = false;
    }
    if (!pabs1 && !pabs2) {
      std::cerr <<"MECOStyleProtonAbsorberMaker: no pabs can be built." << std::endl;
      return;
    }


    //////////////
    // Half length
    //////////////

    // old style codes for pabs1/2 half length
    // double pabs1halflen = (ds2HalfLength + (z0DSup - targetEnd))*0.5;
    // double pabs2halflen  = pabsZHalfLen - pabs1halflen;

    double pabs1halflen = 0, pabs2halflen = 0;
    if (pabs1) {
      pabs1halflen = (targetEndToDS2End - distFromTargetEnd)*0.5;
      if (pabs1halflen > pabsZHalfLen ) pabs1halflen = pabsZHalfLen;
    }
    if (pabs2) {
      pabs2halflen = pabsZHalfLen - pabs1halflen;
    }

    /////////
    // Offset
    /////////

    // old style for pabs1/2 offset
    // double pabs1ZOffset = targetEnd + pabs1halflen ;
    // double pabs2ZOffset = targetEnd + pabs1halflen * 2. + pabs2halflen;

    double pabs1ZOffset = targetEnd + 0.01;
    double pabs2ZOffset = targetEnd + targetEndToDS2End + 0.01;
    if (pabs1) {
      pabs1ZOffset = targetEnd + distFromTargetEnd + pabs1halflen;
    }
    if (pabs2) {
      if (pabs1) {
        pabs2ZOffset = targetEnd + targetEndToDS2End + pabs2halflen;
      }
      else {
        pabs2ZOffset = targetEnd + distFromTargetEnd + pabs2halflen;
      }
    }

    /////////
    // Radius
    /////////

    // old style codes
    // interpolating the outer radius of the DS2 part
    // double r1out1 = ((r2out1 - r1out0)*(pabs1halflen/pabsZHalfLen)) + r1out0;
    // double r1in1  = r1out1 - thick;
    // double r1in0  = r1out0 - thick;
    // double r2in1  = r2out1 - thick;

    double pabs1rIn0 = 0, pabs1rIn1 = 0, pabs1rOut0 = 0, pabs1rOut1 = 0;
    double pabs2rIn0 = 0, pabs2rIn1 = 0, pabs2rOut0 = 0, pabs2rOut1 = 0;

    if (pabs1) {
      pabs1rOut0 = (r2out1 - r1out0)/(pabsZHalfLen*2.)*(pabs1ZOffset-pabs1halflen-targetEnd) + r1out0;
      pabs1rOut1 = (r2out1 - r1out0)/(pabsZHalfLen*2.)*(pabs1ZOffset+pabs1halflen-targetEnd) + r1out0;
      pabs1rIn0  = pabs1rOut0 - thick;
      pabs1rIn1  = pabs1rOut1 - thick;
    }
    if (pabs2) {
      pabs2rOut0 = (r2out1 - r1out0)/(pabsZHalfLen*2.)*(pabs2ZOffset-pabs2halflen-targetEnd) + r1out0;
      pabs2rOut1 = (r2out1 - r1out0)/(pabsZHalfLen*2.)*(pabs2ZOffset+pabs2halflen-targetEnd) + r1out0;
      pabs2rIn0  = pabs2rOut0 - thick;
      pabs2rIn1  = pabs2rOut1 - thick;
    }

    /////////
    // Build
    /////////

    _pabs = auto_ptr<MECOStyleProtonAbsorber>(new MECOStyleProtonAbsorber());

    CLHEP::Hep3Vector pabs1Offset(-1.*solenoidOffset, 0.0, pabs1ZOffset);
    CLHEP::Hep3Vector pabs2Offset(-1.*solenoidOffset, 0.0, pabs2ZOffset);

    _pabs->_parts.push_back( MECOStyleProtonAbsorberPart( 0, pabs1Offset, pabs1rOut0, pabs1rIn0, pabs1rOut1, pabs1rIn1, pabs1halflen, materialName));
    _pabs->_parts.push_back( MECOStyleProtonAbsorberPart( 1, pabs2Offset, pabs2rOut0, pabs2rIn0, pabs2rOut1, pabs2rIn1, pabs2halflen, materialName));

    // global variables
    (_pabs->_pabs1flag) = pabs1;
    (_pabs->_pabs2flag) = pabs2;
    (_pabs->_materialName) = materialName;
    (_pabs->_vdHL) = vdHL;
    (_pabs->_ds2zcenter) = z0DSup;
    (_pabs->_distfromtargetend) = distFromTargetEnd;
    (_pabs->_halflength) = pabsZHalfLen;
    (_pabs->_thickness) = thick;


  }

  void MECOStyleProtonAbsorberMaker::PrintConfig ( ) {

    double pabs1z = (_pabs->part(0)).center().z();
    double pabs2z = (_pabs->part(1)).center().z();
    double pabs1hl = (_pabs->part(0)).halfLength();
    double pabs2hl = (_pabs->part(1)).halfLength();
    std::cout<<"MECOStyleProtonAbsorberMaker Configuration -----------------"<<std::endl;
    std::cout<<" Dist from target end : " << _pabs->distanceFromTargetEnd() << std::endl;
    std::cout<<" vdHL : " << _pabs->virtualDetectorHalfLength() << std::endl;
    std::cout<<" pabs1len (full length) : " << pabs1hl*2. << std::endl;
    std::cout<<" pabs2len (full length) : " << pabs2hl*2. << std::endl;
    std::cout<<" pabs1 extent in Mu2e : " << pabs1z - pabs1hl << ", " << pabs1z << ", " << pabs1z+pabs1hl << std::endl;
    std::cout<<" pabs2 extent in Mu2e : " << pabs2z - pabs2hl << ", " << pabs2z << ", " << pabs2z+pabs2hl << std::endl;
    std::cout<<" [rIn, rOut] : " << std::endl;
    std::cout<<" =  [" << _pabs->part(0).innerRadiusAtStart() << ", " << _pabs->part(0).outerRadiusAtStart() <<"], " << std::endl
             <<"    [" << _pabs->part(0).innerRadiusAtEnd()   << ", " << _pabs->part(0).outerRadiusAtEnd()   <<"], " << std::endl
             <<"    [" << _pabs->part(1).innerRadiusAtStart() << ", " << _pabs->part(1).outerRadiusAtStart() <<"], " << std::endl
             <<"    [" << _pabs->part(1).innerRadiusAtEnd()   << ", " << _pabs->part(1).outerRadiusAtEnd()   <<"]"   << std::endl;
    std::cout<<" Material : " << _pabs->fillMaterial() <<  std::endl;
    std::cout<<" pabs1 is " << ( (_pabs->isAvailable(0)) ? "valid" : "invalid" ) << ", " ;
    std::cout<<" pabs2 is " << ( (_pabs->isAvailable(1)) ? "valid" : "invalid" ) << std::endl ;
  }



} // namespace mu2e
