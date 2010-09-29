//
// Geometry and identifier info about the Calorimeter.
//
//
// $Id: Calorimeter.cc,v 1.5 2010/09/29 19:37:58 logash Exp $
// $Author: logash $ 
// $Date: 2010/09/29 19:37:58 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include "CalorimeterGeom/inc/Calorimeter.hh"

using namespace std;

namespace mu2e { 

  //
  // Convert coordinates from Mu2e frame to local frame of crystal, identified by roid
  //
  CLHEP::Hep3Vector Calorimeter::toCrystalFrame(int roid, 
						CLHEP::Hep3Vector const& pos) const {
    int vaneid = getVaneByRO(roid);
    const Vane & vane = getVane(vaneid);

    CLHEP::Hep3Vector vlocal(-_roHalfThickness,
			     (2*getCrystalRByRO(roid)-_nCrystalR+1)*_crystalHW,
			     (2*getCrystalZByRO(roid)-_nCrystalZ+1)*_crystalHW );

    return *(vane.getRotation())*(pos-vane.getOrigin())-vlocal;
  }

  //
  // Get crystal origin (center) in Mu2e coordinates
  //
  CLHEP::Hep3Vector Calorimeter::getCrystalOriginByRO(int roid) const {

    int vaneid = getVaneByRO(roid);
    const Vane & vane = getVane(vaneid);

    // Crystal center in vane coordinates
    CLHEP::Hep3Vector vlocal(-_roHalfThickness,
			     (2*getCrystalRByRO(roid)-_nCrystalR+1)*_crystalHW,
			     (2*getCrystalZByRO(roid)-_nCrystalZ+1)*_crystalHW );

    return vane.getOrigin() + (vane.getRotation()->inverse())*vlocal;
  }

  //
  // Get crystal axis in Mu2e coordinates - direction from front of 
  // the crystal to readout side
  //
  CLHEP::Hep3Vector Calorimeter::getCrystalAxisByRO(int roid) const {

    int vaneid = getVaneByRO(roid);
    const Vane & vane = getVane(vaneid);

    // Crystal axis in vane coordinates
    CLHEP::Hep3Vector vlocal(1,0,0);

    return (vane.getRotation()->inverse())*vlocal;
  }


} // namespace mu2e
