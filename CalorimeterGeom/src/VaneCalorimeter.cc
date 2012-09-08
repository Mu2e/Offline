//
// Geometry and identifier info about the VaneCalorimeter.
//
//
// $Id: VaneCalorimeter.cc,v 1.1 2012/09/08 02:24:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/09/08 02:24:25 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include "CalorimeterGeom/inc/VaneCalorimeter.hh"

using namespace std;


namespace mu2e {

  //
  // Convert coordinates from Mu2e frame to local frame of crystal, identified by roid
  //
  CLHEP::Hep3Vector VaneCalorimeter::toCrystalFrame(int roid, CLHEP::Hep3Vector const& pos) const 
  {
  
    int vaneid = getVaneByRO(roid);
    const Vane & vane = getVane(vaneid);

    double crystalUnitWidth = _crystalHW + _wrapperThickness + _shellThickness;

    CLHEP::Hep3Vector vlocal(_roHalfThickness,
                             (2*getCrystalRByRO(roid)-_nCrystalR+1)*crystalUnitWidth,
                             (2*getCrystalZByRO(roid)-_nCrystalZ+1)*crystalUnitWidth );

    return *(vane.getRotation())*(pos-vane.getOrigin())-vlocal;
  }

  //---------------------------------------------------------
  //
  // Convert coordinates from Mu2e frame to local frame of the vane, identified by vaneid
  //

  CLHEP::Hep3Vector VaneCalorimeter::toVaneFrame(int vaneid, CLHEP::Hep3Vector const& pos) const 
  {
          const Vane & vane = getVane(vaneid);
          double crystalUnitWidth  = _crystalHW + _wrapperThickness + _shellThickness;
          double crystalUnitLength = _crystalHL + _wrapperThickness;

          CLHEP::Hep3Vector vlocal(_roHalfThickness - crystalUnitLength,
                          (_nCrystalR )*crystalUnitWidth,
                          (_nCrystalZ )*crystalUnitWidth );

          return *(vane.getRotation())*(pos-vane.getOrigin())+vlocal;
  }

  CLHEP::Hep3Vector VaneCalorimeter::fromVaneFrame(int vaneid, CLHEP::Hep3Vector const& pos) const 
  {
            const Vane & vane = getVane(vaneid);

            double crystalUnitWidth  = _crystalHW + _wrapperThickness + _shellThickness;
            double crystalUnitLength = _crystalHL + _wrapperThickness;

            CLHEP::Hep3Vector vlocal(_roHalfThickness - crystalUnitLength,
                            (_nCrystalR )*crystalUnitWidth,
                            (_nCrystalZ )*crystalUnitWidth );

            return (vane.getRotation()->inverse())*(pos-vlocal) + vane.getOrigin();
    }



  //---------------------------------------------------------
  //
  // Get crystal origin (center) in Mu2e coordinates
  //
  CLHEP::Hep3Vector VaneCalorimeter::getCrystalOriginByRO(int roid) const 
  {

    int vaneid = getVaneByRO(roid);
    const Vane & vane = getVane(vaneid);

    double crystalUnitWidth  = _crystalHW + _wrapperThickness + _shellThickness;

    // Crystal center in vane coordinates
    CLHEP::Hep3Vector vlocal(_roHalfThickness,
                             (2*getCrystalRByRO(roid)-_nCrystalR+1)*crystalUnitWidth,
                             (2*getCrystalZByRO(roid)-_nCrystalZ+1)*crystalUnitWidth );

    return vane.getOrigin() + (vane.getRotation()->inverse())*vlocal;
  }

  //
  // Get crystal axis in Mu2e coordinates - direction from front of
  // the crystal to readout side
  //
  CLHEP::Hep3Vector VaneCalorimeter::getCrystalAxisByRO(int roid) const 
  {

    int vaneid = getVaneByRO(roid);
    const Vane & vane = getVane(vaneid);

    // Crystal axis in vane coordinates
    CLHEP::Hep3Vector vlocal(1,0,0);

    return (vane.getRotation()->inverse())*vlocal;
  }


} // namespace mu2e
