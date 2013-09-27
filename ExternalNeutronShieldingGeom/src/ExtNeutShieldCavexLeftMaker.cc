//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexLeftMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexLeft.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldCavexLeft>  ExtNeutShieldCavexLeftMaker::make(const 
							      SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("ExtNeutShieldCavexLeft.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("ExtNeutShieldCavexLeft.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("ExtNeutShieldCavexLeft.X1")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexLeft.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("ExtNeutShieldCavexLeft.X2")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexLeft.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("ExtNeutShieldCavexLeft.X3")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexLeft.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("ExtNeutShieldCavexLeft.X4")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexLeft.Y4")*CLHEP::mm);
    shape.push_back(v4);

    CLHEP::Hep2Vector v5(c.getDouble("ExtNeutShieldCavexLeft.X5")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexLeft.Y5")*CLHEP::mm);
    shape.push_back(v5);

    CLHEP::Hep2Vector v6(c.getDouble("ExtNeutShieldCavexLeft.X6")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexLeft.Y6")*CLHEP::mm);
    shape.push_back(v6);


    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("ExtNeutShieldCavexLeft.centerX")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldCavexLeft.centerY")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldCavexLeft.centerZ")
			   *CLHEP::mm);
			   
    
    double rot( c.getDouble("ExtNeutShieldCavexLeft.rotateX")
		*CLHEP::degree);

    std::unique_ptr<ExtNeutShieldCavexLeft> res(new ExtNeutShieldCavexLeft(
								     shape,
								     leng,
								     mat,
								     site,
								     rot)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldCavexLeft.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
