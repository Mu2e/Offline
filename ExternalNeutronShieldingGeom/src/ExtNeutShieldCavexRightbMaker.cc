//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexRightbMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexRightb.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldCavexRightb>  ExtNeutShieldCavexRightbMaker::make(const 
							      SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("ExtNeutShieldCavexRightb.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("ExtNeutShieldCavexRightb.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("ExtNeutShieldCavexRightb.X1")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexRightb.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("ExtNeutShieldCavexRightb.X2")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexRightb.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("ExtNeutShieldCavexRightb.X3")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexRightb.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("ExtNeutShieldCavexRightb.X4")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexRightb.Y4")*CLHEP::mm);
    shape.push_back(v4);


    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("ExtNeutShieldCavexRightb.centerX")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldCavexRightb.centerY")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldCavexRightb.centerZ")
			   *CLHEP::mm);
			   
    
    double rot( c.getDouble("ExtNeutShieldCavexRightb.rotateX")
		*CLHEP::degree);

    std::unique_ptr<ExtNeutShieldCavexRightb> res(new ExtNeutShieldCavexRightb(
								     shape,
								     leng,
								     mat,
								     site,
								     rot)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldCavexRightb.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
