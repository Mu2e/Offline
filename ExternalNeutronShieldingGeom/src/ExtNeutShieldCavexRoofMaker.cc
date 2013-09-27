//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexRoofMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexRoof.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldCavexRoof>  ExtNeutShieldCavexRoofMaker::make(const 
							      SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("ExtNeutShieldCavexRoof.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("ExtNeutShieldCavexRoof.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("ExtNeutShieldCavexRoof.X1")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexRoof.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("ExtNeutShieldCavexRoof.X2")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexRoof.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("ExtNeutShieldCavexRoof.X3")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexRoof.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("ExtNeutShieldCavexRoof.X4")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexRoof.Y4")*CLHEP::mm);
    shape.push_back(v4);

    CLHEP::Hep2Vector v5(c.getDouble("ExtNeutShieldCavexRoof.X5")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexRoof.Y5")*CLHEP::mm);
    shape.push_back(v5);

    CLHEP::Hep2Vector v6(c.getDouble("ExtNeutShieldCavexRoof.X6")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldCavexRoof.Y6")*CLHEP::mm);
    shape.push_back(v6);


    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("ExtNeutShieldCavexRoof.centerX")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldCavexRoof.centerY")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldCavexRoof.centerZ")
			   *CLHEP::mm);
			   
    
    double rot( c.getDouble("ExtNeutShieldCavexRoof.rotateX")
		*CLHEP::degree);

    std::unique_ptr<ExtNeutShieldCavexRoof> res(new ExtNeutShieldCavexRoof(
								     shape,
								     leng,
								     mat,
								     site,
								     rot)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldCavexRoof.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
