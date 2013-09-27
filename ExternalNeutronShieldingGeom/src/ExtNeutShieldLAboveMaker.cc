//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldLAboveMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldLAbove.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldLAbove>  ExtNeutShieldLAboveMaker::make(const 
							      SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("ExtNeutShieldLAbove.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("ExtNeutShieldLAbove.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("ExtNeutShieldLAbove.X1")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldLAbove.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("ExtNeutShieldLAbove.X2")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldLAbove.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("ExtNeutShieldLAbove.X3")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldLAbove.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("ExtNeutShieldLAbove.X4")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldLAbove.Y4")*CLHEP::mm);
    shape.push_back(v4);

    CLHEP::Hep2Vector v5(c.getDouble("ExtNeutShieldLAbove.X5")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldLAbove.Y5")*CLHEP::mm);
    shape.push_back(v5);

    CLHEP::Hep2Vector v6(c.getDouble("ExtNeutShieldLAbove.X6")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldLAbove.Y6")*CLHEP::mm);
    shape.push_back(v6);


    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("ExtNeutShieldLAbove.centerX")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldLAbove.centerY")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldLAbove.centerZ")
			   *CLHEP::mm);
			   
    
    double rot( c.getDouble("ExtNeutShieldLAbove.rotateX")
		*CLHEP::degree);

    std::unique_ptr<ExtNeutShieldLAbove> res(new ExtNeutShieldLAbove(
								     shape,
								     leng,
								     mat,
								     site,
								     rot)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldLAbove.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
