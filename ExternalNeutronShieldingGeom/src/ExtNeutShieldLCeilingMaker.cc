//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldLCeilingMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldLCeiling.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldLCeiling>  ExtNeutShieldLCeilingMaker::make(const 
							      SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("ExtNeutShieldLCeiling.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("ExtNeutShieldLCeiling.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("ExtNeutShieldLCeiling.X1")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldLCeiling.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("ExtNeutShieldLCeiling.X2")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldLCeiling.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("ExtNeutShieldLCeiling.X3")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldLCeiling.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("ExtNeutShieldLCeiling.X4")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldLCeiling.Y4")*CLHEP::mm);
    shape.push_back(v4);


    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("ExtNeutShieldLCeiling.centerX")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldLCeiling.centerY")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldLCeiling.centerZ")
			   *CLHEP::mm);
			   
    
    double rot( c.getDouble("ExtNeutShieldLCeiling.rotateX")
		*CLHEP::degree);

    std::unique_ptr<ExtNeutShieldLCeiling> res(new ExtNeutShieldLCeiling(
								     shape,
								     leng,
								     mat,
								     site,
								     rot)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldLCeiling.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
