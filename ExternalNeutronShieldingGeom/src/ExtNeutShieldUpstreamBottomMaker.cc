//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstreamBottomMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstreamBottom.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldUpstreamBottom>  ExtNeutShieldUpstreamBottomMaker::make(const 
							      SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("ExtNeutShieldUpstreamBottom.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("ExtNeutShieldUpstreamBottom.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("ExtNeutShieldUpstreamBottom.X1")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstreamBottom.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("ExtNeutShieldUpstreamBottom.X2")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstreamBottom.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("ExtNeutShieldUpstreamBottom.X3")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstreamBottom.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("ExtNeutShieldUpstreamBottom.X4")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstreamBottom.Y4")*CLHEP::mm);
    shape.push_back(v4);


    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("ExtNeutShieldUpstreamBottom.centerX")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldUpstreamBottom.centerY")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldUpstreamBottom.centerZ")
			   *CLHEP::mm);
			   
    
    double rot( c.getDouble("ExtNeutShieldUpstreamBottom.rotateX")
		*CLHEP::degree);

    std::unique_ptr<ExtNeutShieldUpstreamBottom> res(new ExtNeutShieldUpstreamBottom(
								     shape,
								     leng,
								     mat,
								     site,
								     rot)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldUpstreamBottom.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
