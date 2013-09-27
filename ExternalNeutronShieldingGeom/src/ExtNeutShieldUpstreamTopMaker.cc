//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstreamTopMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstreamTop.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldUpstreamTop>  ExtNeutShieldUpstreamTopMaker::make(const 
							      SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("ExtNeutShieldUpstreamTop.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("ExtNeutShieldUpstreamTop.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("ExtNeutShieldUpstreamTop.X1")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstreamTop.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("ExtNeutShieldUpstreamTop.X2")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstreamTop.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("ExtNeutShieldUpstreamTop.X3")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstreamTop.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("ExtNeutShieldUpstreamTop.X4")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstreamTop.Y4")*CLHEP::mm);
    shape.push_back(v4);


    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("ExtNeutShieldUpstreamTop.centerX")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldUpstreamTop.centerY")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldUpstreamTop.centerZ")
			   *CLHEP::mm);
			   
    
    double rot( c.getDouble("ExtNeutShieldUpstreamTop.rotateX")
		*CLHEP::degree);

    std::unique_ptr<ExtNeutShieldUpstreamTop> res(new ExtNeutShieldUpstreamTop(
								     shape,
								     leng,
								     mat,
								     site,
								     rot)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldUpstreamTop.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
