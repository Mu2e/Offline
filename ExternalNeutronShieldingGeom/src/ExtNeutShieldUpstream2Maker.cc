//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream2Maker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream2.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldUpstream2>  ExtNeutShieldUpstream2Maker::make(const 
							      SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("ExtNeutShieldUpstream2.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("ExtNeutShieldUpstream2.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("ExtNeutShieldUpstream2.X1")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream2.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("ExtNeutShieldUpstream2.X2")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream2.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("ExtNeutShieldUpstream2.X3")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream2.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("ExtNeutShieldUpstream2.X4")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream2.Y4")*CLHEP::mm);
    shape.push_back(v4);


    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("ExtNeutShieldUpstream2.centerX")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldUpstream2.centerY")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldUpstream2.centerZ")
			   *CLHEP::mm);
			   
    
    double rot( c.getDouble("ExtNeutShieldUpstream2.rotateX")
		*CLHEP::degree);

    std::unique_ptr<ExtNeutShieldUpstream2> res(new ExtNeutShieldUpstream2(
								     shape,
								     leng,
								     mat,
								     site,
								     rot)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldUpstream2.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
