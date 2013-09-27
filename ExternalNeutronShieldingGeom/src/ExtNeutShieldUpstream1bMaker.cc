//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream1bMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream1b.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldUpstream1b>  ExtNeutShieldUpstream1bMaker::make(const 
							      SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("ExtNeutShieldUpstream1b.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("ExtNeutShieldUpstream1b.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("ExtNeutShieldUpstream1b.X1")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1b.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("ExtNeutShieldUpstream1b.X2")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1b.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("ExtNeutShieldUpstream1b.X3")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1b.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("ExtNeutShieldUpstream1b.X4")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1b.Y4")*CLHEP::mm);
    shape.push_back(v4);

    CLHEP::Hep2Vector v5(c.getDouble("ExtNeutShieldUpstream1b.X5")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1b.Y5")*CLHEP::mm);
    shape.push_back(v5);

    CLHEP::Hep2Vector v6(c.getDouble("ExtNeutShieldUpstream1b.X6")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1b.Y6")*CLHEP::mm);
    shape.push_back(v6);

    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("ExtNeutShieldUpstream1b.centerX")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldUpstream1b.centerY")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldUpstream1b.centerZ")
			   *CLHEP::mm);
			   
    double rot( c.getDouble("ExtNeutShieldUpstream1b.rotateX")
		*CLHEP::degree);

    std::unique_ptr<ExtNeutShieldUpstream1b> res(new ExtNeutShieldUpstream1b(
								     shape,
								     leng,
								     mat,
								     site,
								     rot)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldUpstream1b.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
