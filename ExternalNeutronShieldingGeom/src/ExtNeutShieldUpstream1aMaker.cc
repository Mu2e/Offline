//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream1aMaker.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream1a.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<ExtNeutShieldUpstream1a>  ExtNeutShieldUpstream1aMaker::make(const 
							      SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("ExtNeutShieldUpstream1a.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("ExtNeutShieldUpstream1a.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("ExtNeutShieldUpstream1a.X1")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1a.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("ExtNeutShieldUpstream1a.X2")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1a.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("ExtNeutShieldUpstream1a.X3")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1a.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("ExtNeutShieldUpstream1a.X4")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1a.Y4")*CLHEP::mm);
    shape.push_back(v4);

    CLHEP::Hep2Vector v5(c.getDouble("ExtNeutShieldUpstream1a.X5")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1a.Y5")*CLHEP::mm);
    shape.push_back(v5);

    CLHEP::Hep2Vector v6(c.getDouble("ExtNeutShieldUpstream1a.X6")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1a.Y6")*CLHEP::mm);
    shape.push_back(v6);

    CLHEP::Hep2Vector v7(c.getDouble("ExtNeutShieldUpstream1a.X7")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1a.Y7")*CLHEP::mm);
    shape.push_back(v7);

    CLHEP::Hep2Vector v8(c.getDouble("ExtNeutShieldUpstream1a.X8")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1a.Y8")*CLHEP::mm);
    shape.push_back(v8);

    CLHEP::Hep2Vector v9(c.getDouble("ExtNeutShieldUpstream1a.X9")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1a.Y9")*CLHEP::mm);
    shape.push_back(v9);

    CLHEP::Hep2Vector v10(c.getDouble("ExtNeutShieldUpstream1a.X10")*CLHEP::mm,
		  c.getDouble("ExtNeutShieldUpstream1a.Y10")*CLHEP::mm);
    shape.push_back(v10);

    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("ExtNeutShieldUpstream1a.centerX")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldUpstream1a.centerY")
			   *CLHEP::mm,
			   c.getDouble("ExtNeutShieldUpstream1a.centerZ")
			   *CLHEP::mm);
			   
    
    double rot( c.getDouble("ExtNeutShieldUpstream1a.rotateX")
		*CLHEP::degree);

    std::unique_ptr<ExtNeutShieldUpstream1a> res(new ExtNeutShieldUpstream1a(
								     shape,
								     leng,
								     mat,
								     site,
								     rot)
					     );

    //----------------------------------------------------------------
    if(c.getInt("ExtNeutShieldUpstream1a.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
