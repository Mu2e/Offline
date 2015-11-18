//
// Original author David Norvil Brown

#include <string>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "GeometryService/inc/PSExternalShieldingMaker.hh"
#include "ProductionSolenoidGeom/inc/PSExternalShielding.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //  namespace {

  std::unique_ptr<PSExternalShielding>  PSExternalShieldingMaker::make(const
                                                              SimpleConfig& c)
  {
    // Get the material name
    std::string mat = c.getString("PSExternalShielding.materialName");

    // Get the FULL length in z of the external shield
    double leng = c.getDouble("PSExternalShielding.length")*CLHEP::mm;

    // Now retrieve the vertices and put into a vector
    std::vector<CLHEP::Hep2Vector> shape;
    CLHEP::Hep2Vector v1(c.getDouble("PSExternalShielding.X1")*CLHEP::mm,
                  c.getDouble("PSExternalShielding.Y1")*CLHEP::mm);
    shape.push_back(v1);

    CLHEP::Hep2Vector v2(c.getDouble("PSExternalShielding.X2")*CLHEP::mm,
                  c.getDouble("PSExternalShielding.Y2")*CLHEP::mm);
    shape.push_back(v2);

    CLHEP::Hep2Vector v3(c.getDouble("PSExternalShielding.X3")*CLHEP::mm,
                  c.getDouble("PSExternalShielding.Y3")*CLHEP::mm);
    shape.push_back(v3);

    CLHEP::Hep2Vector v4(c.getDouble("PSExternalShielding.X4")*CLHEP::mm,
                  c.getDouble("PSExternalShielding.Y4")*CLHEP::mm);
    shape.push_back(v4);

    CLHEP::Hep2Vector v5(c.getDouble("PSExternalShielding.X5")*CLHEP::mm,
                  c.getDouble("PSExternalShielding.Y5")*CLHEP::mm);
    shape.push_back(v5);

    CLHEP::Hep2Vector v6(c.getDouble("PSExternalShielding.X6")*CLHEP::mm,
                  c.getDouble("PSExternalShielding.Y6")*CLHEP::mm);
    shape.push_back(v6);

    CLHEP::Hep2Vector v7(c.getDouble("PSExternalShielding.X7")*CLHEP::mm,
                  c.getDouble("PSExternalShielding.Y7")*CLHEP::mm);
    shape.push_back(v7);

    CLHEP::Hep2Vector v8(c.getDouble("PSExternalShielding.X8")*CLHEP::mm,
                  c.getDouble("PSExternalShielding.Y8")*CLHEP::mm);
    shape.push_back(v8);

    // Location of the center of the shape relative to its parent volume
    // (the detector hall).
    CLHEP::Hep3Vector site(c.getDouble("PSExternalShielding.centerX")
                           *CLHEP::mm,
                           c.getDouble("PSExternalShielding.centerY")
                           *CLHEP::mm,
                           c.getDouble("PSExternalShielding.centerZ")
                           *CLHEP::mm);


    std::unique_ptr<PSExternalShielding> res(new PSExternalShielding(
                                                                     shape,
                                                                     leng,
                                                                     mat,
                                                                     site)
                                             );

    //----------------------------------------------------------------
    if(c.getInt("PSExternalShielding.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
