// $Id: PSShieldMaker.cc,v 1.2 2012/03/30 16:31:10 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/30 16:31:10 $
//
// Original author Andrei Gaponenko

#include <algorithm>

#include "cetlib/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ProductionSolenoidGeom/inc/PSShieldMaker.hh"
#include "ProductionSolenoidGeom/inc/PSShield.hh"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e {

  std::auto_ptr<PSShield> PSShieldMaker::make(const SimpleConfig& c, const CLHEP::Hep3Vector& psEndRefPoint)
  {
    std::vector<double> zPlane, rIn, rOut;
    c.getVectorDouble("PSShield.zPlane", zPlane);
    c.getVectorDouble("PSShield.rIn",    rIn, zPlane.size());
    c.getVectorDouble("PSShield.rOut",   rOut, zPlane.size());

    const double distanceToCryoRefZ = c.getDouble("PSShield.distanceToCryoRefZ")*CLHEP::mm;

    // FIXME: replace the condition with std::is_sorted() once Mu2e adopts C++0x.
    if(std::adjacent_find(zPlane.begin(), zPlane.end(), std::greater<double>()) != zPlane.end()) {
      throw cet::exception("GEOM")<<"PSShieldMaker::make(): coordinates in the zPlane vector must be non-decreasing\n";
    }

    // Put the shield at the required distance to the ref plane
    const CLHEP::Hep3Vector shieldOriginInMu2e(psEndRefPoint + CLHEP::Hep3Vector(0,0,  distanceToCryoRefZ - zPlane[0]));



    std::auto_ptr<PSShield> res(new PSShield(
                                             Polycone(zPlane, rIn, rOut,
                                                      shieldOriginInMu2e,
                                                      c.getString("PSShield.materialName")
                                                      ),

                                             c.getHep3Vector("PSShield.cutout.refPoint"),
                                             c.getDouble("PSShield.cutout.r")*CLHEP::mm,
                                             c.getDouble("PSShield.cutout.halfLengh")*CLHEP::mm,
                                             c.getDouble("PSShield.cutout.rotY")*CLHEP::degree
                                             )
                                );

    if(c.getInt("PSShield.verbosityLevel") > 0) {
      std::cout<<*res.get()<<std::endl;
    }

    return res;
  }

} // namespace mu2e
