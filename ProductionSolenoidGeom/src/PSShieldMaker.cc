// $Id: PSShieldMaker.cc,v 1.5 2012/07/15 22:06:18 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:18 $
//
// Original author Andrei Gaponenko

#include <algorithm>
#include <sstream>

#include "cetlib/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ProductionSolenoidGeom/inc/PSShieldMaker.hh"
#include "ProductionSolenoidGeom/inc/PSShield.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  PSShield::Groove PSShieldMaker::readGroove(int i, const SimpleConfig& c) {
      std::ostringstream prefix;
      prefix<<"PSShield.groove"<<1+i<<".";

      return PSShield::Groove(
                              c.getHep3Vector(prefix.str()+"refPoint"),
                              c.getDouble(prefix.str()+"theta")*CLHEP::degree,
                              c.getDouble(prefix.str()+"phi")*CLHEP::degree,
                              c.getDouble(prefix.str()+"r")*CLHEP::mm,
                              c.getDouble(prefix.str()+"halfLengh")*CLHEP::mm
                              );
  }

  std::auto_ptr<PSShield> PSShieldMaker::make(const SimpleConfig& c,
                                              const CLHEP::Hep3Vector& psEndRefPoint,
                                              const CLHEP::Hep3Vector& productionTargetCenter
                                              )
  {
    std::vector<double> zPlane, rIn, rOut;
    c.getVectorDouble("PSShield.zPlane", zPlane);
    c.getVectorDouble("PSShield.rIn",    rIn, zPlane.size());
    c.getVectorDouble("PSShield.rOut",   rOut, zPlane.size());

    // FIXME: replace the condition with std::is_sorted() once Mu2e adopts C++0x.
    if(std::adjacent_find(zPlane.begin(), zPlane.end(), std::greater<double>()) != zPlane.end()) {
      throw cet::exception("GEOM")<<"PSShieldMaker::make(): coordinates in the zPlane vector must be non-decreasing\n";
    }

    // Compute placement of the shield
    const CLHEP::Hep3Vector shieldOriginInMu2e(psEndRefPoint.x(),

                                               psEndRefPoint.y(),

                                               productionTargetCenter.z()
                                               + c.getDouble("PSShield.zOffsetFromProductionTarget")
                                               - zPlane[0]
                                               );

    std::auto_ptr<PSShield> res(new PSShield(Polycone(
                                                      zPlane, rIn, rOut,
                                                      shieldOriginInMu2e,
                                                      c.getString("PSShield.materialName")
                                                      )
                                             )
                                );

    const int nGrooves = c.getInt("PSShield.nGrooves");
    for(int i=0; i<nGrooves; ++i) {
      res->grooves_.push_back(readGroove(i, c));
    }

    if(c.getInt("PSShield.verbosityLevel") > 0) {
      std::cout<<*res.get()<<std::endl;
    }

    return res;
  }

} // namespace mu2e
