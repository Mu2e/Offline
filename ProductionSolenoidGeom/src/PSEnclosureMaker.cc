// $Id: PSEnclosureMaker.cc,v 1.2 2012/03/30 19:18:03 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/30 19:18:03 $
//
// Original author Andrei Gaponenko

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ProductionSolenoidGeom/inc/PSEnclosureMaker.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e {

  std::auto_ptr<PSEnclosure>  PSEnclosureMaker::make(const SimpleConfig& c,
                                                     const CLHEP::Hep3Vector& psEndRefPoint,
                                                     const std::string& psInsideMaterialName)
  {
    const double shellOD = c.getDouble("PSEnclosure.shell.outerDiameter")*CLHEP::mm;
    const double totalLength = c.getDouble("PSEnclosure.length")*CLHEP::mm;
    const double shellThickness = c.getDouble("PSEnclosure.shell.thickness")*CLHEP::mm;
    const double endPlateThickness = c.getDouble("PSEnclosure.endPlate.thickness")*CLHEP::mm;
    const std::string shellMaterialName = c.getString("PSEnclosure.shell.materialName");

    const double shellLength = totalLength - endPlateThickness;

    const CLHEP::Hep3Vector shellOriginInMu2e(psEndRefPoint + CLHEP::Hep3Vector(0,0, -0.5*shellLength));

    std::auto_ptr<PSEnclosure> res(new PSEnclosure(
                                                   // cylindrical shell
                                                   Tube(shellMaterialName,
                                                        shellOriginInMu2e,
                                                        0.5*(shellOD - shellThickness),
                                                        0.5*shellOD,
                                                        0.5*shellLength),

                                                   // vacuum
                                                   Tube(psInsideMaterialName,
                                                        shellOriginInMu2e,
                                                        0.,
                                                        0.5*(shellOD - shellThickness),
                                                        0.5*shellLength
                                                        ),

                                                   // end plate
                                                   Tube(shellMaterialName,
                                                        // The vacuum volume is flush to the PS surface
                                                        psEndRefPoint + CLHEP::Hep3Vector(0,0, -shellLength - 0.5*endPlateThickness),
                                                        0.,
                                                        0.5*shellOD,
                                                        0.5*endPlateThickness
                                                        )
                                                   )

                                   );

    return res;
  }

} // namespace mu2e
