//
// Original author Andrei Gaponenko

#include <sstream>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "GeometryService/inc/PSEnclosureMaker.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  namespace {
    Tube readWindow(int i, const SimpleConfig& c, const CLHEP::Hep3Vector& winRefPoint) {
      std::ostringstream sprefix;
      sprefix<<"PSEnclosure.window"<<1+i<<".";
      const std::string prefix(sprefix.str());

      const double xoff = c.getDouble(prefix+"x")*CLHEP::mm;
      const double yoff = c.getDouble(prefix+"y")*CLHEP::mm;
      const double halfThick = 0.5*c.getDouble(prefix+"thickness")*CLHEP::mm;

      const CLHEP::Hep3Vector windowCenterInMu2e = winRefPoint + CLHEP::Hep3Vector(xoff, yoff, -halfThick);

      return Tube(c.getString(prefix+"materialName"),
                  windowCenterInMu2e,
                  0., // rIn
                  c.getDouble(prefix+"r")*CLHEP::mm,
                  halfThick
                  );
    }
  }


  std::unique_ptr<PSEnclosure>  PSEnclosureMaker::make(const SimpleConfig& c,
                                                     const CLHEP::Hep3Vector& psEndRefPoint)
  {

    const int vers = c.getInt("PSEnclosure.version",1);

    const double totalLength = c.getDouble("PSEnclosure.length")*CLHEP::mm;
    const double shellThickness = c.getDouble("PSEnclosure.shell.thickness")*CLHEP::mm;
    const double endPlateThickness = c.getDouble("PSEnclosure.endPlate.thickness")*CLHEP::mm;
    const std::string shellMaterialName = c.getString("PSEnclosure.shell.materialName");

    std::unique_ptr<PSEnclosure> res = 0;

    if ( vers > 1 ) {
      const double shellODEast = c.getDouble("PSEnclosure.shell.outerDiameterEast")*CLHEP::mm;
      const double shellODWest = c.getDouble("PSEnclosure.shell.outerDiameterWest")*CLHEP::mm;
      const double shellLength = totalLength - endPlateThickness;

      const CLHEP::Hep3Vector shellOriginInMu2e(psEndRefPoint + CLHEP::Hep3Vector(0,0, -0.5*shellLength));

      res = std::unique_ptr<PSEnclosure>(new PSEnclosure(
			    // conical frustrum
			    Cone(0.5*(shellODWest - 2*shellThickness),
				 0.5*shellODWest,
				 0.5*(shellODEast - 2*shellThickness),
				 0.5*shellODEast,
				 0.5*shellLength,
				 0., CLHEP::twopi,
				 shellOriginInMu2e,
				 CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY),
				 shellMaterialName),

			    // end plate
			    Tube(shellMaterialName,
				 // The vacuum volume is flush to the PS surface
				 psEndRefPoint + CLHEP::Hep3Vector(0,0, -shellLength - 0.5*endPlateThickness),
				 0.,
				 0.5*shellODWest,
				 0.5*endPlateThickness
				 )                              
							 )
			    );

      res->setExtraOffset ( c.getDouble("PSEnclosure.v2.extraZOffset",0.0) );

    } else {
      const double shellOD = c.getDouble("PSEnclosure.shell.outerDiameter")*CLHEP::mm;
      const double shellLength = totalLength - endPlateThickness;

      const CLHEP::Hep3Vector shellOriginInMu2e(psEndRefPoint + CLHEP::Hep3Vector(0,0, -0.5*shellLength));

      res = std::unique_ptr<PSEnclosure> (new PSEnclosure(
			    // cylindrical shell
			    Tube(shellMaterialName,
				 shellOriginInMu2e,
				 0.5*(shellOD - 2*shellThickness),
				 0.5*shellOD,
				 0.5*shellLength),

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
    }
    //----------------------------------------------------------------
    const CLHEP::Hep3Vector winRefPoint = psEndRefPoint + CLHEP::Hep3Vector(0, 0, -totalLength);

    const int nwin = c.getInt("PSEnclosure.nWindows");
    for(int i=0; i<nwin; ++i) {
      res->windows_.push_back(readWindow(i, c, winRefPoint));
    }

    //----------------------------------------------------------------
    if(c.getInt("PSEnclosure.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
