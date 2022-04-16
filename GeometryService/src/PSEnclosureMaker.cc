//
// Original author Andrei Gaponenko

#include <sstream>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/GeometryService/inc/PSEnclosureMaker.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSEnclosure.hh"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  void PSEnclosureMaker::readWindow(int i, const SimpleConfig& c, const CLHEP::Hep3Vector& winRefPoint, std::unique_ptr<PSEnclosure>& res) {
      std::ostringstream sprefix;
      sprefix<<"PSEnclosure.window"<<1+i<<".";
      const std::string prefix(sprefix.str());

      const double xoff = c.getDouble(prefix+"x"   )*CLHEP::mm;
      const double yoff = c.getDouble(prefix+"y"   )*CLHEP::mm;
      const double zoff = c.getDouble(prefix+"z",0.)*CLHEP::mm;
      const double halfThick = 0.5*c.getDouble(prefix+"thickness")*CLHEP::mm;
      const CLHEP::Hep3Vector windowCenterInMu2e = winRefPoint + CLHEP::Hep3Vector(xoff, yoff, zoff-halfThick);

      const bool windHasFrame = c.getBool(prefix+"hasFrame",false);
      res->hasFrames_.push_back(windHasFrame);
      const bool windHasFrameOut = c.getBool(prefix+"hasFrameOut",false);
      res->hasFramesOut_.push_back(windHasFrameOut);

      if(windHasFrame) {
        const double fWid = c.getDouble(prefix+"frameRadialWidth")*CLHEP::mm;
        const double fHalfThick = 0.5*c.getDouble(prefix+"frameThickness")*CLHEP::mm;
        const double rFin = c.getDouble(prefix+"r")*CLHEP::mm;//Same as rad of window
        const double rFout = rFin+fWid;

        //center the frame on the window, but pushed in z to avoid overlaps with the pipe
        const CLHEP::Hep3Vector frameCenterInMu2e = windowCenterInMu2e + CLHEP::Hep3Vector(0., 0., halfThick-fHalfThick);

        res->wFramesIn_.push_back( Tube(c.getString(prefix+"frameMaterialName"),
                                        frameCenterInMu2e,
                                        rFin, // rIn
                                        rFout,
                                        fHalfThick
                                        ));
        if(windHasFrameOut) {
          const double fOutHalfThick = 0.5*c.getDouble(prefix+"frameOutThickness")*CLHEP::mm;

          const CLHEP::Hep3Vector frameOutCenterInMu2e = frameCenterInMu2e - CLHEP::Hep3Vector(0., 0., fHalfThick+fOutHalfThick);

          res->wFramesOut_.push_back( Tube(c.getString(prefix+"frameOutMaterialName"),
                                           frameOutCenterInMu2e,
                                           rFin, // same radial parameters as the inside section
                                           rFout,
                                           fOutHalfThick
                                           ));
        } else { //ensure the list indices match the window number
          res->wFramesOut_.push_back(Tube());
        }
      } else { //ensure the list indices match the window number
        res->wFramesIn_.push_back( Tube());
        res->wFramesOut_.push_back(Tube());
      }

      //retrieve the window pipe information
      if(res->version() > 2) {
        const double rin   = c.getDouble(prefix+"pipe.rin")*CLHEP::mm;
        const double rout  = c.getDouble(prefix+"pipe.rout")*CLHEP::mm;
        const double thetax = c.getDouble(prefix+"pipe.thetaX")*CLHEP::degree;
        const double thetay = c.getDouble(prefix+"pipe.thetaY")*CLHEP::degree;
        const double zplate = std::abs(res->endPlatePolycone().zPlanes()[0] - res->endPlatePolycone().zPlanes().back());
        const CLHEP::Hep3Vector origin = windowCenterInMu2e; //set at the same
        //origin to ensure it intersects it.
        //Set the pipe to be long enough to intersect both the
        //plate and the window.
        double zlength = std::abs(res->endPlatePolycone().originInMu2e().z() - windowCenterInMu2e.z()) + zplate;
        zlength /= 0.9*std::abs(std::cos(thetax)*std::cos(thetay));

        CLHEP::HepRotation matrix = CLHEP::HepRotation();
        matrix.rotateX(thetax);
        matrix.rotateY(thetay);

        res->windowPipes_.push_back(Tube(rin, rout, zlength, origin, matrix));
      }
      res->windows_.push_back(Tube(c.getString(prefix+"materialName"),
                                   windowCenterInMu2e,
                                   0., // rIn
                                   c.getDouble(prefix+"r")*CLHEP::mm, //rOut
                                   halfThick
                                   ));
  }  // end of readWindow helper function


  std::unique_ptr<PSEnclosure>  PSEnclosureMaker::make(const SimpleConfig& c,
                                                     const CLHEP::Hep3Vector& psEndRefPoint)
  {

    const int vers = c.getInt("PSEnclosure.version",1);

    const double totalLength = c.getDouble("PSEnclosure.length")*CLHEP::mm;
    const double shellThickness = c.getDouble("PSEnclosure.shell.thickness")*CLHEP::mm;
    const double endPlateThickness = c.getDouble("PSEnclosure.endPlate.thickness")*CLHEP::mm;
    const std::string shellMaterialName = c.getString("PSEnclosure.shell.materialName");
    const double flangeIR = c.getDouble("PSEnclosure.flange.rInner",0.0)*CLHEP::mm;
    const double flangeOR = c.getDouble("PSEnclosure.flange.rOuter",0.0)*CLHEP::mm;
    const double flangeThickness = c.getDouble("PSEnclosure.flange.thickness",0.0)*CLHEP::mm;
    const std::string flangeMaterialName = c.getString("PSEnclosure.flange.materialName","StainlessSteel316L");


    std::unique_ptr<PSEnclosure> res = 0;

    if ( vers > 1 ) {
      const double shellODEast = c.getDouble("PSEnclosure.shell.outerDiameterEast")*CLHEP::mm;
      const double shellODWest = c.getDouble("PSEnclosure.shell.outerDiameterWest")*CLHEP::mm;
      const double shellLength = totalLength - endPlateThickness;

      const CLHEP::Hep3Vector shellOriginInMu2e(psEndRefPoint + CLHEP::Hep3Vector(0,0, -0.5*shellLength));

                            // conical frustrum
      Cone shellCone(0.5*(shellODWest - 2*shellThickness),
                     0.5*shellODWest,
                     0.5*(shellODEast - 2*shellThickness),
                     0.5*shellODEast,
                     0.5*shellLength,
                     0., CLHEP::twopi,
                     shellOriginInMu2e,
                     CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY),
                     shellMaterialName);
      if(vers == 2) {
        res = std::make_unique<PSEnclosure>
          (
           shellCone,
           // end plate
           Tube(shellMaterialName,
                // The vacuum volume is flush to the PS surface
                psEndRefPoint + CLHEP::Hep3Vector(0,0, -shellLength - 0.5*endPlateThickness),
                0.,
                0.5*shellODWest,
                0.5*endPlateThickness
                )
           );
      } else {
        std::vector<double> zPlanes; c.getVectorDouble("PSEnclosure.endPlate.zPlanes", zPlanes);
        std::vector<double> rIns   ; c.getVectorDouble("PSEnclosure.endPlate.rIns"   , rIns   );
        std::vector<double> rOuts  ; c.getVectorDouble("PSEnclosure.endPlate.rOuts"  , rOuts  );
        // double zlength = 0.;
        // if(zPlanes.size() > 0 ) {
        //   zlength = abs(zPlanes[0] - zPlanes[zPlanes.size() - 1]); //planes should be ordered
        // }
        res = std::make_unique<PSEnclosure>
          (shellCone,
           // end plate
           Polycone(zPlanes,
                    rIns,
                    rOuts,
                    // The vacuum volume is flush to the PS surface
                    psEndRefPoint + CLHEP::Hep3Vector(0,0, -shellLength),
                    shellMaterialName
                    )
           );
      }


      res->setExtraOffset ( c.getDouble("PSEnclosure.v2.extraZOffset",0.0) );
      if ( vers > 2 ) res->setFlange ( Tube(flangeMaterialName,
                                       psEndRefPoint + CLHEP::Hep3Vector(0,0,-shellLength - 0.5*flangeThickness),
                                       flangeIR,
                                       flangeOR,
                                       0.5*flangeThickness ) );

    } else {
      // This is version 1
      const double shellOD = c.getDouble("PSEnclosure.shell.outerDiameter")*CLHEP::mm;
      const double shellLength = totalLength - endPlateThickness;

      const CLHEP::Hep3Vector shellOriginInMu2e(psEndRefPoint + CLHEP::Hep3Vector(0,0, -0.5*shellLength));

      res = std::make_unique<PSEnclosure>
        (
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
         );
    }

    //----------------------------------------------------------------
    const CLHEP::Hep3Vector winRefPoint = psEndRefPoint + CLHEP::Hep3Vector(0, 0, -totalLength);

    const int nwin = c.getInt("PSEnclosure.nWindows");
    for(int i=0; i<nwin; ++i) {
      readWindow(i, c, winRefPoint, res);
    }

    //----------------------------------------------------------------
    if(c.getInt("PSEnclosure.verbosityLevel") > 0) {
      std::cout<<*res<<std::endl;
    }

    return res;
  }

} // namespace mu2e
