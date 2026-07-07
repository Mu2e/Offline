//
// Original author Andrei Gaponenko
// Extensive modifications by David N. Brown, Louisville, 2017
#include <algorithm>
#include <sstream>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/GeometryService/inc/PSShieldMaker.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSShield.hh"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"


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

  std::unique_ptr<PSShield> PSShieldMaker::make(const SimpleConfig& c,
                                              const CLHEP::Hep3Vector& psEndRefPoint,
                                              const CLHEP::Hep3Vector& productionTargetCenter
                                              )
  {
    std::unique_ptr<PSShield> res(new PSShield());

    int myVersion = c.getInt("PSShield.version",1);

    res->version_ = myVersion;

    std::vector<double> zPlane;
    c.getVectorDouble("PSShield.zPlane", zPlane);
    if(!std::is_sorted(zPlane.begin(), zPlane.end())) {
      throw cet::exception("GEOM")<<"PSShieldMaker::make(): coordinates in the zPlane vector must be non-decreasing\n";
    }

    // Compute placement of the shield
    const CLHEP::Hep3Vector shieldOriginInMu2e(psEndRefPoint.x(),

                                               psEndRefPoint.y(),

                                               productionTargetCenter.z()
                                               + c.getDouble("PSShield.zOffsetFromProductionTarget")
                                               - zPlane[0]
                                               );
    c.getInt("PSShield.verbosityLevel") > 0 && std::cout << " PSShieldMaker " << __func__ << productionTargetCenter.z()
                                                      << " " << c.getDouble("PSShield.zOffsetFromProductionTarget") << " " << zPlane[0] << std::endl;
    //----------------------------------------------------------------
    // Read in the shells

    const int nShells = c.getInt("PSShield.nShells");
    res->shells_.reserve(nShells);
    for(int ishell = 1; ishell <= nShells; ++ishell) {
      std::ostringstream osinner, osouter;
      osinner << ishell;
      osouter << ishell + 1;

      std::vector<double> rIn, rOut;
      c.getVectorDouble("PSShield.r"+osinner.str(), rIn, zPlane.size());
      c.getVectorDouble("PSShield.r"+osouter.str(), rOut, zPlane.size());
      const std::string material = c.getString("PSShield.material"+osinner.str());

      res->shells_.emplace_back(zPlane, rIn, rOut, shieldOriginInMu2e, material);
    }

    //----------------------------------------------------------------
    // Read in the grooves

    const int nGrooves = c.getInt("PSShield.nGrooves");
    res->grooves_.reserve(nGrooves);
    for(int i=0; i<nGrooves; ++i) {
      res->grooves_.emplace_back(readGroove(i, c));
    }

    if(c.getInt("PSShield.verbosityLevel") > 0) {
      std::cout<<*res.get()<<std::endl;
    }

    //----------------------------------------------------------------
    // Read in the upstream end ring information.
    // DNB, Louisville, Jan 2017
    // There are two ends, so there are two end rings.  Each is a polycone
    if ( myVersion > 1 ) {
      res->endRings_.reserve(4);

      // First do upstream ("front") end ring
      std::vector<double> zPf;
      c.getVectorDouble("PSShield.frontRing.zPlane", zPf);
      if(!std::is_sorted(zPf.begin(), zPf.end())) {
        throw cet::exception("GEOM")<<"PSShieldMaker::make(): coordinates in the frontRing zPlane vector must be non-decreasing\n";
      }
      std::vector<double> rIf, rOf;
      c.getVectorDouble("PSShield.frontRing.rIn", rIf, zPf.size());
      c.getVectorDouble("PSShield.frontRing.rOut", rOf, zPf.size());
      const std::string materialf = c.getString("PSShield.frontRing.material");

      res->endRings_.emplace_back(zPf, rIf, rOf, shieldOriginInMu2e, materialf );

      // Now do downstream ("back") end ring.
      std::vector<double> zPb;
      c.getVectorDouble("PSShield.backRing.zPlane", zPb);
      if(!std::is_sorted(zPb.begin(), zPb.end())) {
        throw cet::exception("GEOM")<<"PSShieldMaker::make(): coordinates in the backRing zPlane vector must be non-decreasing\n";
      }
      std::vector<double> rIb, rOb;
      c.getVectorDouble("PSShield.backRing.rIn", rIb, zPb.size());
      c.getVectorDouble("PSShield.backRing.rOut", rOb, zPb.size());
      const std::string materialb = c.getString("PSShield.backRing.material");

      res->endRings_.emplace_back(zPb, rIb, rOb, shieldOriginInMu2e, materialb );

      // Now the water ring
      std::vector<double> zPw;
      c.getVectorDouble("PSShield.waterRing.zPlane", zPw);
      if(!std::is_sorted(zPw.begin(), zPw.end())) {
        throw cet::exception("GEOM")<<"PSShieldMaker::make(): coordinates in the waterRing zPlane vector must be non-decreasing\n";
      }
      std::vector<double> rIw, rOw;
      c.getVectorDouble("PSShield.waterRing.rIn", rIw, zPw.size());
      c.getVectorDouble("PSShield.waterRing.rOut", rOw, zPw.size());
      const std::string materialw = c.getString("PSShield.waterRing.material");

      res->endRings_.emplace_back(zPw, rIw, rOw, shieldOriginInMu2e, materialw );

      // And the extension of the sheath around the HRS proper
      std::vector<double> zPs;
      c.getVectorDouble("PSShield.sheathRing.zPlane", zPs);
      if(!std::is_sorted(zPs.begin(), zPs.end())) {
        throw cet::exception("GEOM")<<"PSShieldMaker::make(): coordinates in the sheathRing zPlane vector must be non-decreasing\n";
      }
      std::vector<double> rIs, rOs;
      c.getVectorDouble("PSShield.sheathRing.rIn", rIs, zPs.size());
      c.getVectorDouble("PSShield.sheathRing.rOut", rOs, zPs.size());
      const std::string materials = c.getString("PSShield.sheathRing.material");

      res->endRings_.emplace_back(zPs, rIs, rOs, shieldOriginInMu2e, materials );

    } // If version > 1

    //----------------------------------------------------------------
    // Read in the proton beam inlet information.  This is a tube
    // representing the beam pipe extended into HRS.
    // Added by David Norvil Brown, Louisville, March 2015

    std::string pipeMaterial = c.getString("PSShield.inlet.materialName");
    double pipeIR = c.getDouble("PSShield.inlet.innerR")*CLHEP::mm;
    double pipeOR = c.getDouble("PSShield.inlet.outerR")*CLHEP::mm;
    double pipeLength = c.getDouble("PSShield.inlet.length")*CLHEP::mm;
    double angle1 = c.getDouble("PSShield.inlet.angleY")*CLHEP::degree;
    double angle2 = c.getDouble("PSShield.inlet.angleX")*CLHEP::degree;
    bool createBeamPipe = c.getBool("PSShield.inlet.createBeamPipe");

    res->beamAngleY_ = angle1;
    res->beamAngleX_ = angle2;
    res->beamInletCenter_ = c.getHep3Vector("PSShield.inlet.center");
    CLHEP::HepRotation pipeRotat(CLHEP::HepRotation::IDENTITY);
    pipeRotat.rotateY(angle1);
    pipeRotat.rotateX(angle2);

    res->beamInlet_ = Tube(pipeMaterial, shieldOriginInMu2e,
                           pipeIR, pipeOR, pipeLength/2.0, 0, CLHEP::twopi,
                           pipeRotat);

    res->createBeamPipe_ = createBeamPipe;

    return res;
  }

} // namespace mu2e
