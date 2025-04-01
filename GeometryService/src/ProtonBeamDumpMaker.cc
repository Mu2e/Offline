// Andrei Gaponenko, 2011

#include "Offline/GeometryService/inc/ProtonBeamDumpMaker.hh"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

#include "cetlib_except/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

namespace mu2e {

  //================================================================
  std::unique_ptr<ProtonBeamDump> ProtonBeamDumpMaker::make(const SimpleConfig& c,
                                                            const Mu2eHall& hall)
  {

    std::unique_ptr<ProtonBeamDump> dump(new ProtonBeamDump());

    int verbose = c.getInt("protonBeamDump.verbosityLevel", 0);

    // Get relevant Hall solids
    ExtrudedSolid psArea  = hall.getBldgSolid("psArea");        // contains front corners of dump slab
    ExtrudedSolid psWall  = hall.getBldgSolid("psWallUpper");   // contains back corners of dump slab
    ExtrudedSolid psWallS = hall.getBldgSolid("psAreaUpperS");  // contains boundary between front and back
    ExtrudedSolid psWallN = hall.getBldgSolid("psAreaUpperN");  // contains boundary between front and back
    ExtrudedSolid psCeil;
    psCeil = hall.getBldgSolid("psAreaCeilingSW");

    // position
    dump->_coreCenterInMu2e = c.getHep3Vector("protonBeamDump.coreCenterInMu2e");
    const double coreRotY = dump->_coreRotY = c.getDouble("protonBeamDump.coreRotY") * CLHEP::degree;
    dump->_coreRotationInMu2e.rotateY(coreRotY);

    std::vector<double> additionalSteel;
    c.getVectorDouble("protonBeamDump.coreHalfSize", dump->_coreHalfSize, 3);
    c.getVectorDouble("protonBeamDump.neutronCaveHalfSize", dump->_neutronCaveHalfSize, 3);
    c.getVectorDouble("protonBeamDump.mouthHalfSize", dump->_mouthHalfSize, 3);
    c.getVectorDouble("protonBeamDump.additionalSteel", additionalSteel, 3);
    dump->_minCoreShieldingThickness = c.getDouble("protonBeamDump.minCoreShieldingThickness");
    const double coreAirGap = c.getDouble("protonBeamDump.coreAirGap");

    const CLHEP::Hep3Vector& offset = psWall.getOffsetFromMu2eOrigin();
    dump->_shieldingFaceYmin = offset[1] - psWall.getYhalfThickness();
    const CLHEP::Hep3Vector& ceiling = psCeil.getOffsetFromMu2eOrigin();

    // Compute the overall size of the nominal dump front portion
    dump->_frontShieldingHalfSize.resize(3);
    dump->_frontShieldingHalfSize[0] = dump->_coreHalfSize[0] + dump->_minCoreShieldingThickness;
    dump->_frontShieldingHalfSize[1] = (ceiling[1] - psCeil.getYhalfThickness() - dump->_shieldingFaceYmin)/2.0;
    dump->_frontShieldingHalfSize[2] = dump->_coreHalfSize[2] + dump->_neutronCaveHalfSize[2] + dump->_mouthHalfSize[2];
    if(verbose) {
      std::cout<<"ProtonBeamDumpMaker"<<": ProtonBeamDump frontShielding half size = ";
      std::copy(dump->_frontShieldingHalfSize.begin(), dump->_frontShieldingHalfSize.end(), std::ostream_iterator<double>(std::cout, ", "));
      std::cout<<std::endl;
    }

    dump->_frontSteelHalfSize.resize(3);
    dump->_frontSteelHalfSize[0] = additionalSteel[0]/2.0;
    dump->_frontSteelHalfSize[1] = additionalSteel[1]/2.0;
    dump->_frontSteelHalfSize[2] = dump->_coreHalfSize[2];
    dump->_backSteelHalfSize.resize(3);
    dump->_backSteelHalfSize[0] = dump->_frontSteelHalfSize[0];
    dump->_backSteelHalfSize[1] = dump->_frontSteelHalfSize[1];
    dump->_backSteelHalfSize[2] = additionalSteel[2]/2.0 - dump->_frontSteelHalfSize[2]-0.5;
    dump->_coreAirHalfSize.resize(3);
    dump->_coreAirHalfSize[0] = dump->_coreHalfSize[0] + coreAirGap;
    dump->_coreAirHalfSize[1] = dump->_coreHalfSize[1] + coreAirGap/2.0;
    dump->_coreAirHalfSize[2] = dump->_coreHalfSize[2];

    // The offset of the front shielding coordinates w.r.t. the core, along the dump z
    dump->_coreCenterDistanceToReferencePlane = c.getDouble("protonBeamDump.coreCenterDistanceToReferencePlane");
    dump->_coreCenterDistanceToShieldingFace = 2.0*dump->_frontShieldingHalfSize[2] - dump->_coreHalfSize[2];

    const double frontShieldingOffset = dump->_coreCenterDistanceToShieldingFace - dump->_frontShieldingHalfSize[2];
    dump->_frontShieldingCenterInMu2e[0] = dump->_coreCenterInMu2e[0] + frontShieldingOffset*sin(coreRotY);
    dump->_frontShieldingCenterInMu2e[1] = dump->_shieldingFaceYmin   + dump->_frontShieldingHalfSize[1];
    dump->_frontShieldingCenterInMu2e[2] = dump->_coreCenterInMu2e[2] + frontShieldingOffset*cos(coreRotY);

    dump->_backShieldingHalfSize.resize(3);
    dump->_backShieldingHalfSize[0] = dump->_coreHalfSize[0] + dump->_minCoreShieldingThickness;
    dump->_backShieldingHalfSize[1] = psWall.getYhalfThickness();
    dump->_backShieldingHalfSize[2] = dump->_minCoreShieldingThickness/2;

    const double backShieldingOffset = dump->_backShieldingHalfSize[2] + dump->_coreHalfSize[2];
    dump->_backShieldingCenterInMu2e[0] = dump->_coreCenterInMu2e[0] - backShieldingOffset*sin(coreRotY);
    dump->_backShieldingCenterInMu2e[1] = dump->_shieldingFaceYmin    + dump->_backShieldingHalfSize[1];
    dump->_backShieldingCenterInMu2e[2] = dump->_coreCenterInMu2e[2] - backShieldingOffset*cos(coreRotY);

    dump->_mouthCenterInMu2e = dump->_coreCenterInMu2e + dump->_coreRotationInMu2e
      *CLHEP::Hep3Vector(0,0,frontShieldingOffset + dump->_frontShieldingHalfSize[2] - dump->_mouthHalfSize[2]);

    dump->_neutronCaveCenterInMu2e = dump->_coreCenterInMu2e + dump->_coreRotationInMu2e
      *CLHEP::Hep3Vector(0,0, frontShieldingOffset + dump->_frontShieldingHalfSize[2] - 2*dump->_mouthHalfSize[2] - dump->_neutronCaveHalfSize[2]);

    dump->_coreAirCenterInMu2e = dump->_coreCenterInMu2e + dump->_coreRotationInMu2e * CLHEP::Hep3Vector(0,-coreAirGap/2.0,0);

    dump->_frontSteelCenterInMu2e = dump->_coreCenterInMu2e + dump->_coreRotationInMu2e
      *CLHEP::Hep3Vector(0,     dump->_coreAirHalfSize[1] + dump->_frontSteelHalfSize[1],
                                dump->_coreHalfSize[2] - dump->_frontSteelHalfSize[2]);

    dump->_backSteelCenterInMu2e = dump->_coreCenterInMu2e + dump->_coreRotationInMu2e
      *CLHEP::Hep3Vector(0,     dump->_coreAirHalfSize[1] + dump->_backSteelHalfSize[1],
                                dump->_coreHalfSize[2] - 2.0*dump->_frontSteelHalfSize[2] - dump->_backSteelHalfSize[2]);

    //----------------------------------------------------------------
    // Shielding face coordinates
    dump->_shieldingFaceXmin = dump->_frontShieldingCenterInMu2e[0]
      + dump->_frontShieldingHalfSize[2] * sin(coreRotY)
      - dump->_frontShieldingHalfSize[0] * cos(coreRotY)
      ;

    dump->_shieldingFaceXmax = dump->_frontShieldingCenterInMu2e[0]
      + dump->_frontShieldingHalfSize[2] * sin(coreRotY)
      + dump->_frontShieldingHalfSize[0] * cos(coreRotY)
      ;

    dump->_shieldingFaceZatXmin = dump->_frontShieldingCenterInMu2e[2]
      + dump->_frontShieldingHalfSize[2] * cos(coreRotY)
      + dump->_frontShieldingHalfSize[0] * sin(coreRotY)
      ;

    dump->_shieldingFaceZatXmax = dump->_frontShieldingCenterInMu2e[2]
      + dump->_frontShieldingHalfSize[2] * cos(coreRotY)
      - dump->_frontShieldingHalfSize[0] * sin(coreRotY)
      ;

    //----------------------------------------------------------------
    // For the dump shielding we define an extruded solid that extends
    // the side faces of the nominal dump horizontally to the enclosure
    // walls and vertically to the ceiling and floor. In this way we
    // accommodate the slight diference between the enclosure
    // orientation and the nominal dump rotation.  There is a matching
    // solid that fills the area behind the dump core that ends
    // vertically at the floor level of the extinction monitor room

    // Create deltas (of 1.1 mm, with rotation) to avoid overlaps
    double dW = 1.1;
    double ccW = -cos(coreRotY)*dW;
    const CLHEP::Hep2Vector deltaNE( -ccW, -sin(coreRotY)*dW );
    const CLHEP::Hep2Vector deltaNW(  cos(coreRotY)*dW, -sin(coreRotY)*dW );
    const CLHEP::Hep2Vector deltaSW(  cos(coreRotY)*dW,  sin(coreRotY)*dW );
    const CLHEP::Hep2Vector deltaSE( -cos(coreRotY)*dW,  sin(coreRotY)*dW );

    // Extract or compute the appropriate front shielding vertices
    const CLHEP::Hep3Vector psAoffset =  psArea.getOffsetFromMu2eOrigin();
    const CLHEP::Hep3Vector psWoffset =  psWall.getOffsetFromMu2eOrigin();
    const CLHEP::Hep3Vector psSoffset = psWallS.getOffsetFromMu2eOrigin();
    const CLHEP::Hep3Vector psNoffset = psWallN.getOffsetFromMu2eOrigin();

    const auto & psaVertices =  psArea.getVertices();
    const auto & pswVertices =  psWall.getVertices();
    const auto & psSVertices = psWallS.getVertices();
    const auto & psNVertices = psWallN.getVertices();

    double zs = dump->_coreCenterInMu2e[2] - dump->_coreHalfSize[2]*cos(coreRotY)
      + dump->_frontShieldingHalfSize[0]*sin(coreRotY);
    double xs = dump->_coreCenterInMu2e[0] - dump->_coreHalfSize[2]*sin(coreRotY)
      - dump->_frontShieldingHalfSize[0]*cos(coreRotY);
    double zn = dump->_coreCenterInMu2e[2] - dump->_coreHalfSize[2]*cos(coreRotY)
      - dump->_frontShieldingHalfSize[0]*sin(coreRotY);
    double xn = dump->_coreCenterInMu2e[0] - dump->_coreHalfSize[2]*sin(coreRotY)
      + dump->_frontShieldingHalfSize[0]*cos(coreRotY);

    dump->_backShieldingOutline.emplace_back(psNVertices[5][0]+psNoffset[2]+deltaNE[0],psNVertices[5][1]+psNoffset[0]+deltaNE[1]);
    dump->_backShieldingOutline.emplace_back(pswVertices[14][0]+psWoffset[2]+deltaNW[0],pswVertices[14][1]+psWoffset[0]+deltaNW[1]);
    dump->_backShieldingOutline.emplace_back(pswVertices[15][0]+psWoffset[2]+deltaSW[0],pswVertices[15][1]+psWoffset[0]+deltaSW[1]);
    dump->_backShieldingOutline.emplace_back(psSVertices[7][0]+psSoffset[2]+deltaSE[0],psSVertices[7][1]+psSoffset[0]+deltaSE[1]);
    dump->_backShieldingOutline.emplace_back(zs,xs);
    dump->_backShieldingOutline.emplace_back(zn,xn);

    dump->_frontShieldingOutline.emplace_back(psNVertices[5][0]+psNoffset[2]+deltaNE[0],psNVertices[5][1]+psNoffset[0]+deltaNE[1]);
    dump->_frontShieldingOutline.emplace_back(zn,xn);
    dump->_frontShieldingOutline.emplace_back(zs,xs);
    dump->_frontShieldingOutline.emplace_back(psSVertices[7][0]+psSoffset[2]+deltaSE[0],psSVertices[7][1]+psSoffset[0]+deltaSE[1]);
    // Temporary workaround while bldg geom is in flux.
    double tempWorkaround = 10.5;
    dump->_frontShieldingOutline.emplace_back(psaVertices[15][0]+psAoffset[2]+deltaSE[0],psaVertices[15][1]+psAoffset[0]+deltaSE[1]+tempWorkaround);
    dump->_frontShieldingOutline.emplace_back(dump->_shieldingFaceZatXmin,dump->_shieldingFaceXmin);
    dump->_frontShieldingOutline.emplace_back(dump->_shieldingFaceZatXmax,dump->_shieldingFaceXmax);
    dump->_frontShieldingOutline.emplace_back(psaVertices[14][0]+psAoffset[2]+deltaNE[0],psaVertices[14][1]+psAoffset[0]+deltaNE[1]);

    //----------------------------------------------------------------
    if(verbose) {

      for( std::size_t i=0; i<dump->_frontShieldingOutline.size(); i++ ) {
        std::cout<<"front "     <<      dump->_frontShieldingOutline[i][0]
                 <<"\t"         <<      dump->_frontShieldingOutline[i][1]      <<"\n";
      }

      for( std::size_t i=0; i<dump->_backShieldingOutline.size(); i++ ) {
        std::cout<<"back "      <<      dump->_backShieldingOutline[i][0]
                 <<"\t"         <<      dump->_backShieldingOutline[i][1]       <<"\n";
      }

      std::cout<<"ProtonBeamDumpMaker"<<": ProtonBeamDump core center in mu2e = "<<dump->_coreCenterInMu2e<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": shieldingFaceXmin = "<<dump->_shieldingFaceXmin
               <<", Xmax = "<<dump->_shieldingFaceXmax<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": shieldingFaceZatXmin = "<<dump->_shieldingFaceZatXmin
               <<", ZatXmax = "<<dump->_shieldingFaceZatXmax<<std::endl;
      std::cout<<"Front Steel Center in Mu2e:  " << dump->_frontSteelCenterInMu2e[0] << ", " << dump->_frontSteelCenterInMu2e[1] << ", " << dump->_frontSteelCenterInMu2e[2] << std::endl;

    }

    return dump;

  } // make()

} // namespace mu2e
