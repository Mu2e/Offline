// Andrei Gaponenko, 2011

#include "ProtonBeamDumpGeom/inc/ProtonBeamDumpMaker.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

namespace mu2e {

  //================================================================
  std::auto_ptr<ProtonBeamDump> ProtonBeamDumpMaker::make(const SimpleConfig& c,
                                                          double frontShieldingYmin,
                                                          double frontShieldingYmax)
  {

    std::auto_ptr<ProtonBeamDump> dump(new ProtonBeamDump());

    int verbose = c.getInt("protonBeamDump.verbosityLevel", 0);

    // position
    dump->_coreCenterInMu2e = c.getHep3Vector("protonBeamDump.coreCenterInMu2e");
    const double coreRotY = dump->_coreRotY = c.getDouble("protonBeamDump.coreRotY") * CLHEP::degree;
    dump->_coreRotationInMu2e.rotateY(coreRotY);

    c.getVectorDouble("protonBeamDump.coreHalfSize", dump->_coreHalfSize, 3);
    c.getVectorDouble("protonBeamDump.neutronCaveHalfSize", dump->_neutronCaveHalfSize, 3);
    c.getVectorDouble("protonBeamDump.mouthHalfSize", dump->_mouthHalfSize, 3);
    dump->_minCoreShieldingThickness = c.getDouble("protonBeamDump.minCoreShieldingThickness");

    // Compute the overall size
    dump->_frontShieldingHalfSize.resize(3);
    dump->_frontShieldingHalfSize[0] = dump->_coreHalfSize[0] + dump->_minCoreShieldingThickness;
    dump->_frontShieldingHalfSize[1] = (frontShieldingYmax - frontShieldingYmin)/2;
    dump->_frontShieldingHalfSize[2] = c.getDouble("protonBeamDump.frontShieldingThickness")/2;

    // shorthand notation
    const std::vector<double>& frontShieldingHalfSize = dump->_frontShieldingHalfSize;

    if(verbose) {
      std::cout<<"ProtonBeamDumpMaker"<<": ProtonBeamDump frontShielding half size = ";
      std::copy(frontShieldingHalfSize.begin(), frontShieldingHalfSize.end(), std::ostream_iterator<double>(std::cout, ", "));
      std::cout<<std::endl;
    }

    dump->_coreCenterDistanceToShieldingFace = dump->_coreHalfSize[2] + 2*dump->_neutronCaveHalfSize[2] + 2*dump->_mouthHalfSize[2];

    // The offset of the front shielding box center w.r.t. the core, along the dump z
    const double frontShieldingOffset = dump->_coreCenterDistanceToShieldingFace - frontShieldingHalfSize[2];

    dump->_frontShieldingCenterInMu2e[0] = dump->_coreCenterInMu2e[0] + frontShieldingOffset*sin(coreRotY);
    dump->_frontShieldingCenterInMu2e[1] = (frontShieldingYmax + frontShieldingYmin)/2;
    dump->_frontShieldingCenterInMu2e[2] = dump->_coreCenterInMu2e[2] + frontShieldingOffset*cos(coreRotY);

    const CLHEP::Hep3Vector& frontShieldingCenterInMu2e = dump->_frontShieldingCenterInMu2e;

    //----------------------------------------------------------------
    // frontShieldingThickness does not give enough concrete at the back of the core.
    // We need to create an additional box there.
    dump->_backShieldingHalfSize.resize(3);

    dump->_backShieldingHalfSize[0] = dump->_coreHalfSize[0] + dump->_minCoreShieldingThickness;
    dump->_backShieldingHalfSize[1] = dump->_coreHalfSize[1] + dump->_minCoreShieldingThickness;

    // Not mu2e z, but along dump axis
    const double backzmin = -(dump->_minCoreShieldingThickness + dump->_coreHalfSize[2]);
    const double backzmax =  +frontShieldingOffset - dump->_frontShieldingHalfSize[2];

    dump->_backShieldingHalfSize[2] = (backzmax - backzmin)/2;
    dump->_backShieldingCenterInMu2e = dump->_coreCenterInMu2e +
      dump->_coreRotationInMu2e * CLHEP::Hep3Vector(0, 0, (backzmax + backzmin)/2);

    dump->_mouthCenterInMu2e = dump->_coreCenterInMu2e + dump->_coreRotationInMu2e
      *CLHEP::Hep3Vector(0,0,frontShieldingOffset + dump->_frontShieldingHalfSize[2] - dump->_mouthHalfSize[2]);

    dump->_neutronCaveCenterInMu2e = dump->_coreCenterInMu2e + dump->_coreRotationInMu2e
      *CLHEP::Hep3Vector(0,0, frontShieldingOffset + dump->_frontShieldingHalfSize[2] - 2*dump->_mouthHalfSize[2] - dump->_neutronCaveHalfSize[2]);

    //----------------------------------------------------------------
    // Shielding face coordinates
    dump->_shieldingFaceXmin = frontShieldingCenterInMu2e[0]
      + frontShieldingHalfSize[2] * sin(coreRotY)
      - frontShieldingHalfSize[0] * cos(coreRotY)
      ;

    dump->_shieldingFaceXmax = frontShieldingCenterInMu2e[0]
      + frontShieldingHalfSize[2] * sin(coreRotY)
      + frontShieldingHalfSize[0] * cos(coreRotY)
      ;

    dump->_shieldingFaceZatXmin = frontShieldingCenterInMu2e[2]
      + frontShieldingHalfSize[2] * cos(coreRotY)
      + frontShieldingHalfSize[0] * sin(coreRotY)
      ;

    dump->_shieldingFaceZatXmax = frontShieldingCenterInMu2e[2]
      + frontShieldingHalfSize[2] * cos(coreRotY)
      - frontShieldingHalfSize[0] * sin(coreRotY)
      ;

    //----------------------------------------------------------------
    if(verbose) {
      std::cout<<"ProtonBeamDumpMaker"<<": ProtonBeamDump core center in mu2e = "<<dump->_coreCenterInMu2e<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": shieldingFaceXmin = "<<dump->_shieldingFaceXmin
               <<", Xmax = "<<dump->_shieldingFaceXmax<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": shieldingFaceZatXmin = "<<dump->_shieldingFaceZatXmin
               <<", ZatXmax = "<<dump->_shieldingFaceZatXmax<<std::endl;
    }

    return dump;

  } // make()

} // namespace mu2e
