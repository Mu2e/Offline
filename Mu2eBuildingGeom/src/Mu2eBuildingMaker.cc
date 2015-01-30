//
// Construct Mu2e building
//
// $Id: Mu2eBuildingMaker.cc,v 1.17 2013/09/27 17:19:39 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/09/27 17:19:39 $
//
// Original author: Andrei Gaponenko

// Mu2e include files
#include "ConfigTools/inc/SimpleConfig.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/TransportSolenoid.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuildingMaker.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

// C++ include files
#include <iostream>

// Framework include files
#include "cetlib/exception.h"

// CLHEP include files
#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {

  std::unique_ptr<Mu2eBuilding> Mu2eBuildingMaker::make(const SimpleConfig& c,
                                                        const BuildingBasics& basics,
                                                        const Beamline& bl,
                                                        const ProtonBeamDump& dump)
  {
    std::unique_ptr<Mu2eBuilding> b (new Mu2eBuilding(basics));

    using CLHEP::Hep2Vector;

    b->_hallInsideXmin = c.getDouble("hall.insideXmin");
    b->_hallInsideXmax = c.getDouble("hall.insideXmax");
    b->_hallInsideZmax = c.getDouble("hall.insideZmax");

    b->_hallInsideXDSCorner = c.getDouble("hall.insideXDSCorner");
    b->_hallInsideZDSCorner = c.getDouble("hall.insideZDSCorner");

    b->_hallInsideXPSCorner = c.getDouble("hall.insideXPSCorner");
    b->_hallInsideZPSCorner = c.getDouble("hall.insideZPSCorner");

    b->_hallInsideZExtMonUCIWall = c.getDouble("hall.insideZExtMonUCIWall");

    b->_hallWallThickness = c.getDouble("hall.wallThick");
    b->_hallWallExtMonUCIThickness = c.getDouble("hall.wallExtMonUCIThick");

    // Origin used to construct the MECO detector.
    // Magic number to fix:

    // FIXME:::: Use numbers from DetectorSystem! (At least SolenoidOffset?)
    b->_relicMECOOriginInMu2e = CLHEP::Hep3Vector( -3904., 0., 12000.);

    //----------------
    //  Computed stuff

    b->_hallInsideZPStoBeamDumpCorner = dump.frontShieldingCenterInMu2e()[2]
      + (b->_hallInsideXPSCorner - dump.frontShieldingCenterInMu2e()[0])/tan(dump.coreRotY())
      + dump.frontShieldingHalfSize()[0]/sin(dump.coreRotY())
      ;

    //----------------
    //  Check the assumptions used to construct G4 geometry downstream

    if(b->_hallInsideZPStoBeamDumpCorner >= b->_hallInsideZPSCorner) {
      throw cet::exception("GEOM")<<"Mu2eBuildingMaker: can not satisfy (hall.insideXPSCorner, hall.insideZPSCorner)\n";
    }

    if(dump.shieldingFaceZatXmax() <= b->_hallInsideZExtMonUCIWall) {
      throw cet::exception("GEOM")<<"Mu2eBuildingMaker: hallInsideZExtMonUCIWall is too large - conflicts with ProtonBeamDump\n";
    }

    //----------------------------------------------------------------
    // Fragment 1
    b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax(), dump.shieldingFaceZatXmax()));
    b->_hallInsideOutline.push_back(b->_concreteOuterOutline1.back());

    const double concreteZmin = b->hallInsideZExtMonUCIWall() - b->hallWallThickness();
    if(b->hallWallExtMonUCIThickness() < (dump.shieldingFaceZatXmax() - concreteZmin)*tan(dump.coreRotY())) {
      // two-point case
      b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax() - b->hallWallExtMonUCIThickness(),
                                                     dump.shieldingFaceZatXmax() - b->hallWallExtMonUCIThickness()
                                                     /tan(dump.coreRotY())));

      b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax() - b->hallWallExtMonUCIThickness(),
                                                     concreteZmin));
    }
    else { // one-point case
      const double dz = dump.shieldingFaceZatXmax() - concreteZmin;
      b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax() - dz*tan(dump.coreRotY()),
                                                     concreteZmin));
    }
    b->_hallInsideOutline.push_back(Hep2Vector(dump.shieldingFaceXmax(), concreteZmin + b->hallWallThickness()));

    b->_concreteOuterOutline1.push_back(Hep2Vector(b->hallInsideXmax() + b->hallWallThickness(), concreteZmin));
    b->_hallInsideOutline      .push_back(Hep2Vector( b->hallInsideXmax(), concreteZmin + b->hallWallThickness() ) );
    b->_concreteInnerOutlineLowerExt.push_back(Hep2Vector( b->hallInsideXmax(), 
                                                           b->hallInsideZPStoBeamDumpCorner() + b->hallWallThickness()*tan(dump.coreRotY()/2) ) );

    //----------------------------------------------------------------
    // Fragment 2
    b->_concreteOuterOutline2  .push_back(Hep2Vector( b->hallInsideXmax()+b->hallWallThickness(), b->hallInsideZDSCorner()+b->hallWallThickness()) );
    b->_concreteOuterOutlineExt.push_back(Hep2Vector( b->hallInsideXmax()+b->hallWallThickness(), b->hallInsideZDSCorner()+b->hallWallThickness()) );

    b->_hallInsideOutline      .push_back(Hep2Vector(b->hallInsideXmax(), b->hallInsideZDSCorner()));
    b->_concreteInnerOutlineLowerExt.push_back(Hep2Vector(b->hallInsideXmax(), b->hallInsideZDSCorner()));

    //----------------
    b->_concreteOuterOutline2  .push_back(Hep2Vector(b->hallInsideXDSCorner()+b->hallWallThickness(), b->hallInsideZDSCorner()+b->hallWallThickness()) );
    b->_hallInsideOutline      .push_back(Hep2Vector(b->hallInsideXDSCorner()                       , b->hallInsideZDSCorner())                        );

    // Get TS shielding L-above x-extent
    // - hack for now
    const double xExtentOfLowerCeiling = 
      c.getDouble("ExtShieldUpstream.XExtentHack");
    b->_concreteInnerOutlineLowerExt.push_back(Hep2Vector( xExtentOfLowerCeiling, b->hallInsideZDSCorner() ) );

    //----------------
    b->_concreteOuterOutline2.push_back(Hep2Vector(b->hallInsideXDSCorner()+b->hallWallThickness(),b->hallInsideZmax() + b->hallWallThickness()));
    b->_hallInsideOutline    .push_back(Hep2Vector(b->hallInsideXDSCorner()                       , b->hallInsideZmax() )                       );

    //----------------
    b->_concreteOuterOutline2  .push_back(Hep2Vector(b->hallInsideXmin()-b->hallWallThickness(), b->hallInsideZmax()      + b->hallWallThickness()) );
    b->_concreteOuterOutlineExt.push_back(Hep2Vector(b->hallInsideXmin()-b->hallWallThickness(), b->hallInsideZDSCorner() + b->hallWallThickness()) );

    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXmin(), b->hallInsideZmax()));

    //----------------------------------------------------------------
    // Fragment 3

    b->_concreteOuterOutline3  .push_back(Hep2Vector(b->hallInsideXmin()-b->hallWallThickness(), b->hallInsideZPSCorner()-b->hallWallThickness()));
    b->_hallInsideOutline      .push_back(Hep2Vector(b->hallInsideXmin()                       , b->hallInsideZPSCorner())                       );
    b->_concreteInnerOutlineLowerExt.push_back(Hep2Vector( xExtentOfLowerCeiling, b->hallInsideZPSCorner() ) );

    //----------------
    b->_concreteOuterOutline3.push_back(Hep2Vector(b->hallInsideXPSCorner() - b->hallWallThickness(),
                                                   b->hallInsideZPSCorner() - b->hallWallThickness()));

    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXPSCorner(), b->hallInsideZPSCorner()));
    b->_concreteInnerOutlineLowerExt.push_back(Hep2Vector(b->hallInsideXPSCorner(), b->hallInsideZPSCorner() ) );

    //----------------
    b->_concreteOuterOutline3.push_back(Hep2Vector(b->hallInsideXPSCorner() - b->hallWallThickness(),
                                                   b->hallInsideZPStoBeamDumpCorner() + b->hallWallThickness()*tan(dump.coreRotY()/2)
                                                   ));
    
    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXPSCorner(), b->hallInsideZPStoBeamDumpCorner()));
    
    //----------------
    b->_concreteOuterOutline3.push_back(Hep2Vector(dump.shieldingFaceXmin()
                                                   - b->hallWallThickness()*cos(dump.coreRotY())
                                                   ,
                                                   dump.shieldingFaceZatXmin()
                                                   + b->hallWallThickness()*sin(dump.coreRotY())
                                                   ));

    b->_hallInsideOutline.push_back(Hep2Vector(dump.shieldingFaceXmin(), dump.shieldingFaceZatXmin()));
    b->_concreteInnerOutlineLowerExt.push_back(Hep2Vector( b->hallInsideXPSCorner(), 
                                                           b->hallInsideZPStoBeamDumpCorner() + b->hallWallThickness()*tan(dump.coreRotY()/2) ) );
    //----------------------------------------------------------------
    // Beamline slabs
    // - x-extent estimated as edge of TS1 cryostat
    StraightSection const * ts1Cryo = bl.getTS().getTSCryo<StraightSection>( TransportSolenoid::TSRegion::TS1,
                                                                             TransportSolenoid::TSRadialPart::OUT);
    b->_xPosOfSlabEnd = ts1Cryo->getGlobal().x() - ts1Cryo->rOut();

    std::vector<double> xOffset, yThicknesses;
    b->_nBeamlineSlabs = c.getInt( "hall.beamlineSlabs.nSlabs" );
    const double globalOffset = c.getDouble( "hall.beamlineSlabs.xOffsetGlobal" );
    c.getVectorDouble( "hall.beamlineSlabs.xOffset"     , xOffset     , b->_nBeamlineSlabs );
    c.getVectorDouble( "hall.beamlineSlabs.yThicknesses", yThicknesses, b->_nBeamlineSlabs );
    const double zWidth       = c.getDouble( "hall.beamlineSlabs.zWidth" );

    for ( std::size_t iSlab = 0; iSlab < b->_nBeamlineSlabs; iSlab++ ) {
      b->_concreteBeamlineSlabs.emplace_back( 0.5*(b->xPosOfSlabEnd() - xOffset.at(iSlab) - globalOffset ),
                                              0.5*yThicknesses.at(iSlab),
                                              0.5*zWidth );
    }
    
    //----------------------------------------------------------------
    const int diagLevel = c.getInt("world.verbosityLevel", 0);
    if(diagLevel > 0) {
      std::cout << __func__ << " relicMECOOriginInMu2e : " <<  b->_relicMECOOriginInMu2e  << std::endl;
    }

    return b;
  }
}
