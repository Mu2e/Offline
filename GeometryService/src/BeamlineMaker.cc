//
// Construct and return an Beamline.
//
// Original author Peter Shanahan
//                 Kyle Knoepfel (significant updates)
//

// C++ includes
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <utility>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/TransportSolenoid.hh"
#include "BeamlineGeom/inc/Collimator_TS1.hh"
#include "BeamlineGeom/inc/Collimator_TS3.hh"
#include "BeamlineGeom/inc/Collimator_TS5.hh"
#include "GeometryService/inc/BeamlineMaker.hh"
#include "GeomPrimitives/inc/Torus.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

#ifndef __CINT__

namespace mu2e {

  std::unique_ptr<Beamline> BeamlineMaker::make(const SimpleConfig& c) {
    std::unique_ptr<Beamline> res ( new Beamline() );
    BuildBeamline     (c, res.get()  );
    BuildTSCryostat   (c, res.get()  );
    BuildTSCoils      (c, res.get()  );
    BuildTSCollimators(c, &res->_ts  );
    BuildTSVacua      (c, &res->_ts  );
    BuildTSCAs        (c, *(res.get()) );
    BuildTSPolyLining (c, &res->_ts  );
    BuildPbarWindow   (c, &res->_ts  );
    return res;
  }

  void BeamlineMaker::BuildBeamline(const SimpleConfig& c, Beamline* bl) {
    bl->_solenoidOffset = c.getDouble("mu2e.solenoidOffset");
  }

  void BeamlineMaker::BuildTSCryostat(const SimpleConfig& c, Beamline* bl ){

    TransportSolenoid & ts = bl->_ts;
    ts._rTorus   = c.getDouble("ts.rTorus",0.);
    ts._rVac     = c.getDouble("ts.rVac",0.);
    ts._ts1RVac  = c.getDouble("ts.ts1.rVac", ts._rVac);
    ts._ts2RVac  = c.getDouble("ts.ts2.rVac", ts._rVac);
    ts._ts3RVac  = c.getDouble("ts.ts3.rVac", ts._rVac);
    ts._ts4RVac  = c.getDouble("ts.ts4.rVac", ts._rVac);
    ts._ts5RVac  = c.getDouble("ts.ts5.rVac", ts._rVac);

    ts._material = c.getString("ts.materialName");
    ts._downstreamVacuumMaterial = c.getString("ts.downstreamVacuumMaterialName");
    ts._upstreamVacuumMaterial   = c.getString("ts.upstreamVacuumMaterialName");
    ts._thermalShieldMLIMaterial = c.getString("ts.thermalshield.mli.material", "");
    ts._thermalShieldMidMaterial = c.getString("ts.thermalshield.mid.material", "");

    // Parameters for rings (David Norvil Brown added April 1 2015)
    ts._rInRingSide = c.getDouble("ts.rInRingSide");
    ts._rOutRingSide = c.getDouble("ts.rOutRingSide");
    ts._thickRingSide = c.getDouble("ts.thickRingSide");
    ts._rInRing = c.getDouble("ts.rInRing");
    ts._rOutRing = c.getDouble("ts.rOutRing");
    ts._lengthRing = c.getDouble("ts.lengthRing");
    ts._RingMaterial = c.getString("ts.RingMaterialType");

    int nRing = c.getInt("ts.nRing");
    c.getVectorDouble("ts.xRing", ts._xRing, nRing);
    c.getVectorDouble("ts.yRing", ts._yRing, nRing);
    c.getVectorDouble("ts.zRing", ts._zRing, nRing);
    c.getVectorDouble("ts.thetaRing",ts._thetaRing, nRing);

    // - end wall parameters
    ts._build_endWallD2 = c.getBool("ts.tsDendWall2.build", false);  //default to no second component
    ts._rIn_endWallU1 = c.getDouble("ts.tsUendWall1.rIn",c.getDouble("ts.ts1in.rOut") );
    ts._rIn_endWallU2 = c.getDouble("ts.tsUendWall2.rIn",0.);
    ts._rIn_endWallD  = c.getDouble("ts.tsDendWall.rIn", c.getDouble("ts.ts5in.rOut") );
    int nEndWallD2Planes = 0;
    if(ts._build_endWallD2) {
      nEndWallD2Planes = c.getInt("ts.tsDendWall2.planes");
      c.getVectorDouble("ts.tsDendWall2.rIn", ts._rIn_endWallD2, nEndWallD2Planes);
    }

    ts._rOut_endWallU1 = c.getDouble("ts.tsUendWall1.rOut",c.getDouble("ts.ts1out.rIn") );
    ts._rOut_endWallU2 = c.getDouble("ts.tsUendWall2.rOut",0.);
    ts._rOut_endWallD  = c.getDouble("ts.tsDendWall.rOut", c.getDouble("ts.ts5out.rIn") );
    if(ts._build_endWallD2) c.getVectorDouble("ts.tsDendWall2.rOut", ts._rOut_endWallD2, nEndWallD2Planes);

    ts._halfLength_endWallU1 = c.getDouble("ts.tsUendWall1.halfLength");
    ts._halfLength_endWallU2 = c.getDouble("ts.tsUendWall2.halfLength");
    ts._halfLength_endWallD  = c.getDouble("ts.tsDendWall.halfLength" );
    if(ts._build_endWallD2) {
      c.getVectorDouble("ts.tsDendWall2.z", ts._z_endWallD2, nEndWallD2Planes);
      ts._halfLength_endWallD2 = ts._z_endWallD2[nEndWallD2Planes-1]/2.;
    }

    double ts1HalfLength = c.getDouble("ts.ts1.halfLength");
    double ts1CALengthDiff = c.getDouble("ts.ts1.caLengthDiff", 0.0); //distance from end wall it ends
    double ts3HalfLength = bl->solenoidOffset() - ts.torusRadius();
    double ts5HalfLength = c.getDouble("ts.ts5.halfLength");
    double ts5CALengthDiff = c.getDouble("ts.ts5.caLengthDiff", 0.0); //distance from end wall it ends

    double ts1zOffset    = (-ts.torusRadius()-ts1HalfLength );
    double ts5zOffset    = ( ts.torusRadius()+ts5HalfLength );

    // Typedefs for easier use
    typedef TransportSolenoid::TSRegion::enum_type     tsReg_enum;
    typedef TransportSolenoid::TSRadialPart::enum_type tsRad_enum;

    // Bookkeeping index
    // - note that this index mapping is unique as straight section
    // and torus section numbering are always incremented by 2, not 1
    auto MapIndex = [](tsReg_enum iTS, tsRad_enum iRAD) -> unsigned {
      return static_cast<unsigned>(iTS) + static_cast<unsigned>(iRAD);
    };

    // Cryo map parameters - straight sections

    const std::map<unsigned,Tube> straightSectionParams = {
      // Inner straight sections
      { MapIndex(tsReg_enum::TS1,tsRad_enum::IN ), Tube(ts.ts1InnerRadius(),
                                                        c.getDouble("ts.ts1in.rOut"),
                                                        ts1HalfLength,
                                                        CLHEP::Hep3Vector( bl->solenoidOffset(),0.0,ts1zOffset) ) },
      { MapIndex(tsReg_enum::TS3,tsRad_enum::IN ), Tube(ts.ts3InnerRadius(),
                                                        c.getDouble("ts.ts3in.rOut"),
                                                        ts3HalfLength,
                                                        CLHEP::Hep3Vector(),
                                                        CLHEP::HepRotation(CLHEP::HepRotationY(90.0*CLHEP::degree)) ) },
      { MapIndex(tsReg_enum::TS5,tsRad_enum::IN ), Tube(ts.ts5InnerRadius(),
                                                        c.getDouble("ts.ts5in.rOut"),
                                                        ts5HalfLength,
                                                        CLHEP::Hep3Vector(-bl->solenoidOffset(),0.0,ts5zOffset) ) },
      // Outer straight sections
      { MapIndex(tsReg_enum::TS1,tsRad_enum::OUT), Tube(c.getDouble("ts.ts1out.rIn"),
                                                        c.getDouble("ts.ts1out.rOut"),
                                                        ts1HalfLength - ts.endWallU2_halfLength(),
                                                        CLHEP::Hep3Vector( bl->solenoidOffset(),0.0,ts1zOffset-ts.endWallU2_halfLength() ) ) },
      { MapIndex(tsReg_enum::TS3,tsRad_enum::OUT), Tube(c.getDouble("ts.ts3out.rIn"),
                                                        c.getDouble("ts.ts3out.rOut"),
                                                        ts3HalfLength,
                                                        CLHEP::Hep3Vector(),
                                                        CLHEP::HepRotation(CLHEP::HepRotationY(90.0*CLHEP::degree)) ) },
      { MapIndex(tsReg_enum::TS5,tsRad_enum::OUT), Tube(c.getDouble("ts.ts5out.rIn"),
                                                        c.getDouble("ts.ts5out.rOut"),
                                                        ts5HalfLength,
                                                        CLHEP::Hep3Vector(-bl->solenoidOffset(),0.0,ts5zOffset) ) },
    };


    // Cryo map parameters - torus sections

    const std::map<unsigned,Torus> torusSectionParams = {
      // Inner torus sections
      { MapIndex(tsReg_enum::TS2,tsRad_enum::IN ), Torus( ts.torusRadius(),
                                                          ts.ts2InnerRadius(),
                                                          c.getDouble("ts.ts2in.rOut"),
                                                          1.5*CLHEP::pi, CLHEP::halfpi,
                                                          CLHEP::Hep3Vector( ts3HalfLength, 0.,-ts.torusRadius() ),
                                                          CLHEP::HepRotation(CLHEP::HepRotationX(90.0*CLHEP::degree))) },
      { MapIndex(tsReg_enum::TS4,tsRad_enum::IN ), Torus( ts.torusRadius(),
                                                          ts.ts4InnerRadius(),
                                                          c.getDouble("ts.ts4in.rOut"),
                                                          CLHEP::halfpi, CLHEP::halfpi,
                                                          CLHEP::Hep3Vector( -ts3HalfLength, 0.,ts.torusRadius() ),
                                                          CLHEP::HepRotation(CLHEP::HepRotationX(90.0*CLHEP::degree))) },
      // Outer torus sections
      { MapIndex(tsReg_enum::TS2,tsRad_enum::OUT), Torus( ts.torusRadius(),
                                                          c.getDouble("ts.ts2out.rIn"),
                                                          c.getDouble("ts.ts2out.rOut"),
                                                          1.5*CLHEP::pi, CLHEP::halfpi,
                                                          CLHEP::Hep3Vector( ts3HalfLength, 0.,-ts.torusRadius() ),
                                                          CLHEP::HepRotation(CLHEP::HepRotationX(90.0*CLHEP::degree))) },
      { MapIndex(tsReg_enum::TS4,tsRad_enum::OUT), Torus( ts.torusRadius(),
                                                          c.getDouble("ts.ts4out.rIn"),
                                                          c.getDouble("ts.ts4out.rOut"),
                                                          CLHEP::halfpi, CLHEP::halfpi,
                                                          CLHEP::Hep3Vector( -ts3HalfLength, 0.,ts.torusRadius() ),
                                                          CLHEP::HepRotation(CLHEP::HepRotationX(90.0*CLHEP::degree))) },
    };

    // Set cryo map - straight sections
    for ( unsigned iTS = tsReg_enum::TS1 ; iTS <= tsReg_enum::TS5 ; iTS+=2 ) {
      auto its    = (tsReg_enum)iTS;

      for ( unsigned iRAD = tsRad_enum::IN ; iRAD <= tsRad_enum::OUT ; ++iRAD ) {
        auto irad = (tsRad_enum)iRAD;
        auto straightParam = straightSectionParams.find( MapIndex(its,irad) );
	auto strsec = new StraightSection ( straightParam->second );
	if(its == tsReg_enum::TS1)
	  strsec->setLengthDiff(ts1CALengthDiff); //end the cold mass before the endwall
	if(its == tsReg_enum::TS5)
	  strsec->setLengthDiff(ts5CALengthDiff); //end the cold mass before the endwall
        ts._cryoMap[its][irad] = std::unique_ptr<TSSection>( strsec );

      }
    }

    // Set cryo map - torus sections
    for ( unsigned iTS = tsReg_enum::TS2 ; iTS <= tsReg_enum::TS4 ; iTS+=2 ) {
      auto its = (tsReg_enum)iTS;

      for ( unsigned iRAD = tsRad_enum::IN ; iRAD <= tsRad_enum::OUT ; ++iRAD ) {
        auto irad = (tsRad_enum)iRAD;
        auto torusParam = torusSectionParams.find( MapIndex(its,irad ) );

        ts._cryoMap[its][irad] = std::unique_ptr<TSSection>( new TorusSection ( torusParam->second ) );

      }
    }

    if(!c.getBool("ts.thermalshield.build", false)) return;

    // Thermal shielding map parameters - straight sections
    double smallOffset = 1.e-3; //prevent small overlaps in G4 volumes
    const std::map<unsigned,Tube> straightSectionThermalParams = {
      // Inner straight sections
      { MapIndex(tsReg_enum::TS1,tsRad_enum::IN ), Tube(c.getDouble("ts.ts1in.thermalshield.rIn"),
                                                        c.getDouble("ts.ts1in.thermalshield.rOut"),
                                                        ts1HalfLength - ts.endWallU1_halfLength() - smallOffset,
                                                        CLHEP::Hep3Vector( bl->solenoidOffset(),0.0,ts1zOffset + ts.endWallU1_halfLength() + smallOffset) ) },
      { MapIndex(tsReg_enum::TS3,tsRad_enum::IN ), Tube(c.getDouble("ts.ts3in.thermalshield.rIn"),
                                                        c.getDouble("ts.ts3in.thermalshield.rOut"),
                                                        ts3HalfLength,
                                                        CLHEP::Hep3Vector(),
                                                        CLHEP::HepRotation(CLHEP::HepRotationY(90.0*CLHEP::degree)) ) },
      { MapIndex(tsReg_enum::TS5,tsRad_enum::IN ), Tube(c.getDouble("ts.ts5in.thermalshield.rIn"),
                                                        c.getDouble("ts.ts5in.thermalshield.rOut"),
                                                        ts5HalfLength - ts.endWallD_halfLength() - ts.endWallD2_halfLength() - smallOffset,
                                                        CLHEP::Hep3Vector(-bl->solenoidOffset(),0.0,ts5zOffset - ts.endWallD_halfLength()
									  - ts.endWallD2_halfLength() - smallOffset)) },
      // Outer straight sections
      { MapIndex(tsReg_enum::TS1,tsRad_enum::OUT), Tube(c.getDouble("ts.ts1out.thermalshield.rIn"),
                                                        c.getDouble("ts.ts1out.thermalshield.rOut"),
                                                        ts1HalfLength - ts.endWallU1_halfLength() - smallOffset,
                                                        CLHEP::Hep3Vector( bl->solenoidOffset(),0.0,ts1zOffset + ts.endWallU1_halfLength() + smallOffset) ) },
      { MapIndex(tsReg_enum::TS3,tsRad_enum::OUT), Tube(c.getDouble("ts.ts3out.thermalshield.rIn"),
                                                        c.getDouble("ts.ts3out.thermalshield.rOut"),
                                                        ts3HalfLength,
                                                        CLHEP::Hep3Vector(),
                                                        CLHEP::HepRotation(CLHEP::HepRotationY(90.0*CLHEP::degree)) ) },
      { MapIndex(tsReg_enum::TS5,tsRad_enum::OUT), Tube(c.getDouble("ts.ts5out.thermalshield.rIn"),
                                                        c.getDouble("ts.ts5out.thermalshield.rOut"),
                                                        ts5HalfLength - ts.endWallD_halfLength() - ts.endWallD2_halfLength() - smallOffset,
                                                        CLHEP::Hep3Vector(-bl->solenoidOffset(),0.0,ts5zOffset - ts.endWallD_halfLength()
									  - ts.endWallD2_halfLength() - smallOffset) ) },
    };


    // Thermal shield map parameters - torus sections

    const std::map<unsigned,Torus> torusSectionThermalParams = {
      // Inner torus sections
      { MapIndex(tsReg_enum::TS2,tsRad_enum::IN ), Torus( ts.torusRadius(),
                                                          c.getDouble("ts.ts2in.thermalshield.rIn"),
							  c.getDouble("ts.ts2in.thermalshield.rOut"),
                                                          1.5*CLHEP::pi, CLHEP::halfpi,
                                                          CLHEP::Hep3Vector( ts3HalfLength, 0.,-ts.torusRadius() ),
                                                          CLHEP::HepRotation(CLHEP::HepRotationX(90.0*CLHEP::degree))) },
      { MapIndex(tsReg_enum::TS4,tsRad_enum::IN ), Torus( ts.torusRadius(),
                                                          c.getDouble("ts.ts4in.thermalshield.rIn"),
							  c.getDouble("ts.ts4in.thermalshield.rOut"),
                                                          CLHEP::halfpi, CLHEP::halfpi,
                                                          CLHEP::Hep3Vector( -ts3HalfLength, 0.,ts.torusRadius() ),
                                                          CLHEP::HepRotation(CLHEP::HepRotationX(90.0*CLHEP::degree))) },
      // Outer torus sections
      { MapIndex(tsReg_enum::TS2,tsRad_enum::OUT), Torus( ts.torusRadius(),
                                                          c.getDouble("ts.ts2out.thermalshield.rIn"),
							  c.getDouble("ts.ts2out.thermalshield.rOut"),
                                                          1.5*CLHEP::pi, CLHEP::halfpi,
                                                          CLHEP::Hep3Vector( ts3HalfLength, 0.,-ts.torusRadius() ),
                                                          CLHEP::HepRotation(CLHEP::HepRotationX(90.0*CLHEP::degree))) },
      { MapIndex(tsReg_enum::TS4,tsRad_enum::OUT), Torus( ts.torusRadius(),
                                                          c.getDouble("ts.ts4out.thermalshield.rIn"),
							  c.getDouble("ts.ts4out.thermalshield.rOut"),
                                                          CLHEP::halfpi, CLHEP::halfpi,
                                                          CLHEP::Hep3Vector( -ts3HalfLength, 0.,ts.torusRadius() ),
                                                          CLHEP::HepRotation(CLHEP::HepRotationX(90.0*CLHEP::degree))) },
    };

    // Set thermal shield map - straight sections
    for ( unsigned iTS = tsReg_enum::TS1 ; iTS <= tsReg_enum::TS5 ; iTS+=2 ) {
      auto its    = (tsReg_enum)iTS;

      for ( unsigned iRAD = tsRad_enum::IN ; iRAD <= tsRad_enum::OUT ; ++iRAD ) {
        auto irad = (tsRad_enum)iRAD;
        auto straightParam = straightSectionThermalParams.find( MapIndex(its,irad) );
	auto strsec = new StraightSection ( straightParam->second );
        ts._thermalShieldMap[its][irad] = std::unique_ptr<TSSection>( strsec );

      }
    }

    // Set thermal shield map - torus sections
    for ( unsigned iTS = tsReg_enum::TS2 ; iTS <= tsReg_enum::TS4 ; iTS+=2 ) {
      auto its = (tsReg_enum)iTS;

      for ( unsigned iRAD = tsRad_enum::IN ; iRAD <= tsRad_enum::OUT ; ++iRAD ) {
        auto irad = (tsRad_enum)iRAD;
        auto torusParam = torusSectionThermalParams.find( MapIndex(its,irad ) );

        ts._thermalShieldMap[its][irad] = std::unique_ptr<TSSection>( new TorusSection ( torusParam->second ) );

      }
    }

  }

  // The coils assemblies (CA) are approximated by a torus and tubes for now
  void BeamlineMaker::BuildTSCAs (const SimpleConfig& c, Beamline& bl) {

    TransportSolenoid& ts = bl._ts;

    const int  verbosityLevel = c.getInt("ts.cas.verbosityLevel", 0);

    typedef TransportSolenoid::TSCARegion::enum_type   tsCAReg_enum;

    for ( unsigned iTS = tsCAReg_enum::TS1 ; iTS <= tsCAReg_enum::TS5 ; ++iTS ) {

      auto its = static_cast<TransportSolenoid::TSCARegion>(iTS);

      std::vector<double> tmp_Radii;

      verbosityLevel && std::cout << __func__ << " filling radii info for "
                                  << its.name()
                                  << std::endl;

      std::string tsCAname = its.name();
      std::transform(tsCAname.begin(),tsCAname.end(),tsCAname.begin(),
                     [](char c){return std::tolower(c);});
      c.getVectorDouble( tsCAname+".cas.radii", tmp_Radii);
      ts._caRadiiMap[its] = tmp_Radii;

    }

    // getting info about stright sections
    typedef TransportSolenoid::TSRegion::enum_type     tsReg_enum;
    typedef TransportSolenoid::TSRadialPart::enum_type tsRad_enum;

    const StraightSection & ts1 = *(ts.getTSCryo<StraightSection>(tsReg_enum::TS1,tsRad_enum::IN));
    const StraightSection & ts3 = *(ts.getTSCryo<StraightSection>(tsReg_enum::TS3,tsRad_enum::IN));
    const StraightSection & ts5 = *(ts.getTSCryo<StraightSection>(tsReg_enum::TS5,tsRad_enum::IN));

    ts._caMaterial = c.getString("ts.cas.materialName");
    double ts3HalfLength = ts3.getHalfLength();
    double ts3udHalfLength = c.getDouble("ts3ud.cas.halfLength");
    double ts3udgapHalfLength = c.getDouble("ts3udgap.cas.halfLength");
    double ts3uuHalfLength = 0.5*(ts3udHalfLength - ts3udgapHalfLength);
    double ts3ddHalfLength = ts3uuHalfLength;

    // TS1,3,5 will be cones based on the coil radii and positions
    // TS3 will "splice" TS2 & 4

    verbosityLevel && std::cout << __func__ << " ts._caMaterial "  << ts._caMaterial << std::endl;

    // CA map parameters - torus sections
    for ( unsigned iTS = tsCAReg_enum::TS1 ; iTS <= tsCAReg_enum::TS5 ; ++iTS ) {

      auto its = static_cast<TransportSolenoid::TSCARegion>(iTS);

      verbosityLevel && std::cout << __func__ << " loop begin for: "
                                  << its.name()
                                  << std::endl;

      if ( its==tsCAReg_enum::TS1 ) {

        verbosityLevel && std::cout << __func__ << " making "
                                    << its.name()
                                    << std::endl;


        const Cone coneSectionParams( ts.caRadii(its)[0], ts.caRadii(its)[1],
                                      ts.caRadii(its)[2], ts.caRadii(its)[3],
                                      ts1.getHalfLength()-ts.endWallU1_halfLength() - ts1.getLengthDiff()/2.,
                                      0.0, CLHEP::twopi,
                                      ts1.getGlobal()+CLHEP::Hep3Vector(0.0,0.0,ts.endWallU1_halfLength()+ts1.getLengthDiff()/2.),
                                      CLHEP::HepRotation(),
                                      ts._caMaterial);

        verbosityLevel && std::cout << __func__ << " coneSectionParams.materialName() "
                                    << coneSectionParams.materialName()
                                    << std::endl;

        ts._caMap[its] = std::unique_ptr<TSSection>( new ConeSection ( coneSectionParams ) );

      }

      if ( its==tsCAReg_enum::TS5 ) {

        verbosityLevel && std::cout << __func__ << " making "
                                    << its.name()
                                    << std::endl;

        const Cone coneSectionParams( ts.caRadii(its)[0], ts.caRadii(its)[1],
                                      ts.caRadii(its)[2], ts.caRadii(its)[3],
                                      ts5.getHalfLength()-ts.endWallD_halfLength()-ts5.getLengthDiff()/2.,
                                      0.0, CLHEP::twopi,
                                      ts5.getGlobal()-CLHEP::Hep3Vector(0.0,0.0,ts.endWallD_halfLength()+ts5.getLengthDiff()/2.),
                                      CLHEP::HepRotation(),
                                      ts._caMaterial);

        verbosityLevel && std::cout << __func__ << " coneSectionParams.materialName() "
                                    << coneSectionParams.materialName()
                                    << std::endl;

       ts._caMap[its] = std::unique_ptr<TSSection>( new ConeSection ( coneSectionParams ) );

      }

      if ( its==tsCAReg_enum::TS2 || its==tsCAReg_enum::TS4 ) {

        verbosityLevel && std::cout << __func__ << " making "
                                    << its.name()
                                    << std::endl;

        CLHEP::Hep3Vector originInMu2e = (its==tsCAReg_enum::TS2) ?
          CLHEP::Hep3Vector(  ts3HalfLength, 0., -ts.torusRadius() ) :
          CLHEP::Hep3Vector( -ts3HalfLength, 0.,  ts.torusRadius() );

        double phi0 = (its==tsCAReg_enum::TS2) ? 1.5*CLHEP::pi : CLHEP::halfpi;

        const Torus torusSectionParams( ts.torusRadius(),
                                        ts.innerCARadius(its),
                                        ts.outerCARadius(its),
                                        phi0, CLHEP::halfpi,
                                        originInMu2e,
                                        CLHEP::HepRotation(CLHEP::HepRotationX(CLHEP::halfpi)),
                                        ts._caMaterial);

        verbosityLevel && std::cout << __func__
                                    << " torusSectionParams.materialName() "
                                    << torusSectionParams.materialName() << std::endl;

        ts._caMap[its] = std::unique_ptr<TSSection>( new TorusSection ( torusSectionParams ) );

      }

      if ( its==tsCAReg_enum::TS3uu  ) {

        verbosityLevel && std::cout << __func__ << " making "
                                    << its.name()
                                    << std::endl;

        const CLHEP::Hep3Vector ts3uupos( ts3uuHalfLength + ts3udgapHalfLength, 0., 0 );

        const Tube straightSectionParams (ts.innerCARadius(its),
                                          ts.outerCARadius(its),
                                          ts3uuHalfLength,
                                          ts3uupos,
                                          CLHEP::HepRotation(CLHEP::HepRotationY((CLHEP::halfpi))),
                                          0.0, CLHEP::twopi,
                                          ts._caMaterial);

        verbosityLevel && std::cout << __func__ << " straightSectionParams.materialName() "
                                    << straightSectionParams.materialName()
                                    << std::endl;

        ts._caMap[its] = std::unique_ptr<TSSection>( new StraightSection ( straightSectionParams ) );

      }


      if ( its==tsCAReg_enum::TS3dd  ) {

        verbosityLevel && std::cout << __func__ << " making "
                                    << its.name()
                                    << std::endl;

        const CLHEP::Hep3Vector ts3ddpos( -ts3ddHalfLength - ts3udgapHalfLength, 0., 0 );
        const Tube straightSectionParams (ts.innerCARadius(its),
                                          ts.outerCARadius(its),
                                          ts3ddHalfLength,
                                          ts3ddpos,
                                          CLHEP::HepRotation(CLHEP::HepRotationY((CLHEP::halfpi))),
                                          0.0, CLHEP::twopi,
                                          ts._caMaterial);

        verbosityLevel && std::cout << __func__ << " straightSectionParams.materialName() "
                                    << straightSectionParams.materialName()
                                    << std::endl;

        ts._caMap[its] = std::unique_ptr<TSSection>( new StraightSection ( straightSectionParams ) );

      }


      if ( its==tsCAReg_enum::TS3u ||  its==tsCAReg_enum::TS3d  ) {

        verbosityLevel && std::cout << __func__ << " making "
                                    << its.name()
                                    << std::endl;

        verbosityLevel && std::cout << __func__
                                    << " ts3HalfLength       " << ts3HalfLength
                                    << " bl.solenoidOffset() " << bl.solenoidOffset()
                                    << " ts.torusRadius()    " << ts.torusRadius()
                                    << std::endl;

        double halfLength = (ts3HalfLength-ts3udHalfLength)*0.5;

        CLHEP::Hep3Vector originInMu2e = (its==tsCAReg_enum::TS3u) ?
          CLHEP::Hep3Vector( ts3udHalfLength+halfLength,0.0,0.0) :
          CLHEP::Hep3Vector(-ts3udHalfLength-halfLength,0.0,0.0) ;

        const Cone coneSectionParams( ts.caRadii(its)[0], ts.caRadii(its)[1],
                                      ts.caRadii(its)[2], ts.caRadii(its)[3],
                                      halfLength,
                                      0.0, CLHEP::twopi,
                                      originInMu2e,
                                      CLHEP::HepRotation(CLHEP::HepRotationY(CLHEP::halfpi)),
                                      ts._caMaterial);

        verbosityLevel && std::cout << __func__ << " caConsec.materialName() "
                                    << coneSectionParams.materialName()
                                    << std::endl;

        ts._caMap[its] = std::unique_ptr<TSSection>( new ConeSection ( coneSectionParams ) );

      }

      verbosityLevel && std::cout << __func__ << " made "
                                  << its.name()
                                  << std::endl;

    }

  }
  void BeamlineMaker::BuildTSCoils (const SimpleConfig& c, Beamline* bl ) {

    TransportSolenoid * ts = &bl->_ts;

    ts->_coilMaterial = c.getString("ts.coils.material");

    // Loop over TS regions
    for ( unsigned iTS = TransportSolenoid::TSRegion::TS1 ;
          iTS <= TransportSolenoid::TSRegion::TS5 ; ++iTS )
      {
        auto its = (TransportSolenoid::TSRegion)iTS;

        std::vector<double> tmp_rIn, tmp_rOut, tmp_sLength, tmp_xPos, tmp_zPos, tmp_yRotAngle;

        std::ostringstream prefix;
        prefix << "ts" << iTS;
        c.getVectorDouble( prefix.str()+".coils.rIn"      , tmp_rIn       , ts->getNCoils(its) );
        c.getVectorDouble( prefix.str()+".coils.rOut"     , tmp_rOut      , ts->getNCoils(its) );
        c.getVectorDouble( prefix.str()+".coils.sLength"  , tmp_sLength   , ts->getNCoils(its) );
        c.getVectorDouble( prefix.str()+".coils.xPos"     , tmp_xPos      , ts->getNCoils(its) );
        c.getVectorDouble( prefix.str()+".coils.zPos"     , tmp_zPos      , ts->getNCoils(its) );
        c.getVectorDouble( prefix.str()+".coils.yRotAngle", tmp_yRotAngle , ts->getNCoils(its) );

        // Loop over coils per TS region
        ts->_coilMap[its].reserve( ts->getNCoils( its ) );

        for ( unsigned i(0) ; i < ts->getNCoils( its ) ; ++i ) {
          ts->_coilMap[its].emplace_back( tmp_xPos.at(i),  // position
                                          ts->getTSCryo(its,TransportSolenoid::TSRadialPart::IN)->getGlobal().y(),
                                          tmp_zPos.at(i),
                                          tmp_rIn.at(i),   // tube parameters
                                          tmp_rOut.at(i),
                                          tmp_sLength.at(i)*0.5,
                                          CLHEP::HepRotation(CLHEP::HepRotationY(-tmp_yRotAngle.at(i)*CLHEP::degree) ) );
        }

      }

  }

  void BeamlineMaker::BuildTSCollimators (const SimpleConfig& c, TransportSolenoid* ts) {

    // Set collimators
    double coll1HalfLength = c.getDouble("ts.coll1.halfLength");
    double coll3HalfLength = c.getDouble("ts.coll3.halfLength");
    double coll5HalfLength = c.getDouble("ts.coll5.halfLength");
    double coll3Hole       = c.getDouble("ts.coll3.hole");
    double collFlangeHalfLength = c.getDouble("ts.coll.Flange.halfLength");

    CollimatorTS1 & coll1  = ts->_coll1 ;
    CollimatorTS3 & coll31 = ts->_coll31;
    CollimatorTS3 & coll32 = ts->_coll32;
    CollimatorTS5 & coll51 = ts->_coll51;
    CollimatorTS5 & coll52 = ts->_coll52;
    CollimatorTS5 & coll53 = ts->_coll53;

    coll1 .set(coll1HalfLength,CLHEP::Hep3Vector(0.,0.,c.getDouble("ts.coll1.sOffset" ,0.)));
    coll31.set(coll3HalfLength,CLHEP::Hep3Vector(0.,0.,c.getDouble("ts.coll31.sOffset",0.)-coll3HalfLength-coll3Hole/2));
    coll32.set(coll3HalfLength,CLHEP::Hep3Vector(0.,0.,c.getDouble("ts.coll32.sOffset",0.)+coll3HalfLength+coll3Hole/2));
    coll51.set(coll5HalfLength,CLHEP::Hep3Vector(0.,0.,c.getDouble("ts.coll5.sOffset" ,0.)));
    double ts5HL = c.getDouble("ts.ts5.halfLength");
    double coll52HL = c.getDouble("ts.coll52.halfLength", ts5HL);
    coll52.set(coll52HL,CLHEP::Hep3Vector(0.,0.,ts5HL - coll52HL));
    coll53.set(collFlangeHalfLength,CLHEP::Hep3Vector(0.,0., c.getDouble("ts.ts5.halfLength") - collFlangeHalfLength - 2*c.getDouble("ts.tsDendWall.halfLength")));


    // TS1
    coll1._rIn1      = c.getDouble("ts.coll1.innerRadius1",0.);
    coll1._rIn2      = c.getDouble("ts.coll1.innerRadius2",0.);
    coll1._rIn3      = c.getDouble("ts.coll1.innerRadius3",0.);
    coll1._rIn4      = c.getDouble("ts.coll1.innerRadius4",0.);
    coll1._rOut4     = c.getDouble("ts.coll1.outerRadius4",0.);
    coll1._rOut1     = c.getDouble("ts.coll1.outerRadius1",-1.); //if < 0 is ignored
    coll1._material1 = c.getString("ts.coll1.material1Name");
    coll1._material2 = c.getString("ts.coll1.material2Name");
    coll1._material3 = c.getString("ts.coll1.material3Name","None");

    coll1._collarHalfLength = c.getDouble("pbar.coll1Out.halfLength",100.0);
    coll1._collarZ          = c.getDouble("pbar.coll1Out.z",-3104.5); // in Mu2e
    coll1._collarMarginZ    = c.getDouble("pbar.coll1Out.zDiff",0.5);
    coll1._collarrIn        = c.getDouble("pbar.coll1Out.rIn",        120.0);
    coll1._collarphiBegin   = c.getDouble("pbar.coll1Out.phiBegin",   210.0);
    coll1._collarphiDelta   = c.getDouble("pbar.coll1Out.phiDelta",   120.0);

    // TS3
    coll32._hole             = coll31._hole              = c.getDouble("ts.coll3.hole",0.);
    coll32._holeRadius       = coll31._holeRadius        = c.getDouble("ts.coll3.holeRadius",0.);
    coll32._holeHalfHeight   = coll31._holeHalfHeight    = c.getDouble("ts.coll3.holeHalfHeight",0.);
    coll32._holeDisplacement = coll31._holeDisplacement  = c.getDouble("ts.coll3.holeDisplacement",0.);
    coll32._rotationAngle    = coll31._rotationAngle     = c.getDouble("ts.coll3.rotationAngle",0.);
    coll32._material         = coll31._material          = c.getString("ts.coll3.materialName");
    coll32._rOut             = coll31._rOut              = c.getDouble("ts.coll3.outerRadius");

    // For studies of flash mitigation, building a "flashblock", also known
    // as "Rick's" (Coleman) "Tooth"
    coll31._useFlashBlock     = c.getBool("ts.useFlashBlockUp",false);
    coll31._flashBlockHeight  = c.getDouble("ts.flashBlockUp.Height",0.0);
    coll31._flashBlockWidth   = c.getDouble("ts.flashBlockUp.Width",0.0);
    coll31._flashBlockLength  = c.getDouble("ts.flashBlockUp.Length",0.0);
    coll31._flashBlockTO      = c.getDouble("ts.flashBlockUp.TransOffset",0.);
    coll31._flashBlockLO      = c.getDouble("ts.flashBlockUp.LongOffset",0.);
    coll31._flashBlockMaterial= c.getString("ts.flashBlockUp.Material","BronzeC608");

    coll32._useFlashBlock     = c.getBool("ts.useFlashBlockDn",false);
    coll32._flashBlockHeight  = c.getDouble("ts.flashBlockDn.Height",0.0);
    coll32._flashBlockWidth   = c.getDouble("ts.flashBlockDn.Width",0.0);
    coll32._flashBlockLength  = c.getDouble("ts.flashBlockDn.Length",0.0);
    coll32._flashBlockTO      = c.getDouble("ts.flashBlockDn.TransOffset",0.);
    coll32._flashBlockLO      = c.getDouble("ts.flashBlockDn.LongOffset",0.);
    coll32._flashBlockMaterial= c.getString("ts.flashBlockDn.Material","BronzeC608");



    // TS5
    // TS5
    coll51._rIn         = c.getDouble("ts.coll5.Radius1",0.);
    coll51._rOut        = c.getDouble("ts.coll5.Radius2",0.);

    coll52._rIn         = c.getDouble("ts.coll5.Radius3In",coll51._rOut); //default to outer diameter of first layer
    coll52._rOut        = c.getDouble("ts.coll5.Radius3",0.);

    coll53._rIn         = c.getDouble("ts.coll.Flange.Radius1",0.);
    coll53._rOut        = c.getDouble("ts.coll.Flange.Radius2",0.);
    coll53._version     = c.getInt("ts.coll.Flange.version", 1);

    coll51._material    = c.getString("ts.coll5.material1Name");
    coll52._material    = c.getString("ts.coll5.material2Name");
    coll53._material    = c.getString("ts.coll5.material2Name");

  }

  void BeamlineMaker::BuildTSVacua(const SimpleConfig& c, TransportSolenoid* ts ) {

    typedef TransportSolenoid::TSRegion::enum_type     tsReg_enum;
    typedef TransportSolenoid::TSRadialPart::enum_type tsRad_enum;

    const StraightSection * ts1 = ts->getTSCryo<StraightSection>(tsReg_enum::TS1,tsRad_enum::IN);
    const StraightSection * ts3 = ts->getTSCryo<StraightSection>(tsReg_enum::TS3,tsRad_enum::IN);
    const StraightSection * ts5 = ts->getTSCryo<StraightSection>(tsReg_enum::TS5,tsRad_enum::IN);

    // Figure things out for TS1, which is special
    CLHEP::Hep3Vector ts1VacPos = ts1->getGlobal();
    const double zLocCryoEdge   = ts1->getGlobal().z()-ts1->getHalfLength();
    const double zLocCollEdge   = ts->getColl1().getLocal().z()+ts1->getGlobal().z()-ts->getColl1().halfLength();
    const double zVacExtension  = std::abs( zLocCollEdge - zLocCryoEdge );

    const CLHEP::Hep3Vector vacOffset( 0., 0., -zVacExtension*0.5 );
    ts1VacPos = ts1VacPos + vacOffset;

    // Adjust offset of coll1 to be wrt TS1Vacuum (eventual parent
    // volume, not wrt TS1 cryostat)
    ts->_coll1.adjustOffset( -vacOffset );

    // Vacuum map parameters - straight sections
    const std::map<unsigned,Tube> straightSectionParams = {
      { tsReg_enum::TS1, Tube(0., ts->ts1InnerRadius(), ts1->getHalfLength()+zVacExtension*0.5, ts1VacPos ) },
      { tsReg_enum::TS3, Tube(0., ts->ts3InnerRadius(), ts3->getHalfLength(), ts3->getGlobal(), CLHEP::HepRotation(CLHEP::HepRotationY(90.0*CLHEP::degree)) ) },
      { tsReg_enum::TS5, Tube(0., ts->ts5InnerRadius(), ts5->getHalfLength(), ts5->getGlobal() ) }
     };

    // Vacuum map parameters - torus sections
    const std::map<unsigned,Torus> torusSectionParams = {
      { tsReg_enum::TS2, Torus( ts->torusRadius(), 0., ts->ts2InnerRadius(),
                                1.5*CLHEP::pi, CLHEP::halfpi,
                                CLHEP::Hep3Vector( ts3->getHalfLength(), 0.,-ts->torusRadius() ),
                                CLHEP::HepRotation(CLHEP::HepRotationX(90.0*CLHEP::degree))) },
      { tsReg_enum::TS4, Torus( ts->torusRadius(), 0., ts->ts4InnerRadius(),
                                CLHEP::halfpi, CLHEP::halfpi,
                                CLHEP::Hep3Vector( -ts3->getHalfLength(), 0.,ts->torusRadius() ),
                                CLHEP::HepRotation(CLHEP::HepRotationX(90.0*CLHEP::degree))) }
    };

    // Set vacuum map - straight sections
    for ( unsigned iTS = tsReg_enum::TS1 ; iTS <= tsReg_enum::TS5 ; iTS+=2 ) {
      auto its    = (tsReg_enum)iTS;
      auto straightParam = straightSectionParams.find( its );
      ts->_vacuumMap[its] = std::unique_ptr<TSSection>( new StraightSection ( straightParam->second ) );
    }

    // Set cryo map - torus sections
    for ( unsigned iTS = tsReg_enum::TS2 ; iTS <= tsReg_enum::TS4 ; iTS+=2 ) {
      auto its = (tsReg_enum)iTS;
      auto torusParam    = torusSectionParams.find( its );
      ts->_vacuumMap[its] = std::unique_ptr<TSSection>( new TorusSection ( torusParam->second ) );
    }

  }

    void BeamlineMaker::BuildTSPolyLining(const SimpleConfig& c, TransportSolenoid* ts ) {

    typedef TransportSolenoid::TSRegion::enum_type     tsReg_enum;
    typedef TransportSolenoid::TSRadialPart::enum_type tsRad_enum;

    TorusSection    const * ts2 = ts->getTSCryo<TorusSection>( tsReg_enum::TS2, tsRad_enum::IN );
    TorusSection    const * ts4 = ts->getTSCryo<TorusSection>( tsReg_enum::TS4, tsRad_enum::IN );

    // Vacuum map parameters - torus sections
    const std::map<unsigned,Torus> torusSectionParams = {
      { tsReg_enum::TS2, Torus( ts->torusRadius(),
                                c.getDouble("ts.polyliner.rIn",0.),
                                c.getDouble("ts.polyliner.rOut",ts->ts2InnerRadius()),
                                ts2->phiStart()+5*CLHEP::degree,
                                ts2->deltaPhi()-2*5*CLHEP::degree,
                                ts2->getGlobal() ) },
      //                                *ts2->getRotation() ) },
      { tsReg_enum::TS4, Torus( ts->torusRadius(),
                                c.getDouble("ts.polyliner.rIn",0.),
                                c.getDouble("ts.polyliner.rOut",ts->ts4InnerRadius()),
                                ts4->phiStart()+5*CLHEP::degree,
                                ts4->deltaPhi()-2*5*CLHEP::degree,
                                ts4->getGlobal() ) }
      //                                *ts4->getRotation() ) }
    };

    // Set cryo map - torus sections
    for ( unsigned iTS = tsReg_enum::TS2 ; iTS <= tsReg_enum::TS4 ; iTS+=2 ) {
      auto its = (tsReg_enum)iTS;
      auto torusParam    = torusSectionParams.find( its );
      ts->_polyLiningMap[its] = std::unique_ptr<TorusSection>( new TorusSection ( torusParam->second ) );
      ts->_polyLiningMap[its]->setMaterial( c.getString("ts.polyliner.materialName") );
    }

  }

  void BeamlineMaker::BuildPbarWindow(const SimpleConfig& c, TransportSolenoid* ts){

    PbarWindow & pbarWindow ( ts->_pbarWindow );
    pbarWindow._version  = c.getInt("pbar.version",1);
    pbarWindow._shape    = c.getString("pbar.Type","disk");
    pbarWindow._material = c.getString("pbar.materialName");
    pbarWindow.set( c.getDouble("pbar.halfLength"),
                    CLHEP::Hep3Vector() );
    pbarWindow._rOut     = ts->ts3InnerRadius();

    // Parameters for wedge
    pbarWindow._y0       = c.getDouble("pbarwedge.y0",0.);
    pbarWindow._y1       = c.getDouble("pbarwedge.y1",0.);
    pbarWindow._dz0      = c.getDouble("pbarwedge.dz0",0.);
    pbarWindow._dz1      = c.getDouble("pbarwedge.dz1",0.);
    pbarWindow._wedgeZOffset = c.getDouble("pbarwedge.zOffset",0.);

    // for version 3 (made of strips)
    pbarWindow._diskRadius = c.getDouble("pbar.diskradius", 250.0);
    pbarWindow._nStrips    = c.getInt("pbarwedge.nStrips", 0);
    pbarWindow._width      = c.getDouble("pbarwedge.width",0.0);
    pbarWindow._stripThickness=c.getDouble("pbarwedge.stripThickness",0.0);
    if ( pbarWindow._nStrips > 0 ) {
      c.getVectorDouble("pbarwedge.stripHeights",pbarWindow._stripHeights,pbarWindow._nStrips);
    }
    // for version 4 (variable thickness strips)
    if (pbarWindow._nStrips > 0 && pbarWindow._version == 4) {
      c.getVectorDouble("pbarwedge.stripThicknesses",pbarWindow._stripThicknesses,pbarWindow._nStrips);
    }
    // allow for the wedge and window to have different materials
    pbarWindow._wedgeMaterial = c.getString("pbarwedge.wedgeMaterial",pbarWindow._material); //default to window material
  }

} // namespace mu2e

#endif
