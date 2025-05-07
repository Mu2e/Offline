// Andrei Gaponenko, 2011

#include "Offline/GeometryService/inc/ProtonBeamDumpMaker.hh"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <iostream>
#include <limits>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

#include "cetlib_except/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

namespace mu2e {

  namespace {

    CLHEP::Hep2Vector toOutlinePlane(CLHEP::Hep3Vector mu2ePoint) {
      return CLHEP::Hep2Vector(mu2ePoint.z(), mu2ePoint.x());
    }

    CLHEP::Hep3Vector toMu2eSpace(CLHEP::Hep2Vector outlinePoint, double y) {
      return CLHEP::Hep3Vector(outlinePoint.y(), y, outlinePoint.x());
    }

    struct Line {
      CLHEP::Hep2Vector A;
      CLHEP::Hep2Vector B;
      Line(CLHEP::Hep2Vector Ain, CLHEP::Hep2Vector Bin) : A{Ain}, B{Bin} {}
    };

    CLHEP::Hep2Vector intersection(Line L1, Line L2) {
      const auto AB = L1.A - L1.B;
      const auto DC = L2.B - L2.A;
      const auto AC = L1.A - L2.A;

      const double det = DC.x()*AB.y() - DC.y()*AB.x();
      const double tRHS = DC.x() * AC.y() - DC.y()*AC.x();

      if(std::abs(tRHS) < std::abs(det) * std::numeric_limits<double>::max()) { // can divide
        const double t = tRHS / det;
        return L1.A - t * AB;
      }
      else {
        throw cet::exception("GEOM")<<"ProtonBeamDumpMaker: internal error, det is too small\n";
      }
    }

  } // end of anonymous namespace

  //================================================================
  std::unique_ptr<ProtonBeamDump> ProtonBeamDumpMaker::make(const SimpleConfig& c,
                                                            const Mu2eHall& hall)
  {
    using CLHEP::Hep3Vector;
    using CLHEP::Hep2Vector;

    std::unique_ptr<ProtonBeamDump> dump(new ProtonBeamDump());

    int verbose = c.getInt("protonBeamDump.verbosityLevel", 0);

    dump->_coreCenterInMu2e = c.getHep3Vector("protonBeamDump.coreCenterInMu2e");
    const double coreRotY = dump->_coreRotY = c.getDouble("protonBeamDump.coreRotY") * CLHEP::degree;
    dump->_coreRotationInMu2e.rotateY(coreRotY);
    //
    // Dump core position and rotation have been initialized, from
    // this point on mu2eToBeamDump_position() and other coordinate
    // conversion methods can be used.

    c.getVectorDouble("protonBeamDump.coreHalfSize", dump->_coreHalfSize, 3);
    c.getVectorDouble("protonBeamDump.neutronCaveHalfSize", dump->_neutronCaveHalfSize, 3);
    c.getVectorDouble("protonBeamDump.mouthHalfSize", dump->_mouthHalfSize, 3);

    const double coreAirSideGap = c.getDouble("protonBeamDump.coreAirSideGap");
    const double coreAirTopGap = c.getDouble("protonBeamDump.coreAirTopGap");
    const double coreAirBottomGap = c.getDouble("protonBeamDump.coreAirBottomGap");

    dump->_coreAirHalfSize.resize(3);
    dump->_coreAirHalfSize[0] = dump->_coreHalfSize[0] + coreAirSideGap;
    dump->_coreAirHalfSize[1] = dump->_coreHalfSize[1] + (coreAirTopGap+coreAirBottomGap)/2;
    dump->_coreAirHalfSize[2] = dump->_coreHalfSize[2];

    dump->_neutronCaveCenterInMu2e =
      dump->beamDumpToMu2e_position(Hep3Vector(0,
                                               0,
                                               dump->_coreHalfSize[2] + dump->_neutronCaveHalfSize[2]));

    dump->_mouthCenterInMu2e =
      dump->beamDumpToMu2e_position(Hep3Vector(0,
                                               0,
                                               dump->_coreHalfSize[2]
                                               + 2*dump->_neutronCaveHalfSize[2]
                                               + dump->_mouthHalfSize[2]
                                               ));

    dump->_coreAirCenterInMu2e =
      dump->beamDumpToMu2e_position(Hep3Vector(0,
                                               (coreAirTopGap-coreAirBottomGap)/2,
                                               0
                                               ));


    //----------------------------------------------------------------
    // Get relevant Hall solids

    const ExtrudedSolid& psArea  = hall.getBldgSolid("psArea");         // bottom of dump
    const ExtrudedSolid& psWall  = hall.getBldgSolid("psWallUpper");    // contains back corners of dump slab
    const ExtrudedSolid& psCeil  = hall.getBldgSolid("psAreaCeilingSW");// top of dump and ExtMon cutout

    // The scale of an artificial gap, in mm, between touching volumes
    // needed to eliminate G4-reported volume overlaps that are
    // artifacts of rounding errors.
    const double ovlgap = 0.5;

    // unit vectors in the beam dump reference frame
    const Hep2Vector nx = toOutlinePlane(dump->beamDumpToMu2e_momentum(Hep3Vector(1,0,0)));
    const Hep2Vector nz = toOutlinePlane(dump->beamDumpToMu2e_momentum(Hep3Vector(0,0,1)));

    //----------------------------------------------------------------
    // For the dump shielding we define an extruded solid that extends
    // the side faces of the nominal dump horizontally to the enclosure
    // walls and vertically to the ceiling and floor. In this way we
    // accommodate the slight diference between the enclosure
    // orientation and the nominal dump rotation.

    // The vertical extent of beam dump concrete is from psArea concrete to the ceiling.
    const double concreteYmin = psArea.getOffsetFromMu2eOrigin().y() + psArea.getYhalfThickness();
    const double concreteYmax = psCeil.getOffsetFromMu2eOrigin().y() - psCeil.getYhalfThickness();

    dump->_dumpConcreteHalfHeight = (concreteYmax - concreteYmin)/2;
    dump->_dumpConcreteCenterInMu2e = Hep3Vector(0, (concreteYmax+concreteYmin)/2, 0);

    // Beam dump outline is contained by psWallUpper vertices 13-14-15-16,
    // with the back of the dump following the 14-15 line.
    // The front corners are at the intersections of straight lines
    // through psWallUpper vertices 13-14 or 15-16 and the beam dump face.

    const double dumpCoreCenterToConcreteFaceDistance =
      dump->_coreHalfSize[2]
      + 2 * dump->_neutronCaveHalfSize[2]
      + 2*dump->_mouthHalfSize[2];

    // Beam dump face line in the plan view
    Line faceLine(
                  // Concrete face point on dump axis
                  toOutlinePlane(dump->beamDumpToMu2e_position(Hep3Vector(0, 0, dumpCoreCenterToConcreteFaceDistance))),
                  // Move 10m sideway.  An arbitrary point to define the line.
                  toOutlinePlane(dump->beamDumpToMu2e_position(Hep3Vector(10000., 0, dumpCoreCenterToConcreteFaceDistance)))
                  );

    const auto& pswVertices =  psWall.getVertices();
    const Hep2Vector& psWoff2 = toOutlinePlane(psWall.getOffsetFromMu2eOrigin());

    dump->_dumpConcreteOutline.emplace_back(pswVertices[14]+psWoff2 + ovlgap*(-nx+nz));
    dump->_dumpConcreteOutline.emplace_back(pswVertices[15]+psWoff2 + ovlgap*( nx+nz));
    dump->_dumpConcreteOutline.emplace_back(intersection(faceLine, Line(pswVertices[15] + psWoff2, pswVertices[16] + psWoff2)) + ovlgap*( nx));
    dump->_dumpConcreteOutline.emplace_back(intersection(faceLine, Line(pswVertices[13] + psWoff2, pswVertices[14] + psWoff2)) + ovlgap*(-nx));

    // The ExtMon room subtraction.
    // The the back is the same as the concrete back.  Add a margin for the G4 boolean.
    const double booleanVolumeMargin = 50; // mm
    dump->_extMonSubtractionOutline.emplace_back(pswVertices[14]+psWoff2 + booleanVolumeMargin*( nx-nz));
    dump->_extMonSubtractionOutline.emplace_back(pswVertices[15]+psWoff2 + booleanVolumeMargin*(-nx-nz));

    // The front of the subtraction is lined up with psCeil vertices 6
    // and 7.  This line cuts through the bulk of the concrete and the
    // position of the line should be exact (no margin), but the ends
    // should protrude sideways. The protrusion is taken care of by
    // using the vertices that are spaced wider than the beam dump
    // concrete width.

    const auto& pscVertices =  psCeil.getVertices();
    const Hep2Vector& psCoff2 = toOutlinePlane(psCeil.getOffsetFromMu2eOrigin());
    dump->_extMonSubtractionOutline.emplace_back(pscVertices[6]+psCoff2);
    dump->_extMonSubtractionOutline.emplace_back(pscVertices[7]+psCoff2);

    // The bottom of the cutout should be exact.  It lines up with psWallUpper
    const double extMonCutoutYmin = psWall.getOffsetFromMu2eOrigin().y() + psWall.getYhalfThickness();
    // The top should go  above the top of the beam dump concrete
    const double extMonCutoutYmax = concreteYmax + booleanVolumeMargin;

    dump->_extMonSubtractionHalfHeight = (extMonCutoutYmax - extMonCutoutYmin)/2;
    dump->_extMonSubtractionCenterInMu2e = Hep3Vector(0, (extMonCutoutYmax + extMonCutoutYmin)/2, 0);


    //----------------------------------------------------------------
    // The extra steel on top of the core

    // Top steel clearances are from the side walls, which coincide
    // with the sides of the dump concrete outline.  The sides are
    // nominally parallel to the beam dump axis, so it should not
    // matter what outline vertex is used.  We look at the both
    // vertices and select the more constraining number to account for
    // imperfections.

    const double beamLeftWallX /*in the beam dump coordinates*/
      = std::max(
                 dump->mu2eToBeamDump_position(toMu2eSpace(dump->_dumpConcreteOutline[1], 0)).x(),
                 dump->mu2eToBeamDump_position(toMu2eSpace(dump->_dumpConcreteOutline[2], 0)).x()
                 );

    const double beamRightWallX /*in the beam dump coordinates*/
      = std::min(
                 dump->mu2eToBeamDump_position(toMu2eSpace(dump->_dumpConcreteOutline[0], 0)).x(),
                 dump->mu2eToBeamDump_position(toMu2eSpace(dump->_dumpConcreteOutline[3], 0)).x()
                 );

    const double beamDumpBackZ /*in the beam dump coordinates*/
      = std::max(
                 dump->mu2eToBeamDump_position(toMu2eSpace(dump->_dumpConcreteOutline[0], 0)).z(),
                 dump->mu2eToBeamDump_position(toMu2eSpace(dump->_dumpConcreteOutline[1], 0)).z()
                 );

    const double topSteelFrontZ = dump->_coreHalfSize[2]; // reaches out to the front of the core

    //----------------
    std::vector<double> topSteelFlatWallClearance;
    c.getVectorDouble("protonBeamDump.topSteelFlat.wallClearance", topSteelFlatWallClearance);

    const double topSteelFlatXmin = beamLeftWallX + topSteelFlatWallClearance[0];
    const double topSteelFlatXmax = beamRightWallX - topSteelFlatWallClearance[1];

    dump->_topSteelFlatHalfSize.resize(3);
    dump->_topSteelFlatHalfSize[0] = (topSteelFlatXmax - topSteelFlatXmin)/2;
    dump->_topSteelFlatHalfSize[1] = c.getDouble("protonBeamDump.topSteelFlat.thickness")/2;
    dump->_topSteelFlatHalfSize[2] = (topSteelFrontZ - beamDumpBackZ)/2;

    dump->_topSteelFlatCenterInMu2e =
      dump->beamDumpToMu2e_position(Hep3Vector(
                                               (topSteelFlatXmax + topSteelFlatXmin)/2,
                                               dump->_coreHalfSize[1] + coreAirTopGap + dump->_topSteelFlatHalfSize[1],
                                               dump->_coreHalfSize[2] - dump->_topSteelFlatHalfSize[2]
                                               ));

    //----------------
    std::vector<double> topSteelScallopedWallClearance;
    c.getVectorDouble("protonBeamDump.topSteelScalloped.wallClearance", topSteelScallopedWallClearance);

    const double topSteelScallopedXmin = beamLeftWallX + topSteelScallopedWallClearance[0];
    const double topSteelScallopedXmax = beamRightWallX - topSteelScallopedWallClearance[1];

    dump->_topSteelScallopedHalfSize.resize(3);
    dump->_topSteelScallopedHalfSize[0] = (topSteelScallopedXmax - topSteelScallopedXmin)/2;
    dump->_topSteelScallopedHalfSize[1] = c.getDouble("protonBeamDump.topSteelScalloped.thickness")/2;
    dump->_topSteelScallopedHalfSize[2] = (topSteelFrontZ - beamDumpBackZ)/2;

    dump->_topSteelScallopedCenterInMu2e =
      dump->beamDumpToMu2e_position(Hep3Vector(
                                               (topSteelScallopedXmax + topSteelScallopedXmin)/2,

                                               dump->_coreHalfSize[1]
                                               + coreAirTopGap
                                               + 2*dump->_topSteelFlatHalfSize[1]
                                               + dump->_topSteelScallopedHalfSize[1],

                                               dump->_coreHalfSize[2] - dump->_topSteelScallopedHalfSize[2]
                                               ));

    dump->_scallopDistanceToCollimator = c.getDouble("protonBeamDump.topSteelScalloped.distanceToCollimator");

    //----------------------------------------------------------------
    if(verbose) {
      std::cout<<"ProtonBeamDumpMaker"<<": ProtonBeamDump core center in mu2e = "<<dump->_coreCenterInMu2e<<std::endl;

      std::cout<<"core halfsize:  "
               << dump->_coreHalfSize[0]<<", "
               << dump->_coreHalfSize[1]<<", "
               << dump->_coreHalfSize[2]<<", "
               << std::endl;

      for( std::size_t i=0; i<dump->_dumpConcreteOutline.size(); i++ ) {
        // whitespace or tab separated print for the ease of use with gnuplot
        std::cout<<"dumpConcreteOutline  "
                 <<dump->_dumpConcreteOutline[i][0]<<"\t"
                 <<dump->_dumpConcreteOutline[i][1]<<"\n";
      }
      std::cout<<"dumpConcreteHalfHeight = "<<dump->_dumpConcreteHalfHeight<<std::endl;
      std::cout<<"dumpConcreteCenterInMu2e = "<<dump->_dumpConcreteCenterInMu2e<<std::endl;

      for( std::size_t i=0; i<dump->_extMonSubtractionOutline.size(); i++ ) {
        // whitespace or tab separated print for the ease of use with gnuplot
        std::cout<<"extMonSubtractionOutline  "
                 <<dump->_extMonSubtractionOutline[i][0]<<"\t"
                 <<dump->_extMonSubtractionOutline[i][1]<<"\n";
      }
      std::cout<<"extMonSubtractionHalfHeight = "<<dump->_extMonSubtractionHalfHeight<<std::endl;
      std::cout<<"extMonSubtractionCenterInMu2e = "<<dump->_extMonSubtractionCenterInMu2e<<std::endl;

      std::cout<<"dump mouth halfsize:  "
               << dump->_mouthHalfSize[0]<<", "
               << dump->_mouthHalfSize[1]<<", "
               << dump->_mouthHalfSize[2]<<", "
               << std::endl;

      std::cout<<"dump neutronCave halfsize:  "
               << dump->_neutronCaveHalfSize[0]<<", "
               << dump->_neutronCaveHalfSize[1]<<", "
               << dump->_neutronCaveHalfSize[2]<<", "
               << std::endl;

      std::cout<<"Flat steel halfsize:  "
               << dump->_topSteelFlatHalfSize[0]<<", "
               << dump->_topSteelFlatHalfSize[1]<<", "
               << dump->_topSteelFlatHalfSize[2]<<", "
               << std::endl;

      std::cout<<"Flat steel Center in Mu2e:  " << dump->_topSteelFlatCenterInMu2e << std::endl;

      std::cout<<"Scalloped steel halfsize:  "
               << dump->_topSteelScallopedHalfSize[0]<<", "
               << dump->_topSteelScallopedHalfSize[1]<<", "
               << dump->_topSteelScallopedHalfSize[2]<<", "
               << std::endl;

      std::cout<<"Scalloped steel Center in Mu2e:  " << dump->_topSteelScallopedCenterInMu2e << std::endl;

      std::cout<<"Scalloped steel distance to collimator:  " << dump->_scallopDistanceToCollimator << std::endl;
    }

    return dump;

  } // make()

} // namespace mu2e
