// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuildingMaker.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
//#define AGDEBUG(stuff)

namespace mu2e {

  //================================================================
  ExtMonFNALBuilding::CollimatorExtMonFNAL
  ExtMonFNALBuildingMaker::readCollimatorExtMonFNAL(const std::string& name,
                                                    double zLength,
                                                    double angleH,
                                                    double angleV,
                                                    const SimpleConfig& c) {

    ExtMonFNALBuilding::CollimatorExtMonFNAL col;

    col._name = name;
    col._horizontalLength = zLength;
    c.getVectorDouble("extMonFNAL."+name+".channelWidth",            col._channelWidth           , 2);
    c.getVectorDouble("extMonFNAL."+name+".channelHeigh",            col._channelHeight          , 2);
    c.getVectorDouble("extMonFNAL."+name+".alignmentPlugRadius",     col._alignmentPlugRadius    , 2);
    c.getVectorDouble("extMonFNAL."+name+".alignmentHoleRClearance", col._alignmentHoleRClearance, 2);
    col._radiusTransitiondZ = c.getDouble("extMonFNAL."+name+".radiusTransitiondZ");
    col._angleH = angleH;
    col._angleV = angleV;

    return col;
  }

  //================================================================
  ExtMonFNALBuilding::FilterMagnetExtMonFNAL
  ExtMonFNALBuildingMaker::readFilterMagnetExtMonFNAL(const SimpleConfig& c) {
    ExtMonFNALBuilding::FilterMagnetExtMonFNAL mag;
    c.getVectorDouble("extMonFNAL.magnet.outerHalfSize", mag._outerHalfSize, 3);
    mag._apertureWidth = c.getDouble("extMonFNAL.magnet.apertureWidth") * CLHEP::mm;
    mag._apertureHeight = c.getDouble("extMonFNAL.magnet.apertureHeight") * CLHEP::mm;
    mag._fieldStrength = c.getDouble("extMonFNAL.magnet.fieldStrength") * CLHEP::tesla;
    return mag;
  }

  //================================================================
  std::auto_ptr<ExtMonFNALBuilding> ExtMonFNALBuildingMaker::make(const SimpleConfig& c, const ProtonBeamDump& dump) {
    using CLHEP::Hep3Vector;

    std::auto_ptr<ExtMonFNALBuilding> emfb(new ExtMonFNALBuilding());

    int verbose = c.getInt("extMonFNAL.verbosityLevel");

    emfb->roomInsideFullHeight_ = c.getDouble("extMonFNAL.room.insideFullHeight");
    emfb->roomInsideYmin_ = dump.coreCenterInMu2e()[1] + dump.backShieldingHalfSize()[1];
    emfb->roomInsideYmax_ = emfb->roomInsideYmin_ + emfb->roomInsideFullHeight_;

    const double wallThickness = emfb->roomWallThickness_ = c.getDouble("extMonFNAL.room.wall.thickness");
    emfb->roomFloorThickness_ = c.getDouble("extMonFNAL.room.floor.thickness");
    emfb->roomCeilingThickness_ = c.getDouble("extMonFNAL.room.ceiling.thickness");

    emfb->magnetRoomLength_ = c.getDouble("extMonFNAL.magnetRoomLength");
    emfb->coll2ShieldingDumpXmin_ = c.getDouble("extMonFNAL.collimator2.shielding.dumpXmin");
    emfb->coll2ShieldingDumpXmax_ = c.getDouble("extMonFNAL.collimator2.shielding.dumpXmax");

    //----------------------------------------------------------------
    // We need to rotate the extruded solid coordinates by +90 degrees to bring it to the horizontal plane.
    // No other rotation is needed since the outline is already in Mu2e (x,z), i.e. rorated with the dump
    emfb->roomRefPointInMu2e_ = Hep3Vector(0, (emfb->roomInsideYmax() + emfb->roomInsideYmin())/2, 0);
    emfb->roomRotationInMu2e_ = CLHEP::HepRotationX(+90.*CLHEP::degree);

    //----------------------------------------------------------------
    emfb->_filterEntranceOffsetX = c.getDouble("extMonFNAL.entranceOffsetX") * CLHEP::mm;
    emfb->_filterEntranceOffsetY = c.getDouble("extMonFNAL.entranceOffsetY") * CLHEP::mm;

    const double angleH = emfb->_filterAngleH =  c.getDouble("extMonFNAL.angleH") * CLHEP::radian;
    const double entranceAngleV = emfb->_filterEntranceAngleV =  c.getDouble("extMonFNAL.entranceAngleV") * CLHEP::radian;

    const double pNominal = emfb->_extMonFNAL_nominalMomentum = c.getDouble("extMonFNAL.nominalMomentum") * CLHEP::MeV;

    emfb->_filterMagnet = readFilterMagnetExtMonFNAL(c);

    const double col1zLength = 2*dump.frontShieldingHalfSize()[2];
    const double col2zLength = c.getDouble("extMonFNAL.collimator2.shielding.thickness");

    emfb->_collimator1 = readCollimatorExtMonFNAL("collimator1", col1zLength, angleH, entranceAngleV, c);
    emfb->_collimator2 = readCollimatorExtMonFNAL("collimator2",
                                                  col2zLength,
                                                  angleH,
                                                  entranceAngleV - 2 * emfb->_filterMagnet.trackBendHalfAngle(pNominal),
                                                  c);

    //----------------------------------------------------------------
    // Compute the placement of filter elements

    // collimator1
    const Hep3Vector collimator1CenterInDump(emfb->filterEntranceOffsetX()
                                             + 0.5*col1zLength*tan(emfb->collimator1().angleH()),

                                             emfb->filterEntranceOffsetY()
                                             + 0.5*col1zLength*tan(emfb->collimator1().angleV())/cos(emfb->filterAngleH()),

                                             dump.coreCenterDistanceToShieldingFace() - dump.frontShieldingHalfSize()[2]
                                             );

    emfb->_collimator1CenterInMu2e = dump.beamDumpToMu2e_position(collimator1CenterInDump);

    {
      CLHEP::HepRotation tmp(CLHEP::HepRotation::IDENTITY);
      tmp.rotateX(entranceAngleV).rotateY(-angleH);
      emfb->_collimator1RotationInMu2e = dump.coreRotationInMu2e() * tmp;
    }

    //----------------------------------------------------------------
    // Position the magnet.  Keep B field horizontal.  Then there is
    // no bend in the XZ plane and the magnet center should be on the
    // continuation of the collimator1 center line.
    //
    // The magnet should also have the same rotation in the projection
    // on XZ as the collimator(s).
    //
    // The reference trajectory: a positive particle with the nominal
    // momentum, travelling parallel to collimator1 axis at the bottom
    // center of collimator1 aperture, should enter magnet at the
    // bottom center of the magnet aperture, bend down, and exit a the
    // bottom center of magnet aperture.

    const double magnetAngleV =
      emfb->filterEntranceAngleV() - emfb->filterMagnet().trackBendHalfAngle(emfb->extMonFNAL_nominalMomentum());

    {
      CLHEP::HepRotation tmp(CLHEP::HepRotation::IDENTITY);
      tmp.rotateX(magnetAngleV).rotateY(-angleH);
      emfb->_filterMagnetRotationInMu2e = dump.coreRotationInMu2e() * tmp;
    }

    // The z position of the (center,bottom,center) point of the magnet aperture
    // in dump coordinates
    const double magnetRefZ = dump.coreCenterDistanceToShieldingFace()
      - 2*dump.frontShieldingHalfSize()[2]
      - c.getDouble("extMonFNAL.magnet.refDistanceToUpstreamWall");

    // Z of magnet entrance at the (center,bottom) of aperture (in dump coords)
    const double magnetEntranceZ = magnetRefZ
      + emfb->filterMagnet().outerHalfSize()[2] * cos(magnetAngleV) * cos(emfb->filterAngleH());

    // Height of the entrance point, dump coords
    const double magnetEntranceY = collimator1CenterInDump[1]
      + (collimator1CenterInDump[2] - magnetEntranceZ)*tan(emfb->filterEntranceAngleV())/cos(emfb->filterAngleH())

      // the entrance point should be on the trajectory of reference
      // particle at the *bottom* of the collimator1 channel - shift for that:
      - 0.5 * emfb->collimator1().channelHeight()[emfb->collimator1().channelHeight().size()-1]/cos(emfb->collimator1().angleV())
      ;

    AGDEBUG("magnetEntranceY = "<<magnetEntranceY<<" = "<<emfb->_collimator1CenterInMu2e[1]<<" - "<<dump.coreCenterInMu2e()[1]
            <<" + "<<(emfb->_collimator1CenterInMu2e[2] - (magnetEntranceZ+dump.coreCenterInMu2e()[2]))*tan(emfb->filterEntranceAngleV())/cos(emfb->filterAngleH())
            <<" - "<<0.5 * emfb->collimator1().channelHeight()[emfb->collimator1().channelHeight().size()-1]/cos(emfb->collimator1().angleV())
            );

    // X of the entrance point, dump coords
    const double magnetEntranceX = emfb->_filterEntranceOffsetX
      + tan(emfb->filterAngleH()) * (dump.coreCenterDistanceToShieldingFace() - magnetEntranceZ);

    // length from magnet entrance to magnet center in projection on (XZ)
    const double projLength =
      (emfb->filterMagnet().outerHalfSize()[2]* cos(magnetAngleV) - 0.5*emfb->filterMagnet().apertureHeight()*sin(magnetAngleV));

    // Compute the center of the magnet
    const Hep3Vector magnetCenterInDump
      (magnetEntranceX + sin(emfb->filterAngleH()) * projLength
       ,
       magnetEntranceY + sin(magnetAngleV)*emfb->filterMagnet().outerHalfSize()[2] + cos(magnetAngleV)*0.5*emfb->filterMagnet().apertureHeight()
       ,
       magnetEntranceZ - cos(emfb->filterAngleH()) * projLength
       );

    AGDEBUG("magnetCenterInDump = "<<magnetCenterInDump);

    emfb->_filterMagnetCenterInMu2e = dump.beamDumpToMu2e_position(magnetCenterInDump);

    //----------------------------------------------------------------
    // collimator2: position so that the reference trajectory enters it at the bottom center
    // of the aperture

    const double magnetRoomLength = c.getDouble("extMonFNAL.magnetRoomLength");

    // dump coords
    const double col2CenterZ = dump.coreCenterDistanceToShieldingFace()
      - 2*dump.frontShieldingHalfSize()[2]
      - magnetRoomLength
      - 0.5*emfb->_collimator2.horizontalLength()
      ;

    // dump coords
    const double col2CenterX = emfb->_filterEntranceOffsetX +
      tan(emfb->filterAngleH()) * (dump.coreCenterDistanceToShieldingFace() - col2CenterZ);

    // Z of the exit point from the magnet, at the (center,bottom) of aperture (in dump coords)
    const double magnetExitZ = magnetEntranceZ
      - 2*emfb->filterMagnet().outerHalfSize()[2] * cos(magnetAngleV) * cos(emfb->filterAngleH());

    // Height of the exit point (in dump coords)
    const double magnetExitY = magnetEntranceY + sin(magnetAngleV) * 2*emfb->_filterMagnet._outerHalfSize[2];



    const double col2CenterY = magnetExitY + (magnetExitZ - col2CenterZ)*tan(emfb->_collimator2.angleV())/cos(emfb->filterAngleH())
      + 0.5 * emfb->_collimator2._channelHeight[0]/cos(emfb->_collimator2.angleV());

    emfb->_collimator2CenterInMu2e = dump.beamDumpToMu2e_position(Hep3Vector(col2CenterX, col2CenterY, col2CenterZ));

    {
      CLHEP::HepRotation tmp(CLHEP::HepRotation::IDENTITY);
      tmp.rotateX(emfb->_collimator2._angleV).rotateY(-angleH);
      emfb->_collimator2RotationInMu2e = dump.coreRotationInMu2e() * tmp;
    }

    //----------------------------------------------------------------
    if(verbose) {
      std::cout<<"ExtMonFNALBuildingMaker"<<": Filter entrance offsets  = ("<<emfb->_filterEntranceOffsetX<<", "<<emfb->_filterEntranceOffsetY<<")"<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filter nominal momentum = "<<emfb->extMonFNAL_nominalMomentum()/CLHEP::GeV<<" GeV/c"<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filterAngleH = "<<emfb->filterAngleH()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filterAngleH in Mu2e, degrees= "<<(dump.coreRotY() - emfb->filterAngleH())/CLHEP::degree<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filter half bend angle  = "<<emfb->filterMagnet().trackBendHalfAngle(emfb->extMonFNAL_nominalMomentum())<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filter.angleV = "<<emfb->filterEntranceAngleV()
               <<", c1.angleV  = "<<emfb->collimator1().angleV()
               <<", magnet.angleV = "<<magnetAngleV
               <<", c2.angleV() = "<<emfb->collimator2().angleV()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator1CenterInMu2e = "<<emfb->_collimator1CenterInMu2e<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator1.horizontalLength = "<<emfb->_collimator1._horizontalLength<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": filterMagnetCenterInMu2e = "<<emfb->_filterMagnetCenterInMu2e<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator2CenterInMu2e = "<<emfb->_collimator2CenterInMu2e<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator2.horizontalLength = "<<emfb->_collimator2._horizontalLength<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": ExtMonFNALBuilding::filterEntranceInMu2e() = "<<emfb->filterEntranceInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": ExtMonFNALBuilding::filterExitInMu2e() = "<<emfb->filterExitInMu2e()<<std::endl;
    }

    //----------------------------------------------------------------
    std::vector<double> roomInsideX, roomInsideZ;
    c.getVectorDouble("extMonFNAL.roomInsideX", roomInsideX);
    c.getVectorDouble("extMonFNAL.roomInsideZ", roomInsideZ);
    if(roomInsideZ.size() != roomInsideX.size()) {
      throw cet::exception("GEOM")<<"ExtMonFNALBuildingMaker: ERROR: "
                                  <<"different sizes for extMonFNAL.roomInsideX (="<<roomInsideX.size()
                                  <<") and extMonFNAL.roomInsideZ (="<<roomInsideZ.size()<<")"
        ;
    }

    if(roomInsideX.size() < 3) {
      throw cet::exception("GEOM")<<"ExtMonFNALBuildingMaker: ERROR: "
                                  <<"too few points in extMonFNAL.roomInsideX, roomInsideZ: "
                                  <<" need at least 3 points to define the outline, got "
                                  <<roomInsideX.size()
        ;
    }

    // Convert inside outline to Mu2e coords
    for(unsigned i=0; i<roomInsideX.size(); ++i) {
      // Convert inputs to Mu2e coordinates
      Hep3Vector pointInDump(roomInsideX[i],
                             0,
                             roomInsideZ[i] + dump.coreCenterDistanceToShieldingFace() - 2*dump.frontShieldingHalfSize()[2]);

      Hep3Vector pointInMu2e(dump.beamDumpToMu2e_position(pointInDump));
      emfb->roomInsideOutline_.push_back(CLHEP::Hep2Vector(pointInMu2e.x(), pointInMu2e.z()));
    }

    // Add the wall thickness as appropriate to the inside outline to get the outside outline
    // then convert to Mu2e
    for(unsigned i=0; i<roomInsideX.size(); ++i) {

      Hep3Vector insidePointInDump
        (roomInsideX[i],
         0,
         roomInsideZ[i] + dump.coreCenterDistanceToShieldingFace() - 2*dump.frontShieldingHalfSize()[2]);

      // For ExtMonFNAL area all the angles are straight in the civil
      // drawing and the walls are parallel to the dump coordinate
      // axes, so we'll just handle this case.
      // The "regular" inside points need to be shifted in one of the four directions
      // (xz) \in { (++), (+-), (-+), (--) }
      // The first and the last points are special: we know to shift them by (-0) and (+0).

      // All the shifts are in the units of wallThickness.
      int xdir(0), zdir(0);
      if(i==0) {
        xdir = -1;
        zdir = 0;
      }
      else if(i==(roomInsideX.size()-1)) {
        xdir = +1;
        zdir = 0;
      }
      else {
        // This is a regular point - need to decide among the four directions

        // First check the "wall parallel to dump coordinate axes" assumption.
        // The doubles here are not computed but read from user input,
        // there is no rounding and the "==" comparison should work fine.
        if((roomInsideX[i-1] != roomInsideX[i]) && (roomInsideZ[i-1] != roomInsideZ[i])) {
          throw cet::exception("GEOM")<<"ExtMonFNALBuildingMaker: ERROR:"
                                      <<" the wall between the points "<<(i-1)<<" and "<<i
                                      <<" of the room outline extMonFNAL.roomInsideX, roomInsideZ"
                                      <<" is not parallel to a dump coordinate axis."
            ;
        }
        if((roomInsideX[i+1] != roomInsideX[i]) && (roomInsideZ[i+1] != roomInsideZ[i])) {
          throw cet::exception("GEOM")<<"ExtMonFNALBuildingMaker: ERROR:"
                                      <<" the wall between the points "<<i<<" and "<<(i+1)
                                      <<" of the room outline extMonFNAL.roomInsideX, roomInsideZ"
                                      <<" is not parallel to a dump coordinate axis."
            ;
        }

        // Decide the direction of the shift
        if(roomInsideZ[i] < roomInsideZ[i-1]) {
          xdir = -1;
          zdir = (roomInsideX[i+1] < roomInsideX[i]) ? +1 : -1;
        }
        else if(roomInsideZ[i] > roomInsideZ[i-1]) {
          xdir = +1;
          zdir = (roomInsideX[i+1] < roomInsideX[i]) ? +1 : -1;
        }
        else { // the (i-1)--(i) wall is parallel to the X axis
          if(roomInsideX[i] < roomInsideX[i-1]) {
            zdir = +1;
            xdir = (roomInsideZ[i+1] < roomInsideZ[i]) ? -1 : +1;
          }
          else if(roomInsideX[i] > roomInsideX[i-1]) {
            zdir = -1;
            xdir = (roomInsideZ[i+1] < roomInsideZ[i]) ? -1 : +1;
          }
        }
      } // special vs regular point

      //std::cout<<"# AG: point "<<i<<", (xdir, zdir): "<<xdir<<" "<<zdir<<std::endl;
      Hep3Vector outsidePointInDump = insidePointInDump + wallThickness * Hep3Vector(xdir, 0, zdir);

      // Convert inputs to Mu2e coordinates
      Hep3Vector outsidePointInMu2e(dump.beamDumpToMu2e_position(outsidePointInDump));
      emfb->wallOutsideOutline_.push_back(CLHEP::Hep2Vector(outsidePointInMu2e.x(), outsidePointInMu2e.z()));

      // debug:
      if(false) {
        std::cout<<"extMonFNAL_outlines"
                 <<" "<<insidePointInDump.x()
                 <<" "<<insidePointInDump.z()
                 <<" "<<outsidePointInDump.x()
                 <<" "<<outsidePointInDump.z()
                 <<" "<<emfb->roomInsideOutline_[i].x()
                 <<" "<<emfb->roomInsideOutline_[i].y()
                 <<" "<<emfb->wallOutsideOutline_[i].x()
                 <<" "<<emfb->wallOutsideOutline_[i].y()
                 <<std::endl;
      }
    } // for(points)

    //----------------------------------------------------------------
    const Hep3Vector coll2ShieldingCenterInDump(
                                                (emfb->coll2ShieldingDumpXmax() + emfb->coll2ShieldingDumpXmin())/2
                                                ,
                                                dump.backShieldingHalfSize()[1] + 0.5*emfb->roomInsideFullHeight()
                                                ,
                                                dump.coreCenterDistanceToShieldingFace() - 2*dump.frontShieldingHalfSize()[2]
                                                -(emfb->magnetRoomLength()+0.5*emfb->collimator2().horizontalLength())
                                                );

    emfb->coll2ShieldingCenterInMu2e_ = dump.beamDumpToMu2e_position(coll2ShieldingCenterInDump);
    emfb->coll2ShieldingRotationInMu2e_ = dump.coreRotationInMu2e();

    emfb->coll2ShieldingHalfSize_.resize(3);
    emfb->coll2ShieldingHalfSize_[0] = 0.5*(emfb->coll2ShieldingDumpXmax_ - emfb->coll2ShieldingDumpXmin_);
    emfb->coll2ShieldingHalfSize_[1] = 0.5*emfb->roomInsideFullHeight();
    emfb->coll2ShieldingHalfSize_[2] = 0.5*col2zLength;

    return emfb;

  } // make()

} // namespace mu2e
