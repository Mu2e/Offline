// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuildingMaker.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnetMaker.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  namespace {
    std::vector<double> d2r(std::vector<double> v) {
      for(auto& x : v) { x/=2; }
      return v;
    }
  }

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

    {
      std::vector<double> tmp;
      c.getVectorDouble("extMonFNAL."+name+".channelDiameter", tmp, 2);
      col._channelRadius = d2r(tmp);
    }

    {
      std::vector<double> tmp;
      c.getVectorDouble("extMonFNAL."+name+".alignmentPlugDiameter", tmp, 2);
      col._alignmentPlugRadius = d2r(tmp);
    }

    {
      std::vector<double> tmp;
      c.getVectorDouble("extMonFNAL."+name+".alignmentHoleDiameter", tmp, 2);
      col._alignmentHoleRadius = d2r(tmp);
    }

    col._radiusTransitiondZ = c.getDouble("extMonFNAL."+name+".radiusTransitiondZ");
    col._angleH = angleH;
    col._angleV = angleV;

    return col;
  }

  //================================================================
  std::unique_ptr<ExtMonFNALBuilding> ExtMonFNALBuildingMaker::make(const SimpleConfig& c,
                                                	    	    const Mu2eHall& hall,
								    const ProtonBeamDump& dump) {
    using CLHEP::Hep3Vector;
    using CLHEP::Hep2Vector;

    std::unique_ptr<ExtMonFNALBuilding> emfb(new ExtMonFNALBuilding());

    int verbose = c.getInt("extMonFNAL.verbosityLevel");
 
  // Get relevant Hall solid
    ExtrudedSolid extMonRoom = hall.getBldgSolid("extMon");
    const CLHEP::Hep3Vector& offset = extMonRoom.getOffsetFromMu2eOrigin();
  // Get front and back corner coordinates of north wall
    const auto & roomVertices = extMonRoom.getVertices();
    const double zfront = roomVertices[0][0]+offset[2];
    const double zback = roomVertices[1][0]+offset[2];
    const double xfront = roomVertices[0][1]+offset[0];
    const double xback = roomVertices[1][1]+offset[0];
    const double roomLength = sqrt((zfront-zback)*(zfront-zback)+(xfront-xback)*(xfront-xback));


    emfb->roomInsideFullHeight_ = 2*extMonRoom.getYhalfThickness();
    emfb->roomInsideYmin_ = offset[1] - extMonRoom.getYhalfThickness();
    emfb->roomInsideYmax_ = emfb->roomInsideYmin_ + emfb->roomInsideFullHeight_;

    const double col2zLength = c.getDouble("extMonFNAL.collimator2.shielding.thickness");
    const double magnetRoomLength = c.getDouble("extMonFNAL.magnetRoomLength");
    const double dxdL = (xback-xfront)/roomLength;
    const double dzdL = (zback-zfront)/roomLength;
    const double shieldwidth = 2*dump.backShieldingHalfSize()[1];
    emfb->coll2ShieldingOutline_.emplace_back(zfront+dzdL*magnetRoomLength,
					      xfront+dxdL*magnetRoomLength);
    emfb->coll2ShieldingOutline_.emplace_back(zfront+dzdL*(magnetRoomLength+col2zLength),
					      xfront+dxdL*(magnetRoomLength+col2zLength));
    emfb->coll2ShieldingOutline_.emplace_back(zfront+dzdL*(magnetRoomLength+col2zLength)+shieldwidth*sin(dump.coreRotY()),
					      xfront+dxdL*(magnetRoomLength+col2zLength)-shieldwidth*cos(dump.coreRotY()));
    emfb->coll2ShieldingOutline_.emplace_back(zfront+dzdL*magnetRoomLength+shieldwidth*sin(dump.coreRotY()),
					      xfront+dxdL*magnetRoomLength-shieldwidth*cos(dump.coreRotY()));

    emfb->magnetRoomLength_ = magnetRoomLength;
    emfb->coll2ShieldingDumpXmin_ = c.getDouble("extMonFNAL.collimator2.shielding.dumpXmin");
    emfb->coll2ShieldingDumpXmax_ = c.getDouble("extMonFNAL.collimator2.shielding.dumpXmax");

    //----------------------------------------------------------------
    // We need to rotate the extruded solid coordinates by +90 degrees to bring it to the horizontal plane.
    // No other rotation is needed since the outline is already in Mu2e (x,z), i.e. rorated with the dump
//    emfb->roomRefPointInMu2e_ = Hep3Vector(0, (emfb->roomInsideYmax() + emfb->roomInsideYmin())/2, 0);
//    emfb->roomRotationInMu2e_ = CLHEP::HepRotationX(+90.*CLHEP::degree);

    //----------------------------------------------------------------
    emfb->_filterEntranceOffsetX = c.getDouble("extMonFNAL.entranceOffsetX") * CLHEP::mm;
    emfb->_filterEntranceOffsetY = c.getDouble("extMonFNAL.entranceOffsetY") * CLHEP::mm;

    const double angleH = emfb->_filterAngleH =  c.getDouble("extMonFNAL.angleH") * CLHEP::radian;
    const double entranceAngleV = emfb->_filterEntranceAngleV =  c.getDouble("extMonFNAL.entranceAngleV") * CLHEP::radian;

    const double col1zLength = 2*dump.frontShieldingHalfSize()[2];
    emfb->_collimator1 = readCollimatorExtMonFNAL("collimator1", col1zLength, angleH, entranceAngleV, c);

    //----------------------------------------------------------------
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

    const Hep3Vector collimator1ExitInDump(emfb->filterEntranceOffsetX()
                                           + 1.*col1zLength*tan(emfb->collimator1().angleH()),

                                           emfb->filterEntranceOffsetY()
                                           + 1.*col1zLength*tan(emfb->collimator1().angleV())/cos(emfb->filterAngleH()),

                                           dump.coreCenterDistanceToShieldingFace() - col1zLength
                                           );

    //----------------------------------------------------------------
    // Filter magnet positioning

    // the distance between the exit point of the reference trajectory from the upstream wall
    // and its entrance to the magnet, on the magnet physical face.
    const double filterMagToColl = c.getDouble("extMonFNAL.filter.magnet.distanceToUpstreamWall")
      / (cos(angleH) * cos(entranceAngleV));

    // The point where the reference trajectory crosses the magnet face, in the dump coordinates
    const CLHEP::Hep3Vector refTrajFMEntranceInDump = collimator1ExitInDump +
      Hep3Vector(0,0, -filterMagToColl).rotateX(entranceAngleV).rotateY(-angleH);

    const double pNominal = c.getDouble("extMonFNAL.filter.nominalMomentum") * CLHEP::MeV;
    emfb->_filterMagnet = ExtMonFNALMagnetMaker::read(c, "extMonFNAL.filter.magnet",
                                                      emfb->_collimator1RotationInMu2e,
                                                      dump.beamDumpToMu2e_position(refTrajFMEntranceInDump),
                                                      pNominal);

    const Hep3Vector magnetRefInDump = dump.mu2eToBeamDump_position(emfb->_filterMagnet.refPointInMu2e());

    //----------------------------------------------------------------
    // collimator2

   emfb->_collimator2 = readCollimatorExtMonFNAL("collimator2",
                                                  col2zLength,
                                                  angleH,
                                                  entranceAngleV - 2 * emfb->_filterMagnet.trackBendHalfAngle(pNominal),
                                                  c);

    const double col2CenterZdump = dump.coreCenterDistanceToShieldingFace()
      - 2*dump.frontShieldingHalfSize()[2]
      - emfb->magnetRoomLength_
      - 0.5*emfb->_collimator2.horizontalLength()
      ;

    const double col2CenterXdump = magnetRefInDump[0]
      +  (magnetRefInDump[2] - col2CenterZdump)*tan(emfb->filterAngleH());

    const double col2CenterYdump = magnetRefInDump[1]
      +  (magnetRefInDump[2] - col2CenterZdump)*tan(emfb->_collimator2.angleV())/cos(emfb->filterAngleH());

    emfb->_collimator2CenterInMu2e = dump.beamDumpToMu2e_position(Hep3Vector(col2CenterXdump, col2CenterYdump, col2CenterZdump));
    emfb->_collimator2RotationInMu2e = emfb->filterMagnet().outRotationInMu2e();

    //----------------------------------------------------------------
    if(verbose) {
      std::cout<<"ExtMonFNALBuildingMaker"<<": Filter entrance offsets  = ("<<emfb->_filterEntranceOffsetX<<", "<<emfb->_filterEntranceOffsetY<<")"<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filter nominal momentum = "<<emfb->filterMagnet().nominalMomentum()/CLHEP::GeV<<" GeV/c"<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filterAngleH = "<<emfb->filterAngleH()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filterAngleH in Mu2e, degrees= "<<(dump.coreRotY() - emfb->filterAngleH())/CLHEP::degree<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filter half bend angle  = "<<emfb->filterMagnet().trackBendHalfAngle(emfb->filterMagnet().nominalMomentum())<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filter.angleV = "<<emfb->filterEntranceAngleV()
               <<", c1.angleV  = "<<emfb->collimator1().angleV()
               <<", c2.angleV() = "<<emfb->collimator2().angleV()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator1CenterInMu2e = "<<emfb->_collimator1CenterInMu2e<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator1.horizontalLength = "<<emfb->_collimator1._horizontalLength<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator1 exit in Mu2e = "<< dump.beamDumpToMu2e_position(collimator1ExitInDump)<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": ref traj entrace to filter magnet in Mu2e = "<< dump.beamDumpToMu2e_position(refTrajFMEntranceInDump)<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": filterMagnet().refPointInMu2e() = "<<emfb->_filterMagnet.refPointInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filterMagnet().geometricCenterInMu2e() = "<<emfb->_filterMagnet.geometricCenterInMu2e()<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": filterMagnet().magnetRotationInMu2e() = "<<emfb->_filterMagnet.magnetRotationInMu2e()<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator2CenterInMu2e = "<<emfb->_collimator2CenterInMu2e<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator2.horizontalLength = "<<emfb->_collimator2._horizontalLength<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": ExtMonFNALBuilding::filterEntranceInMu2e() = "<<emfb->filterEntranceInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": ExtMonFNALBuilding::filterExitInMu2e() = "<<emfb->filterExitInMu2e()<<std::endl;
    }
/*
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

    if(roomInsideZ[0] != 0.) {
      throw cet::exception("GEOM")<<"ExtMonFNALBuildingMaker: ERROR: roomInsideZ[0] must be 0";
    }
    if(roomInsideZ.back() != 0.) {
      throw cet::exception("GEOM")<<"ExtMonFNALBuildingMaker: ERROR: the last entry in roomInsideZ must be 0";
    }

    // Convert inside outline to Mu2e coords
    for(unsigned i=0; i<roomInsideX.size(); ++i) {
      // Convert inputs to Mu2e coordinates
      Hep3Vector pointInDump(roomInsideX[i],
                             0,
                             roomInsideZ[i] + dump.coreCenterDistanceToShieldingFace() - 2*dump.frontShieldingHalfSize()[2]);

      Hep3Vector pointInMu2e(dump.beamDumpToMu2e_position(pointInDump));
      emfb->roomInsideOutline_.push_back(Hep2Vector(pointInMu2e.x(), pointInMu2e.z()));
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
      emfb->wallOutsideOutline_.push_back(Hep2Vector(outsidePointInMu2e.x(), outsidePointInMu2e.z()));

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
    // Floor outline differs from the ceiling because the back dump
    // core back shielding protrudes there.
    //
    // Need to add new first/last points, or shift the existing ones.
    if(true) {
      emfb->floorOutsideOutline_ = emfb->wallOutsideOutline_;

      Hep3Vector shiftInDump(0, 0, -2*dump.backShieldingHalfSize()[2]);
      Hep3Vector shiftInMu2e = dump.coreRotationInMu2e() * shiftInDump;
      Hep2Vector shift(shiftInMu2e.x(), shiftInMu2e.z());

      Hep3Vector fsRef(0,
                       0,
                       dump.coreCenterDistanceToShieldingFace() - 2*dump.frontShieldingHalfSize()[2]
                       );

      // Start of the path:
      if(roomInsideX[0] - wallThickness < -dump.backShieldingHalfSize()[0]) { // Need to add points

        Hep3Vector startFirst = dump.beamDumpToMu2e_position
          (fsRef + Hep3Vector(-dump.backShieldingHalfSize()[0], 0, -2*dump.backShieldingHalfSize()[2]));

        Hep3Vector startSecond =  dump.beamDumpToMu2e_position
          (fsRef + Hep3Vector(-dump.backShieldingHalfSize()[0], 0, 0));

        emfb->floorOutsideOutline_.insert(emfb->floorOutsideOutline_.begin(),
                                          Hep2Vector(startSecond.x(), startSecond.z()));

        emfb->floorOutsideOutline_.insert(emfb->floorOutsideOutline_.begin(),
                                          Hep2Vector(startFirst.x(), startFirst.z()));
      }
      else { // need to shift
        emfb->floorOutsideOutline_.front() += shift;
      }

      // End of the path:
      if(roomInsideX.back() + wallThickness > +dump.backShieldingHalfSize()[0]) { // Need to add points
        Hep3Vector endLast = dump.beamDumpToMu2e_position
          (fsRef + Hep3Vector(+dump.backShieldingHalfSize()[0], 0, -2*dump.backShieldingHalfSize()[2]));

        Hep3Vector endSecond =  dump.beamDumpToMu2e_position
          (fsRef + Hep3Vector(+dump.backShieldingHalfSize()[0], 0, 0));

        emfb->floorOutsideOutline_.push_back(Hep2Vector(endSecond.x(), endSecond.z()));
        emfb->floorOutsideOutline_.push_back(Hep2Vector(endLast.x(), endLast.z()));
      }
      else { // need to shift
        emfb->floorOutsideOutline_.back()  += shift;
      }
    }

    // debug:
    if(false) {
      for(unsigned i = 0; i < emfb->floorOutsideOutline_.size(); ++i) {
        std::cout<<"extMonFNAL_floor "
                 <<" "<<emfb->floorOutsideOutline_[i].x()
                 <<" "<<emfb->floorOutsideOutline_[i].y()
                 <<std::endl;
      }
    }
*/
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
