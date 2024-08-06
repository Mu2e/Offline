// Andrei Gaponenko, 2011

#include "Offline/GeometryService/inc/ExtMonFNALBuildingMaker.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

#include "Offline/GeometryService/inc/ExtMonFNALMagnetMaker.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>
#include <numeric>

#include "cetlib_except/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  namespace {
    std::vector<double> d2r(std::vector<double> v) {
      for(auto& x : v) { x/=2; }
      return v;
    }
  }

  //================================================================
  ExtMonFNALCollimator
  ExtMonFNALBuildingMaker::readCollimatorExtMonFNAL(const std::string& name,
                                                    double angleH,
                                                    double angleV,
                                                    const SimpleConfig& c) {

    ExtMonFNALCollimator col;

    col._name = name;

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

    c.getVectorDouble("extMonFNAL."+name+".alignmentPlugInnerShellThickness", col._alignmentPlugInnerShellThickness, 2);
    c.getVectorDouble("extMonFNAL."+name+".alignmentPlugOuterShellThickness", col._alignmentPlugOuterShellThickness, 2);

    {
      std::vector<double> tmp;
      c.getVectorDouble("extMonFNAL."+name+".shotLinerInnerDiameter", tmp, 2);
      col._shotLinerInnerRadius = d2r(tmp);
    }

    c.getVectorDouble("extMonFNAL."+name+".shotLinerInnerThickness", col._shotLinerInnerThickness, 2);

    col._shotLinerOuterRadius = c.getDouble("extMonFNAL."+name+".shotLinerOuterDiameter")/2;
    col._shotLinerOuterThickness = c.getDouble("extMonFNAL."+name+".shotLinerOuterThickness");
    col._length = c.getDouble("extMonFNAL."+name+".length");
    col._radiusTransitiondZ = c.getDouble("extMonFNAL."+name+".radiusTransitiondZ");
    col._angleH_inBeamDump = angleH;
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
    ExtrudedSolid extMonRoomWall = hall.getBldgSolid("extMonExteriorWall");
    ExtrudedSolid extMonRoomWallE = hall.getBldgSolid("extMonExteriorWallE");
    const CLHEP::Hep3Vector& offset = extMonRoomWall.getOffsetFromMu2eOrigin();
    // Get corner coordinates of extinction monitor room
    const auto & roomVertices = extMonRoomWall.getVertices();
    const auto & roomVerticesE = extMonRoomWallE.getVertices();
    const double xfront = roomVertices[0][1]+offset[0];
    const double zfront = roomVertices[0][0]+offset[2];
    const double xback = roomVertices[1][1]+offset[0];
    const double zback = roomVertices[1][0]+offset[2];
    const double xScorner = roomVerticesE[1][1]+offset[0];
    const double zScorner = roomVerticesE[1][0]+offset[2];
    const double roomLength = sqrt((zfront-zback)*(zfront-zback)+(xfront-xback)*(xfront-xback));
    const double roomWidth = sqrt((zfront-zScorner)*(zfront-zScorner)+(xfront-xScorner)*(xfront-xScorner));


    emfb->roomInsideFullHeight_ = 2*extMonRoomWall.getYhalfThickness();
    emfb->roomInsideYmin_ = offset[1] - extMonRoomWall.getYhalfThickness();
    emfb->roomInsideYmax_ = emfb->roomInsideYmin_ + emfb->roomInsideFullHeight_;

    const double col2zLength = c.getDouble("extMonFNAL.collimator2.shielding.thickness");
    const double magnetRoomLength = c.getDouble("extMonFNAL.magnetRoomLength");
    const double dxdL = (xback-xfront)/roomLength;
    const double dzdL = (zback-zfront)/roomLength;
    const double shieldwidth = c.getDouble("extMonFNAL.collimator2.shielding.width");
    emfb->coll2ShieldingOutline_.emplace_back(zfront+dzdL*magnetRoomLength,
                                              xfront+dxdL*magnetRoomLength);
    emfb->coll2ShieldingOutline_.emplace_back(zfront+dzdL*(magnetRoomLength+col2zLength),
                                              xfront+dxdL*(magnetRoomLength+col2zLength));
    emfb->coll2ShieldingOutline_.emplace_back(zfront+dzdL*(magnetRoomLength+col2zLength)-shieldwidth*dxdL,
                                              xfront+dxdL*(magnetRoomLength+col2zLength)+shieldwidth*dzdL);
    emfb->coll2ShieldingOutline_.emplace_back(zfront+dzdL*magnetRoomLength-shieldwidth*dxdL,
                                              xfront+dxdL*magnetRoomLength+shieldwidth*dzdL);

    //----------------------------------------------------------------
    static CLHEP::HepRotation shieldingRot(CLHEP::HepRotation::IDENTITY);
    shieldingRot.rotateX( 90*CLHEP::degree);
    shieldingRot.rotateZ( 90*CLHEP::degree);

    emfb->coll2ShieldingRotationInMu2e_ = shieldingRot;
    emfb->coll2ShieldingCenterInMu2e_ = CLHEP::Hep3Vector(0, (emfb->roomInsideYmin()+emfb->roomInsideYmax())/2.0, 0);
    //----------------------------------------------------------------

    emfb->magnetRoomLength_ = magnetRoomLength;

    // hand stacked shielding sizes and locations in magnet room
    const double steelLength = c.getDouble("extMonFNAL.steelLength");
    const double steelwidthN = c.getDouble("extMonFNAL.steelWidthN");
    const double steelwidthS = c.getDouble("extMonFNAL.steelWidthS");
    const double layerHeight = c.getDouble("extMonFNAL.layerHeight");
    emfb->shieldingNCenterInMu2e_[0] = xfront + (dxdL*steelLength + dzdL*steelwidthN)/2.0;
    emfb->shieldingNCenterInMu2e_[1] = emfb->roomInsideYmin_ + layerHeight/2.0;
    emfb->shieldingNCenterInMu2e_[2] = zfront + (dzdL*steelLength - dxdL*steelwidthN)/2.0;
    const double safety_gap=0.1; // 100 microns
    emfb->shieldingSCenterInMu2e_[0] = xScorner + (dxdL*steelLength - dzdL*steelwidthS)/2.0+safety_gap; // SDF add safety gap
    emfb->shieldingSCenterInMu2e_[1] = emfb->roomInsideYmin_ + layerHeight/2.0;
    emfb->shieldingSCenterInMu2e_[2] = zScorner + (dzdL*steelLength + dxdL*steelwidthS)/2.0-7.*safety_gap; // SDF add safety gap

    emfb->shieldingBCenterInMu2e_[0] = xfront + (dxdL*(magnetRoomLength+steelLength) + dzdL*roomWidth)/2.0+safety_gap; // SDF add safety gap
    emfb->shieldingBCenterInMu2e_[1] = emfb->roomInsideYmin_ + layerHeight/2.0;
    emfb->shieldingBCenterInMu2e_[2] = zfront + (dzdL*(magnetRoomLength+steelLength) - dxdL*roomWidth)/2.0;

    emfb->shieldingNHalfSize_.resize(3);
    emfb->shieldingNHalfSize_[0] = steelwidthN/2.0;
    emfb->shieldingNHalfSize_[1] = layerHeight/2.0;
    emfb->shieldingNHalfSize_[2] = steelLength/2.0;

    emfb->shieldingSHalfSize_.resize(3);
    emfb->shieldingSHalfSize_[0] = steelwidthS/2.0;
    emfb->shieldingSHalfSize_[1] = layerHeight/2.0;
    emfb->shieldingSHalfSize_[2] = steelLength/2.0;

    emfb->shieldingBHalfSize_.resize(3);
    emfb->shieldingBHalfSize_[0] = roomWidth/2.0-safety_gap;// SDF add safety gap
    emfb->shieldingBHalfSize_[1] = layerHeight/2.0;
    emfb->shieldingBHalfSize_[2] = (magnetRoomLength-steelLength)/2.0-2;

    CLHEP::HepRotation nonrotation(CLHEP::HepRotation::IDENTITY);
    emfb->shieldingRotationInMu2e_ = nonrotation.rotateY(asin(dxdL) * CLHEP::radian);

    emfb->_HVACductRadius = c.getDouble("extMonFNAL.HVACductRadius");
    emfb->_HVACductHalfLength = shieldwidth/2;
    emfb->HVACductCenterInMu2e_[0] = xfront + dxdL*(magnetRoomLength+shieldwidth/2) + dzdL*emfb->_HVACductRadius;
    emfb->HVACductCenterInMu2e_[1] = emfb->roomInsideYmin_ + emfb->roomInsideFullHeight_ - emfb->_HVACductRadius;
    emfb->HVACductCenterInMu2e_[2] = zfront + dzdL*(magnetRoomLength+shieldwidth/2) - dxdL*emfb->_HVACductRadius;

    //----------------------------------------------------------------
    const double filterEntranceOffsetX = c.getDouble("extMonFNAL.entranceOffsetX") * CLHEP::mm;
    const double filterEntranceOffsetY = c.getDouble("extMonFNAL.entranceOffsetY") * CLHEP::mm;

    const double angleH = c.getDouble("extMonFNAL.angleH") * CLHEP::radian;
    const double entranceAngleV = c.getDouble("extMonFNAL.entranceAngleV") * CLHEP::radian;

    const double col1zLength = 2*dump.frontShieldingHalfSize()[2];
    emfb->filter_.collimator1_ = readCollimatorExtMonFNAL("collimator1", angleH, entranceAngleV, c);

    //----------------------------------------------------------------
    // collimator1
    const double referenceLength = dump.coreCenterDistanceToReferencePlane()    - dump.coreCenterDistanceToShieldingFace()
      + dump.frontShieldingHalfSize()[2];
    const Hep3Vector collimator1CenterInDump(filterEntranceOffsetX
                                             + referenceLength*tan(angleH),

                                             filterEntranceOffsetY
                                             + referenceLength*tan(entranceAngleV)/cos(angleH),

                                             dump.coreCenterDistanceToShieldingFace() - dump.frontShieldingHalfSize()[2]
                                             );

    emfb->filter_.collimator1_._centerInMu2e = dump.beamDumpToMu2e_position(collimator1CenterInDump);
    emfb->filter_.collimator1_.setFromDumpAngles(angleH, entranceAngleV, dump);

    const Hep3Vector collimator1ExitInDump(filterEntranceOffsetX
                                           + 1.*col1zLength*tan(emfb->filter().collimator1().angleH_inBeamDump()),

                                           filterEntranceOffsetY
                                           + 1.*col1zLength*tan(emfb->filter().collimator1().angleV())/cos(emfb->filter().collimator1().angleH_inBeamDump()),

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
    emfb->filter_.nominalMomentum_ = pNominal;
    emfb->filter_.magnet_ = ExtMonFNALMagnetMaker::read(c, "extMonFNAL.filter.magnet",
                                                        emfb->filter().collimator1().rotationInMu2e(),
                                                        dump.beamDumpToMu2e_position(refTrajFMEntranceInDump),
                                                        pNominal);

    const Hep3Vector magnetRefInDump = dump.mu2eToBeamDump_position(emfb->filter_.magnet_.refPointInMu2e());

    //----------------------------------------------------------------
    // collimator2

    emfb->filter_.collimator2_ = readCollimatorExtMonFNAL("collimator2",
                                                          angleH,
                                                          entranceAngleV - 2 * emfb->filter_.magnet_.trackBendHalfAngle(pNominal),
                                                          c);

    // We want to position the collimator inside the concrete created from
    // the "coll2ShieldingOutline_" above.    The outline vertices
    // are coordinates in the Mu2e (z,x) plane.   Find their 2D geometric
    // center, then bring it into 3D and convert from Mu2e to dump system coordinates.

    const auto coll2cog2d = std::accumulate(emfb->coll2ShieldingOutline_.begin(),
                                            emfb->coll2ShieldingOutline_.end(),
                                            Hep2Vector(0,0))
      / emfb->coll2ShieldingOutline_.size();

    // point in Mu2e coordinates
    const Hep3Vector coll2cog3d = Hep3Vector(coll2cog2d.y(), 0, coll2cog2d.x())
      + emfb->coll2ShieldingCenterInMu2e_ ;

    // point in dump coordinates
    const auto coll2cogInDump = dump.mu2eToBeamDump_position(coll2cog3d);

    const double col2CenterZdump = coll2cogInDump.z();

    const double col2CenterXdump = magnetRefInDump[0]
      +  (magnetRefInDump[2] - col2CenterZdump)*tan(emfb->filter().collimator1().angleH_inBeamDump());

    const double col2CenterYdump = magnetRefInDump[1]
      +  (magnetRefInDump[2] - col2CenterZdump)*tan(emfb->filter().collimator2().angleV())/cos(emfb->filter().collimator1().angleH_inBeamDump());

    emfb->filter_.collimator2_._centerInMu2e = dump.beamDumpToMu2e_position(Hep3Vector(col2CenterXdump, col2CenterYdump, col2CenterZdump));
    emfb->filter_.collimator2_._rotationInMu2e = emfb->filter().magnet().outRotationInMu2e();

    //----------------------------------------------------------------
    if(verbose) {
      std::cout<<"ExtMonFNALBuildingMaker"<<": Filter entrance offsets  = ("<<filterEntranceOffsetX<<", "<<filterEntranceOffsetY<<")"<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filter nominal momentum = "<<emfb->filter().magnet().nominalMomentum()/CLHEP::GeV<<" GeV/c"<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filterAngleH = "<<emfb->filter().collimator1().angleH_inBeamDump()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filterAngleH in Mu2e, degrees= "<<(dump.coreRotY() - emfb->filter().collimator1().angleH_inBeamDump())/CLHEP::degree<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filter half bend angle  = "<<emfb->filter().magnet().trackBendHalfAngle(emfb->filter().nominalMomentum())<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filter.angleV = "<<emfb->filter().collimator1().angleV()
               <<", c1.angleV  = "<<emfb->filter().collimator1().angleV()
               <<", c2.angleV() = "<<emfb->filter().collimator2().angleV()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator1CenterInMu2e = "<<emfb->filter().collimator1().centerInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator1 exit in Mu2e = "<< dump.beamDumpToMu2e_position(collimator1ExitInDump)<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": ref traj entrace to filter magnet in Mu2e = "<< dump.beamDumpToMu2e_position(refTrajFMEntranceInDump)<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": filterMagnet().refPointInMu2e() = "<<emfb->filter().magnet().refPointInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": filterMagnet().geometricCenterInMu2e() = "<<emfb->filter().magnet().geometricCenterInMu2e()<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": filterMagnet().magnetRotationInMu2e() = "<<emfb->filter().magnet().magnetRotationInMu2e()<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": collimator2CenterInMu2e = "<<emfb->filter().collimator2().centerInMu2e()<<std::endl;

      std::cout<<"ExtMonFNALBuildingMaker"<<": ExtMonFNALBuilding::filterEntranceInMu2e() = "<<emfb->filter().entranceInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALBuildingMaker"<<": ExtMonFNALBuilding::filterExitInMu2e() = "<<emfb->filter().exitInMu2e()<<std::endl;
    }

    //----------------------------------------------------------------
    // Detector room box

    Hep2Vector corner1(emfb->coll2ShieldingOutline_[1] + Hep2Vector(emfb->coll2ShieldingCenterInMu2e_[2], emfb->coll2ShieldingCenterInMu2e_[0]));
    Hep2Vector corner2(extMonRoomWall.getVertices()[1] + Hep2Vector(offset[2], offset[0]));
    Hep2Vector corner3(extMonRoomWall.getVertices()[2] + Hep2Vector(offset[2], offset[0]));

    // Sanity check.  Did we use the correct indexes to get room corner, which should be 90 degrees?
    //
    // The default tolerance is order 10^{-14} while we have a worse than 10^{-6} precision
    // on the concrete coordinates.   So we have have to relax the tolerance.
    Hep2Vector dv21(corner2-corner1);
    dv21.setTolerance(1.e-4);
    if(!dv21.isOrthogonal(corner3-corner2)) {
      throw cet::exception("GEOM")<<"ExtMonFNALBuildingMaker::make(): "
                                  <<" Error: unexpected non-right angle when computing ExtMon detector geometry: "
                                  <<" non-orthogonalizy = "<<dv21.howOrthogonal(corner3-corner2)
                                  <<"\n";
    }

    // Compute detector room parameters.  Due to the finite floating
    // point precision we get spurious volume overlaps when using the
    // exact box size and position.  Use an ad-hoc margin to slightly
    // shrink the box we place.
    const double detectorRoomMargin = 0.1*CLHEP::mm;
    emfb->detectorRoomHalfSize_ = std::vector<double>{
      0.5*(corner3-corner2).mag() - detectorRoomMargin,
      extMonRoomWall.getYhalfThickness() - detectorRoomMargin,
      0.5*(corner2-corner1).mag() - detectorRoomMargin
    };

    Hep2Vector detectorRoomCenter((corner1 + corner3)/2);
    emfb->detectorRoomCenterInMu2e_ = Hep3Vector(detectorRoomCenter[1] , (emfb->roomInsideYmin()+emfb->roomInsideYmax())/2, detectorRoomCenter[0]);

    emfb->detectorRoomRotationInMu2e_.rotateY(asin( (corner1[1]-corner2[1])/(corner1-corner2).mag() ));

    //----------------------------------------------------------------
    return emfb;

  } // make()

} // namespace mu2e
