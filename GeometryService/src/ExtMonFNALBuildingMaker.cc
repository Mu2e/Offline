// Andrei Gaponenko, 2011

#include "Offline/GeometryService/inc/ExtMonFNALBuildingMaker.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

#include "Offline/GeometryService/inc/ExtMonFNALMagnetMaker.hh"
#include "Offline/GeometryService/inc/ExtMonFNALFilterMaker.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>
#include <numeric>

#include "cetlib_except/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  //================================================================
  std::unique_ptr<ExtMonFNALBuilding> ExtMonFNALBuildingMaker::make(const SimpleConfig& c,
                                                                    const Mu2eHall& hall,
                                                                    const ProtonBeamDump& dump) {
    using CLHEP::Hep3Vector;
    using CLHEP::Hep2Vector;

    std::unique_ptr<ExtMonFNALBuilding> emfb(new ExtMonFNALBuilding());

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
    // FilterMaker needs information on coll2ShieldingOutline, call it
    // after it is initialized
    emfb->filter_ = ExtMonFNALFilterMaker::read(c, *emfb, dump);

    //----------------------------------------------------------------
    return emfb;

  } // make()

} // namespace mu2e
