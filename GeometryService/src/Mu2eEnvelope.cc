// Andrei Gaponenko, 2012

#include "GeometryService/inc/Mu2eEnvelope.hh"

#include <limits>
#include <algorithm>
#include <iterator>
#include <utility>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Vector/TwoVector.h"

#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"

namespace mu2e {

  namespace {
    // a helper to compute the bounds of an extruded solid
    template<class Point>
    std::pair<double,double> getMinMax(const std::vector<Point>& outline, double (Point::*coord)() const) {
      typename std::vector<Point>::const_iterator i = outline.begin();
      assert(i != outline.end());
      double a = ((*i).*coord)();
      double b = a;
      for(; i!=outline.end(); ++i) {
        a = std::min(a, ((*i).*coord)());
        b = std::max(b, ((*i).*coord)());
      }
      return std::make_pair(a, b);
    }
  }

  //----------------------------------------------------------------
  Mu2eEnvelope::Mu2eEnvelope()
    : xmin_(0), xmax_(0), ymin_(0), ymax_(0), zmin_(0), zmax_(0)
  {}

  //----------------------------------------------------------------
  Mu2eEnvelope::Mu2eEnvelope(const Mu2eBuilding& building,
                             const ProtonBeamDump& dump,
                             const ExtMonFNALBuilding& emfb)
    : xmin_(0), xmax_(0), ymin_(0), ymax_(0), zmin_(0), zmax_(0)
  {
    const std::pair<double,double> emfXlimits = getMinMax(emfb.wallOutsideOutline(), &CLHEP::Hep2Vector::x);
    const std::pair<double,double> emfZlimits = getMinMax(emfb.wallOutsideOutline(), &CLHEP::Hep2Vector::y);

    zmin_ =
      std::min(
               building.hallInsideZExtMonUCIWall() - building.hallWallThickness()
               ,
               emfZlimits.first - emfb.dirtOverheadHorizontalMargin()
               );

    zmax_ = building.hallInsideZmax() + building.hallWallThickness();

    xmin_ =
      std::min(
               building.hallInsideXmin() - building.hallWallThickness()
               ,
               emfXlimits.first - emfb.dirtOverheadHorizontalMargin()
               );

    xmax_ =
      std::max(
               building.hallInsideXmax() + building.hallWallThickness()
               ,
               emfXlimits.second + emfb.dirtOverheadHorizontalMargin()
               );

    ymin_ =
      std::min(
               building.hallInsideYmin() - building.hallFloorThickness()
               ,
               emfb.roomInsideYmin() - emfb.roomFloorThickness()
               );

    const double dirtFormalYmax  = std::max(
                                            emfb.roomInsideYmax() + emfb.roomCeilingThickness() + emfb.dirtOverheadThickness()
                                            ,
                                            dump.frontShieldingCenterInMu2e()[1] + dump.frontShieldingHalfSize()[1]
                                            );

    ymax_ =
      std::max(
               dirtFormalYmax
               ,
               building.hallInsideYmax() + building.hallCeilingThickness()
               // + building.dirtOverburdenDepth() + 2*building.dirtCapHalfHeight()
               );
  }

  //================================================================
  std::ostream& operator<<(std::ostream& os, const Mu2eEnvelope& env) {
    return os<<"Mu2eEnvelope(xmin="<<env.xmin()
	     <<",xmax="<<env.xmax()
	     <<",ymin="<<env.ymin()
	     <<",ymax="<<env.ymax()
	     <<",zmin="<<env.zmin()
	     <<",zmax="<<env.zmax()
	     <<" )";
  }

  //================================================================

} // namespace mu2e
