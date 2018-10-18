// Andrei Gaponenko, 2012

#include "GeometryService/inc/Mu2eEnvelope.hh"

#include <limits>
#include <algorithm>
#include <iterator>
#include <utility>
#include <cmath>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/TwoVector.h"

#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

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
    : xmin_(), xmax_(), ymin_(), ymax_(), zmin_(), zmax_()
  {}

  //----------------------------------------------------------------
  Mu2eEnvelope::Mu2eEnvelope(const Mu2eHall& building,
			     const SimpleConfig& config)
    : xmin_(), xmax_(), ymin_(), ymax_(), zmin_(), zmax_()
  {

    // Determine envelope from building solids
    for ( const auto& volPair : building.getBldgSolids() ) {
      const ExtrudedSolid& vol = volPair.second;
      
      const auto xMinMax = getMinMax( vol.getVertices(), &CLHEP::Hep2Vector::y );
      xmin_ = std::min( xmin_, xMinMax.first  + vol.getOffsetFromMu2eOrigin().x() );
      xmax_ = std::max( xmax_, xMinMax.second + vol.getOffsetFromMu2eOrigin().x() );

      const auto zMinMax = getMinMax( vol.getVertices(), &CLHEP::Hep2Vector::x );
      zmin_ = std::min( zmin_, zMinMax.first  + vol.getOffsetFromMu2eOrigin().z() );
      zmax_ = std::max( zmax_, zMinMax.second + vol.getOffsetFromMu2eOrigin().z() );

      const double yCenter = vol.getOffsetFromMu2eOrigin().y();
      ymin_ = std::min( ymin_ , yCenter-vol.getYhalfThickness() );
      ymax_ = std::max( ymax_ , yCenter+vol.getYhalfThickness() );      

    }
     
    // Add extra room for dirt margins -- cannot be zero!

    if ( config.getDouble( "world.dirt.minimalMargin.zmin" ) <= 0. ||
         config.getDouble( "world.dirt.minimalMargin.zmax" ) <= 0. || 
         config.getDouble( "world.dirt.minimalMargin.xmin" ) <= 0. || 
         config.getDouble( "world.dirt.minimalMargin.xmax" ) <= 0. ) {
      throw cet::exception("GEOM") << "All world margins must be greater than 0.";
    }


    zmin_ -= config.getDouble( "world.dirt.minimalMargin.zmin" );
    zmax_ += config.getDouble( "world.dirt.minimalMargin.zmax" );

    xmin_ -= config.getDouble( "world.dirt.minimalMargin.xmin" );
    xmax_ += config.getDouble( "world.dirt.minimalMargin.xmax" );

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
