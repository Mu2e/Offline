#include "Offline/KinKalGeom/inc/CaloFrontPanelPipes.hh"
#include <algorithm>
#include <cmath>
namespace mu2e {
  namespace KKGeom {
    namespace {
      // chord length through the walls of a tube whose axis is at transverse distance dist from a line along z
      double wallChord(double dist, double rout, double rin) {
        if(dist >= rout) return 0.0;
        double chord = 2.0*std::sqrt(rout*rout - dist*dist);
        if(dist < rin) chord -= 2.0*std::sqrt(rin*rin - dist*dist);
        return chord;
      }
    }

    double CaloFrontPanelPipes::wallDepth(KinKal::VEC3 const& pos) const {
      // the pipes don't overlap, so the closest centerline is the only one that can be crossed
      double depth(0.0);
      for(auto const& arc : arcs_){
        // transform into the frame of the unmirrored arc
        double qx = arc.xsign_*(pos.X() - arc.cx_);
        double qy = arc.ysign_*(pos.Y() - arc.cy_);
        double phi = std::atan2(qy,qx) - arc.phi0_;
        phi -= 2.0*M_PI*std::floor(phi/(2.0*M_PI)); // [0,2pi)
        double dist;
        if(phi <= arc.dphi_){
          dist = std::fabs(std::sqrt(qx*qx+qy*qy) - arc.rtor_);
        } else { // beyond the arc ends: distance to the nearest end
          double dx0 = qx - arc.rtor_*std::cos(arc.phi0_), dy0 = qy - arc.rtor_*std::sin(arc.phi0_);
          double dx1 = qx - arc.rtor_*std::cos(arc.phi0_+arc.dphi_), dy1 = qy - arc.rtor_*std::sin(arc.phi0_+arc.dphi_);
          dist = std::sqrt(std::min(dx0*dx0+dy0*dy0,dx1*dx1+dy1*dy1));
        }
        depth = std::max(depth,wallChord(dist,arc.rout_,arc.rin_));
      }
      for(auto const& line : lines_){
        double lx = line.x1_ - line.x0_, ly = line.y1_ - line.y0_;
        double frac = std::clamp(((pos.X()-line.x0_)*lx + (pos.Y()-line.y0_)*ly)/(lx*lx+ly*ly),0.0,1.0);
        double dx = pos.X() - (line.x0_ + frac*lx), dy = pos.Y() - (line.y0_ + frac*ly);
        depth = std::max(depth,wallChord(std::sqrt(dx*dx+dy*dy),line.rout_,line.rin_));
      }
      return depth;
    }
  }
}
