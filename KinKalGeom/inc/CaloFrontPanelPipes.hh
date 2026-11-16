//
//  Cooling pipes embedded in a calorimeter disk front panel, modeled as thin-walled tubes whose centerlines lie in a
//  single z plane. KinKal can't intersect tori, so the pipes are represented by their centerlines (arcs and straight
//  segments) in that plane: a track crossing the plane within a tube radius of a centerline crosses the pipe walls.
//  The material depth along z is the chord through the tube walls at the crossing point, which is then treated as a
//  thin shell by the extrapolation.
//  original author: Sophie Middleton (2026)
//
#ifndef KinKalGeom_CaloFrontPanelPipes_hh
#define KinKalGeom_CaloFrontPanelPipes_hh
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <memory>
#include <string>
#include <vector>
namespace mu2e {
  namespace KKGeom {
    class CaloFrontPanelPipes {
      public:
        using AnnPtr = std::shared_ptr<KinKal::Annulus>;
        // circular arc centerline: points center + (xsign*rtor*cos(phi), ysign*rtor*sin(phi)) for phi in [phi0,phi0+dphi].
        // The signs reproduce the pi rotations used to place mirrored copies of the same torus.
        struct Arc {
          double cx_, cy_, xsign_, ysign_, rtor_, phi0_, dphi_;
          double rout_, rin_; // tube outer and inner radius
        };
        // straight centerline from (x0,y0) to (x1,y1)
        struct Line {
          double x0_, y0_, x1_, y1_;
          double rout_, rin_;
        };
        CaloFrontPanelPipes(SurfaceId const& sid, AnnPtr const& plane, std::string const& material) :
          sid_(sid), plane_(plane), material_(material) {}
        void addArc(Arc const& arc) { arcs_.push_back(arc); }
        void addLine(Line const& line) { lines_.push_back(line); }
        auto const& surfaceId() const { return sid_; }
        auto const& plane() const { return plane_; }
        auto const& material() const { return material_; }
        auto const& arcs() const { return arcs_; }
        auto const& lines() const { return lines_; }
        // material depth along z (mm) through the pipe walls at a position in the pipe plane; 0 if no pipe is crossed
        double wallDepth(KinKal::VEC3 const& pos) const;
      private:
        SurfaceId sid_;
        AnnPtr plane_; // plane containing the pipe centerlines, bounded by the pipe extent
        std::string material_; // pipe wall material
        std::vector<Arc> arcs_;
        std::vector<Line> lines_;
    };
  }
}
#endif
