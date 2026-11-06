#include "Offline/KinKalGeom/inc/KinKalGeomCheck.hh"
#include "Offline/KinKalGeom/inc/DetectorSolenoid.hh"
#include "Offline/KinKalGeom/inc/CRV.hh"
#include "KinKal/Geometry/Rectangle.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <algorithm>
#include <cmath>
#include <utility>

namespace mu2e {
  namespace KKGeomCheck {
    namespace {
      // intersect [a0,a1] with [b0,b1] in place; false if empty
      bool clip(double& a0, double& a1, double b0, double b1) { a0 = std::max(a0,b0); a1 = std::min(a1,b1); return a1 > a0; }
      // interval where |c + s*d| < h along one axis
      bool slabInterval(double c, double d, double h, double& s0, double& s1) {
        if(std::fabs(d) < 1.0e-12) return std::fabs(c) < h;
        double t0 = (-h-c)/d, t1 = (h-c)/d;
        if(t0 > t1) std::swap(t0,t1);
        return clip(s0,s1,t0,t1);
      }
      // interval where the radial distance from the axis is below r
      bool radialInterval(VEC3 const& q, VEC3 const& dq, double r, double& s0, double& s1) {
        double const a = dq.Mag2(), b = q.Dot(dq), c = q.Mag2() - r*r;
        if(a < 1.0e-12) return c < 0.0;
        double const disc = b*b - a*c;
        if(disc <= 0.0) return false;
        double const sq = std::sqrt(disc);
        return clip(s0,s1,(-b-sq)/a,(-b+sq)/a);
      }
    }

    bool Element::inside(VEC3 const& point) const {
      VEC3 const d = point - center;
      double const w = d.Dot(wdir);
      if(std::fabs(w) > whalf) return false;
      if(shape == box) return std::fabs(d.Dot(udir)) <= uhalf && std::fabs(d.Dot(vdir)) <= vhalf;
      double const r = (d - w*wdir).R();
      return r >= rin && r <= rout;
    }

    VEC3 Element::sample(std::mt19937& rng) const {
      std::uniform_real_distribution<double> flat(-1.0,1.0);
      if(shape == box) return center + (uhalf*flat(rng))*udir + (vhalf*flat(rng))*vdir + (whalf*flat(rng))*wdir;
      double const r = std::sqrt(rin*rin + 0.5*(1.0+flat(rng))*(rout*rout-rin*rin));
      double const phi = M_PI*flat(rng);
      return center + (r*std::cos(phi))*udir + (r*std::sin(phi))*vdir + (whalf*flat(rng))*wdir;
    }

    double Element::volume() const {
      return shape == box ? 8.0*uhalf*vhalf*whalf : 2.0*M_PI*(rout*rout-rin*rin)*whalf;
    }

    void Element::bounds(VEC3& lo, VEC3& hi) const {
      VEC3 ext;
      if(shape == box){
        ext = VEC3(uhalf*std::fabs(udir.X()) + vhalf*std::fabs(vdir.X()) + whalf*std::fabs(wdir.X()),
                   uhalf*std::fabs(udir.Y()) + vhalf*std::fabs(vdir.Y()) + whalf*std::fabs(wdir.Y()),
                   uhalf*std::fabs(udir.Z()) + vhalf*std::fabs(vdir.Z()) + whalf*std::fabs(wdir.Z()));
      } else {
        auto rad = [this](double a){ return rout*std::sqrt(std::max(0.0,1.0-a*a)) + whalf*std::fabs(a); };
        ext = VEC3(rad(wdir.X()),rad(wdir.Y()),rad(wdir.Z()));
      }
      lo = center - ext; hi = center + ext;
    }

    bool cylinderInterval(VEC3 const& center, VEC3 const& axis, double radius, double halflen,
        VEC3 const& start, VEC3 const& dir, double smin, double smax, double& s0, double& s1) {
      s0 = smin; s1 = smax;
      VEC3 const d = start - center;
      double const z = d.Dot(axis), dz = dir.Dot(axis);
      if(!slabInterval(z,dz,halflen,s0,s1)) return false;
      return radialInterval(d - z*axis, dir - dz*axis, radius, s0, s1);
    }

    double Element::pathLength(VEC3 const& start, VEC3 const& dir, double smin, double smax) const {
      VEC3 const d = start - center;
      double s0 = smin, s1 = smax;
      if(shape == box){
        if(!slabInterval(d.Dot(udir),dir.Dot(udir),uhalf,s0,s1)) return 0.0;
        if(!slabInterval(d.Dot(vdir),dir.Dot(vdir),vhalf,s0,s1)) return 0.0;
        if(!slabInterval(d.Dot(wdir),dir.Dot(wdir),whalf,s0,s1)) return 0.0;
        return s1 - s0;
      }
      double const z = d.Dot(wdir), dz = dir.Dot(wdir);
      if(!slabInterval(z,dz,whalf,s0,s1)) return 0.0;
      VEC3 const q = d - z*wdir, dq = dir - dz*wdir;
      double o0 = s0, o1 = s1;
      if(!radialInterval(q,dq,rout,o0,o1)) return 0.0;
      double len = o1 - o0;
      double i0 = o0, i1 = o1;
      if(rin > 0.0 && radialInterval(q,dq,rin,i0,i1)) len -= i1 - i0;
      return len;
    }

    Elements elements(KinKalGeom const& kkg) {
      Elements elems;
      if(kkg.DS()){
        for(auto const& cyl : kkg.DS()->materialCylinders()){
          Element el;
          el.name = cyl.sid_.name();
          el.material = cyl.material_;
          el.kind = Element::dsCylinder;
          el.shape = Element::tube;
          el.center = cyl.surface_->center();
          el.wdir = cyl.surface_->axis().Unit();
          el.udir = (std::fabs(el.wdir.X()) < 0.9 ? VEC3(1.0,0.0,0.0) : VEC3(0.0,1.0,0.0)).Cross(el.wdir).Unit();
          el.vdir = el.wdir.Cross(el.udir);
          el.whalf = cyl.surface_->halfLength();
          el.rin = cyl.surface_->radius() - 0.5*cyl.thickness_;
          el.rout = cyl.surface_->radius() + 0.5*cyl.thickness_;
          elems.push_back(el);
        }
      }
      auto const& planes = kkg.passiveMaterialPlanes();
      for(size_t ip = 0; ip < planes.size(); ++ip){
        auto const& plane = planes[ip];
        Element el;
        el.name = plane.sid_.name() + "_" + std::to_string(plane.sid_.index());
        el.material = plane.material_;
        el.kind = Element::passivePlane;
        el.center = plane.surface_->center();
        el.udir = plane.surface_->uDirection();
        el.vdir = plane.surface_->vDirection();
        el.wdir = plane.surface_->normal();
        el.whalf = 0.5*plane.thickness_;
        el.planes.push_back(ip);
        if(auto rect = std::dynamic_pointer_cast<KinKal::Rectangle>(plane.surface_)){
          el.shape = Element::box;
          el.uhalf = rect->uHalfLength();
          el.vhalf = rect->vHalfLength();
          // merge with the next plane if it is the second Gauss plane of the same slab (offset 2t/sqrt(3) along the normal)
          if(ip+1 < planes.size()){
            auto const& next = planes[ip+1];
            auto nrect = std::dynamic_pointer_cast<KinKal::Rectangle>(next.surface_);
            VEC3 const sep = next.surface_->center() - el.center;
            double const gauss = 2.0*plane.thickness_/std::sqrt(3.0);
            if(nrect && next.material_ == plane.material_ && std::fabs(next.thickness_-plane.thickness_) < 1.0e-6
                && std::fabs(std::fabs(next.surface_->normal().Dot(el.wdir))-1.0) < 1.0e-9
                && std::fabs(nrect->uHalfLength()-el.uhalf) < 1.0e-6 && std::fabs(nrect->vHalfLength()-el.vhalf) < 1.0e-6
                && (sep - sep.Dot(el.wdir)*el.wdir).R() < 1.0e-6 && std::fabs(std::fabs(sep.Dot(el.wdir))-gauss) < 1.0e-6){
              el.center += 0.5*sep;
              el.whalf = plane.thickness_;
              el.planes.push_back(++ip);
            }
          }
        } else if(auto ann = std::dynamic_pointer_cast<KinKal::Annulus>(plane.surface_)){
          el.shape = Element::tube;
          el.rin = ann->innerRadius();
          el.rout = ann->outerRadius();
        } else {
          continue;
        }
        elems.push_back(el);
      }
      if(kkg.CRV()){
        for(size_t isect = 0; isect < kkg.CRV()->sectors().size(); ++isect){
          auto const& sector = kkg.CRV()->sectors()[isect];
          Element el;
          el.name = sector.sname_;
          el.sector = int(isect);
          el.kind = Element::crvSector;
          el.shape = Element::box;
          el.center = sector.sector_->center();
          el.udir = sector.sector_->uDirection();
          el.vdir = sector.sector_->vDirection();
          el.wdir = sector.sector_->normal();
          el.uhalf = sector.sector_->uHalfLength();
          el.vhalf = sector.sector_->vHalfLength();
          el.whalf = sector.whw_;
          elems.push_back(el);
        }
      }
      return elems;
    }

    std::vector<Overlap> overlaps(Elements const& elems, unsigned nsample, unsigned seed, double minVolume, bool skipSectorPairs) {
      std::vector<Overlap> found;
      std::vector<VEC3> lo(elems.size()), hi(elems.size());
      for(size_t ie = 0; ie < elems.size(); ++ie) elems[ie].bounds(lo[ie],hi[ie]);
      std::mt19937 rng(seed);
      for(size_t ie = 0; ie < elems.size(); ++ie){
        for(size_t je = ie+1; je < elems.size(); ++je){
          if(skipSectorPairs && elems[ie].kind == Element::crvSector && elems[je].kind == Element::crvSector) continue;
          if(lo[ie].X() >= hi[je].X() || lo[je].X() >= hi[ie].X() || lo[ie].Y() >= hi[je].Y() || lo[je].Y() >= hi[ie].Y()
              || lo[ie].Z() >= hi[je].Z() || lo[je].Z() >= hi[ie].Z()) continue;
          // sample the smaller element
          auto const& a = elems[ie].volume() < elems[je].volume() ? elems[ie] : elems[je];
          auto const& b = &a == &elems[ie] ? elems[je] : elems[ie];
          unsigned nin = 0;
          for(unsigned is = 0; is < nsample; ++is) if(b.inside(a.sample(rng))) ++nin;
          double const vol = a.volume()*nin/double(nsample);
          if(nin > 0 && vol > minVolume) found.push_back(Overlap{ie,je,vol});
        }
      }
      return found;
    }

    std::vector<SectorPair> sectorLayerSums(KinKalGeom const& kkg, double tol, double step) {
      std::vector<SectorPair> found;
      if(!kkg.CRV()) return found;
      auto const& sectors = kkg.CRV()->sectors();
      for(size_t is = 0; is < sectors.size(); ++is){
        auto const& r1 = *sectors[is].sector_;
        for(size_t js = is+1; js < sectors.size(); ++js){
          auto const& r2 = *sectors[js].sector_;
          if(std::fabs(std::fabs(r1.normal().Dot(r2.normal()))-1.0) > 1.0e-9) continue;
          VEC3 const sep = r2.center() - r1.center();
          if(std::fabs(sep.Dot(r1.normal())) > tol) continue;
          // shared area in the axes of r1
          auto range = [&](VEC3 const& axis, double half1, double& lo, double& hi){
            double const half2 = r2.uHalfLength()*std::fabs(r2.uDirection().Dot(axis)) + r2.vHalfLength()*std::fabs(r2.vDirection().Dot(axis));
            lo = std::max(-half1,sep.Dot(axis)-half2); hi = std::min(half1,sep.Dot(axis)+half2);
            return hi > lo;
          };
          double ulo, uhi, vlo, vhi;
          if(!range(r1.uDirection(),r1.uHalfLength(),ulo,uhi) || !range(r1.vDirection(),r1.vHalfLength(),vlo,vhi)) continue;
          SectorPair pair{is,js,2.0,0.0,0};
          // lines along the normal and tilted by 45 degrees either way across the bars, every 'step' across the shared band
          for(double slope : {0.0,1.0,-1.0}){
            VEC3 const dir = (r1.normal() + slope*r1.vDirection()).Unit();
            for(double v = vlo + 0.5*step; v < vhi; v += step){
              VEC3 const pos = r1.center() + 0.5*(ulo+uhi)*r1.uDirection() + v*r1.vDirection();
              double const sum = sectors[is].layerFraction(pos,dir) + sectors[js].layerFraction(pos,dir);
              pair.minsum = std::min(pair.minsum,sum);
              pair.maxsum = std::max(pair.maxsum,sum);
              ++pair.nlines;
            }
          }
          if(pair.nlines > 0) found.push_back(pair);
        }
      }
      return found;
    }

    std::vector<Coplanar> coplanarPlanes(KinKalGeom const& kkg, double tol) {
      std::vector<Coplanar> found;
      auto const& planes = kkg.passiveMaterialPlanes();
      for(size_t ip = 0; ip < planes.size(); ++ip){
        auto r1 = std::dynamic_pointer_cast<KinKal::Rectangle>(planes[ip].surface_);
        if(!r1) continue;
        for(size_t jp = ip+1; jp < planes.size(); ++jp){
          auto r2 = std::dynamic_pointer_cast<KinKal::Rectangle>(planes[jp].surface_);
          if(!r2 || planes[jp].material_ != planes[ip].material_ || std::fabs(planes[jp].thickness_-planes[ip].thickness_) > 1.0e-6) continue;
          if(std::fabs(std::fabs(r1->normal().Dot(r2->normal()))-1.0) > 1.0e-9) continue;
          VEC3 const sep = r2->center() - r1->center();
          if(std::fabs(sep.Dot(r1->normal())) > tol) continue;
          // in-plane extents of r2 along the axes of r1
          auto overlap1d = [&](VEC3 const& axis, double half1){
            double const half2 = r2->uHalfLength()*std::fabs(r2->uDirection().Dot(axis)) + r2->vHalfLength()*std::fabs(r2->vDirection().Dot(axis));
            double const c2 = sep.Dot(axis);
            return std::min(half1,c2+half2) - std::max(-half1,c2-half2);
          };
          double const uo = overlap1d(r1->uDirection(),r1->uHalfLength());
          double const vo = overlap1d(r1->vDirection(),r1->vHalfLength());
          if(uo > -tol && vo > -tol) found.push_back(Coplanar{ip,jp,uo,vo});
        }
      }
      return found;
    }
  }
}
