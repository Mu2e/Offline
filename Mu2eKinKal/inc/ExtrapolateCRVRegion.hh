// Predicate to extrapolate through the CRV region, crossing BOTH the CRV sector planes AND the
// per-sector aluminium strongback (passive) planes in ONE interleaved, time-ordered pass.
//
// WHY THIS EXISTS
// Each CRV module has a thin (~12.7 mm) Al strongback on the tracker-facing side of its scintillator
// stack. The sector (veto) plane is modeled at the stack center, so the strongback plane sits ~60 mm
// (about a stack half-thickness) toward the tracker from it; going outward the surfaces alternate
// strongback -> sector for every module. KinKal xings can only be ADDED at the trajectory front/back (an interior xing cannot
// be inserted), so these surfaces must be crossed in strict spatial/time order. The previous code
// crossed them in two SEPARATE passes -- all strongbacks (ExtrapolatePlanes), then all sectors
// (ExtrapolateCRV). In a STACKED-sector geometry (the extracted KPP modules EX/T1/T2 sit ~150 mm
// apart in y) the strongback pass advanced the trajectory front PAST the lower sectors, leaving
// them stranded in the interior where they could no longer be added -- CRV_T1 was lost on every
// track, and EX/T2 were under-counted. Searching both surface types together, this predicate stops
// the extrapolation at the NEXT crossing of EITHER type, so the KKExtrap loop adds them strictly
// outward (T1_SB -> T1 -> EX_SB -> EX -> ...) and nothing is stranded.
//
// RUN-2 SAFETY: run-2 geometry presents at most one CRV sector (and its one strongback) per
// extrapolation direction, so the outward order is unchanged (strongback then sector) and no
// stacking occurs. Concrete passive planes are NOT handled here (they stay in
// extrapolatePassiveMaterial); only the CRV strongbacks are merged with the sectors. Sectors and
// strongback planes share the same minvnorm floor (a near-zero guard against the normvel divide in
// time_to_sector, not a physics cut): grazing incidence is not pre-filtered by a velocity threshold,
// it is left to fail at the actual surface intersection (newinter.good()), same as any other surface.
#ifndef Mu2eKinKal_ExtrapolateCRVRegion_hh
#define Mu2eKinKal_ExtrapolateCRVRegion_hh
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/Rectangle.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "Offline/KinKalGeom/inc/CRV.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::TimeRange;
  using KKGeom::CRV;
  using KinKal::Intersection;
  using KinKal::timeDirSign;
  class ExtrapolateCRVRegion {
    public:
      using CRVSV = std::vector<KKGeom::KKCRVSector>;
      using MaterialPlane = KinKalGeom::PassiveMaterialPlane;
      using MaterialPlaneCollection = std::vector<MaterialPlane>;
      // a crossing of either a CRV sector or a strongback passive plane
      struct RegionXing {
        Intersection inter_;
        bool isSector_ = false;
        int sectorIdx_ = -1;                    // index into sectors_ (sector crossings)
        double whw_ = 0;                        // sector module half-width (sector crossings)
        MaterialPlane const* plane_ = nullptr;  // the strongback plane (plane crossings)
        RegionXing() {}
        RegionXing(Intersection const& inter, int isect, double whw) :
          inter_(inter), isSector_(true), sectorIdx_(isect), whw_(whw) {}
        RegionXing(Intersection const& inter, MaterialPlane const& plane) :
          inter_(inter), isSector_(false), plane_(&plane) {}
      };
      using RegionXingCollection = std::vector<RegionXing>;
      struct sortXings { // sort in the extrapolation time direction
        TimeDir tdir_;
        bool operator () (RegionXing const& a, RegionXing const& b) {
          return tdir_ == TimeDir::forwards ? a.inter_.time_ < b.inter_.time_ : a.inter_.time_ > b.inter_.time_;
        }
        sortXings(TimeDir tdir) : tdir_(tdir) {}
      };

      ExtrapolateCRVRegion(double maxdt, double maxdtstep, double dptol, double intertol,
          double minv, CRV const& crv, MaterialPlaneCollection const& planes, int debug=0) :
        maxDt_(maxdt), maxDtStep_(maxdtstep), dptol_(dptol), intertol_(intertol),
        minvnorm_(minv), debug_(debug), sectors_(crv.sectors()), planes_(planes) {}
      double maxDt() const { return maxDt_; }
      double maxDtStep() const { return maxDtStep_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      auto const& sectors() const { return sectors_; }
      auto const& sector(size_t isect) const { return sectors_[isect]; }
      auto const& intersections() const { return inters_; }
      int debug() const { return debug_; }
      void reset() const { inters_.clear(); }
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
    private:
      double maxDt_ = -1;
      double maxDtStep_ = -1;
      double dptol_ = 1e10;
      double intertol_ = 1e10;
      double minvnorm_ = 1e-5;
      int debug_ = 0;
      CRVSV sectors_;
      MaterialPlaneCollection const& planes_;
      mutable RegionXingCollection inters_;
  };

  template <class KTRAJ> bool ExtrapolateCRVRegion::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    // Stop the extrapolation at the NEXT crossing of any CRV sector OR strongback plane. The
    // KKExtrap loop adds whichever crossing(s) this piece produced (time-sorted, outermost first)
    // then re-extrapolates, so the surfaces are crossed strictly outward and none is left interior.
    bool retval(false);
    reset();
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    static const double epsilon(1.0e-6);
    static const double coincidenttol(1.0e-3); // ns; crossings closer than this are the same physical point (a seam)
    auto stime = tdir == TimeDir::forwards ? ktraj.range().begin()+epsilon : ktraj.range().end()-epsilon;
    auto etime = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    TimeRange trange(stime,etime,false);
    auto pos = ktraj.position3(stime);
    auto vel = ktraj.velocity(stime);
    double piece_span = fabs(stime - etime);
    // helper: has a crossing at this time already been cached (a coplanar seam)?
    auto coincident = [&](double t){ for(auto const& p : inters_) if(fabs(p.inter_.time_ - t) < coincidenttol) return true; return false; };
    // --- CRV sectors ---
    for(size_t isect = 0; isect < sectors_.size(); ++isect ){
      auto const& sector = sectors_[isect];
      double normvel = vel.Dot(sector.sector_->normal())*timeDirSign(tdir);
      double signed_perp = (sector.sector_->center()-pos).Dot(sector.sector_->normal());
      double time_to_sector = signed_perp / normvel;
      // skip planes already behind this piece's evaluation point (already crossed/added); grazing
      // incidence is not pre-filtered here, it is left to fail at the intersection attempt below
      if(fabs(normvel) < minvnorm_ || time_to_sector < 0 )continue;
      auto newinter = KinKal::intersect(fittraj,*sector.sector_,trange,intertol_,tdir);
      if(newinter.good()){
        if(!coincident(newinter.time_)) inters_.emplace_back(newinter,(int)isect,sector.whw_);
      } else if(newinter.onsurface_ && newinter.inbounds_) {
        retval |= trange.beyond(newinter.time_,tdir); // crossing just beyond this piece: keep extrapolating
      } else if(!newinter.onsurface_ && tdir == TimeDir::backwards) {
        if(time_to_sector > piece_span) retval = true; // far sector beyond this (short) piece: keep going
      }
    }
    // --- strongback passive planes (minvnorm floor, like ExtrapolatePlanes) ---
    for(auto const& plane : planes_){
      double normvel = vel.Dot(plane.surface_->normal())*timeDirSign(tdir);
      double signed_perp = (plane.surface_->center()-pos).Dot(plane.surface_->normal());
      double time_to_plane = signed_perp / normvel;
      if(fabs(normvel) < minvnorm_ || time_to_plane < 0 )continue;
      auto newinter = KinKal::intersect(fittraj,*plane.surface_,trange,intertol_,tdir);
      if(newinter.good()){
        if(!coincident(newinter.time_)) inters_.emplace_back(newinter,plane);
      } else if(newinter.onsurface_ && newinter.inbounds_) {
        retval |= trange.beyond(newinter.time_,tdir);
      } else if(!newinter.onsurface_ && tdir == TimeDir::backwards) {
        if(time_to_plane > piece_span) retval = true;
      }
    }
    // sort the crossings in the time direction so the KKExtrap loop adds them outermost-first
    std::sort(inters_.begin(),inters_.end(),sortXings(tdir));
    if(!inters_.empty()) return false; // stop: this piece crossed the CRV region; record now
    return retval;
  }
}
#endif
