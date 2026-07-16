// Common scaffolding shared by the surface-crossing extrapolation predicates (ExtrapolateCylinders,
// ExtrapolatePlanes): dt/tolerance bookkeeping, intersection storage, and time-ordering. Each
// predicate's needsExtrapolation search is geometry-specific and NOT shared here: ExtrapolateCylinders
// marches in bounded sub-windows to avoid stepping over closely-spaced concentric shells, while
// ExtrapolatePlanes gates directly on normal velocity/distance over the whole piece. Forcing those two
// search strategies into one shared body would obscure both, so only the bookkeeping below is common.
#ifndef Mu2eKinKal_ExtrapolateSurface_hh
#define Mu2eKinKal_ExtrapolateSurface_hh

#include "KinKal/General/TimeDir.hh"

#include <vector>

namespace mu2e {
  using KinKal::TimeDir;

  template <class SurfaceIntersection> class ExtrapolateSurface {
    public:
      using SurfaceIntersectionCollection = std::vector<SurfaceIntersection>;

      struct SortIntersections {
        TimeDir tdir_;
        bool operator () (SurfaceIntersection const& lhs, SurfaceIntersection const& rhs) {
          return tdir_ == TimeDir::forwards ? lhs.inter_.time_ < rhs.inter_.time_ : lhs.inter_.time_ > rhs.inter_.time_;
        }
        SortIntersections(TimeDir tdir) : tdir_(tdir) {}
      };

      double maxDt() const { return maxDt_; }
      double maxDtStep() const { return maxDtStep_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      auto const& intersections() const { return inters_; }
      int debug() const { return debug_; }
      void reset() const { inters_.clear(); }

    protected:
      ExtrapolateSurface(double maxdt, double maxdtstep, double dptol, double intertol, double minv, int debug) :
        maxDt_(maxdt), maxDtStep_(maxdtstep), dptol_(dptol), intertol_(intertol), minvnorm_(minv), debug_(debug) {}

      double maxDt_ = -1;
      double maxDtStep_ = -1;
      double dptol_ = 1e10;
      double intertol_ = 1e10;
      double minvnorm_ = 1e-5;
      int debug_ = 0;
      mutable SurfaceIntersectionCollection inters_;
  };
}

#endif
