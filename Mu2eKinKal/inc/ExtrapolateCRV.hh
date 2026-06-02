// predicate to extrapolate to CRV
#ifndef Mu2eKinKal_ExtrapolateCRV_hh
#define Mu2eKinKal_ExtrapolateCRV_hh
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/Rectangle.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "Offline/KinKalGeom/inc/CRV.hh"
#include "cetlib_except/exception.h"
#include <memory>
#include <vector>
#include <limits>
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::TimeRange;
  using KKGeom::CRV;
  using KinKal::Rectangle;
  using KinKal::Intersection;
  class ExtrapolateCRV {
    public:
      using CRVSV = std::vector<KKGeom::KKCRVSector>;
      struct CRVI {
        Intersection inter_;
        double whw_ = 0;
        int isect_ = -1;
        CRVI(){}
        CRVI(Intersection const& inter,double whw, int isect) : inter_(inter),whw_(whw),isect_(isect) {}
      };
      using CRVIV = std::vector<CRVI>;

      struct sortIntersections { // sort in time direction
        TimeDir tdir_;
        bool operator () (CRVI const& inter1, CRVI const& inter2) {
          return tdir_ == TimeDir::forwards ? inter1.inter_.time_ < inter2.inter_.time_ : inter1.inter_.time_ > inter2.inter_.time_;
        }
        sortIntersections(TimeDir tdir) : tdir_(tdir) {}
      };

      ExtrapolateCRV(double maxdt, double maxdtstep, double dptol, double intertol, double minv, CRV const& crv, int debug=0) :
        maxDt_(maxdt), maxDtStep_(maxdtstep), dptol_(dptol), intertol_(intertol), minvnorm_(minv), debug_(debug), sectors_(crv.sectors()) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; }
      double maxDtStep() const { return maxDtStep_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      // CRV specific
      auto const& sectors() const { return sectors_; }
      auto const& sector(size_t isect) const { return sectors_[isect]; }
      auto const& intersections() const { return inters_; }
      int debug() const { return debug_; }
      void reset() const { inters_.clear(); }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
    private:
      double maxDt_ = -1; // maximum extrapolation time
      double maxDtStep_ = -1; // maximum extrapolation time step in a single iteration
      double dptol_ = 1e10; // fractional momentum tolerance
      double intertol_ = 1e10; // intersection tolerance (mm)
      double minvnorm_ = 1e-5; // minimum vel normal (outwards) to plane
      int debug_ = 0; // debug level
      CRVSV sectors_; // sectors cache
      mutable CRVIV inters_; // cache of most recent intersections
  };

  template <class KTRAJ> bool ExtrapolateCRV::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    // we are answering the question: did the segment last added to this extrapolated track cross a CRV sector or not?
    // if so, stop extrapolating (for now). If not, and if we're still heading towards at least 1 CRV sector , keep going.
    // cache the intersection if it's found, so it can be used without recomputing it.
    bool retval(false);
    reset();
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    // Make sure the range is positive definite
    static const double epsilon(1.0e-6);
    auto stime = tdir == TimeDir::forwards ? ktraj.range().begin()+epsilon : ktraj.range().end()-epsilon;
    auto etime = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    TimeRange trange(stime,etime,false);
    // test from the start of the end piece
    auto pos = ktraj.position3(stime);
    auto vel = ktraj.velocity(stime);
    if(debug_ > 4)std::cout << "CRV extrap tdir " << tdir << " stime " << stime << " etime " << etime << " vel " << vel << " pos " << pos << std::endl;
    for(size_t isect = 0; isect < sectors_.size(); ++isect ){
      auto const& sector = sectors_[isect];
      double normvel = vel.Dot(sector.sector_->normal())*timeDirSign(tdir); // sign by extrapolation direction
      double sdist = (sector.sector_->center()-pos).Dot(vel)*timeDirSign(tdir);
      if(debug_ > 4)std::cout << "CRV extrap normvel " << normvel << " sdist " << sdist << std::endl;
      // stop if horizontal or plane is past the current point.
      if(fabs(normvel) < minvnorm_ || sdist < 0 )continue;
      // try to intersect
      auto newinter = KinKal::intersect(fittraj,*sector.sector_,trange,intertol_,tdir);
      if(debug_ > 3)std::cout << "CRV " << newinter  << std::endl;
      if(newinter.good()){
        inters_.emplace_back(newinter,sector.whw_,(int)isect);
        if(debug_ > 1)std::cout << "Good CRV " <<  newinter << " sector " << sector.sname_ << std::endl;
      } else if ( newinter.onsurface_ && newinter.inbounds_) { // inbounds might be too strict for CentralHelix tracks, will need to check TODO
        retval |= trange.beyond(newinter.time_,tdir);
        // there's a potential intersection, but the trajectory hasn't gotten there yet. Tell the track to keep extending
        if(trange.beyond(newinter.time_,tdir) && debug_ > 2)std::cout << "Potential CRV " <<  newinter << std::endl;
      }
    }
    // sort intersections in the time direction
    sortIntersections isort(tdir);
    std::sort(inters_.begin(),inters_.end(),isort);
    return retval;
  }
}
#endif
