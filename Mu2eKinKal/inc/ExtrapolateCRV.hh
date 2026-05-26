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
      using RecPtr = std::shared_ptr<KinKal::Rectangle>;
      ExtrapolateCRV() : maxDt_(-1.0), dptol_(1e10), intertol_(1e10), minvnorm_(1e-5),
      debug_(0){}

      ExtrapolateCRV(double maxdt, double dptol, double intertol, double minv, CRV const& crv, int debug=0) :
        maxDt_(maxdt), dptol_(dptol), intertol_(intertol), minvnorm_(minv), debug_(debug) {
          for(auto const& sector : crv.sectors() ) {
            sectors_.push_back(sector.sector_);
          }
        }
      // interface for extrapolation
      double maxDt() const { return maxDt_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      // CRV specific
      auto const& sectors() const { return sectors_; }
      auto const& intersection() const { return inter_; }
      auto const& sector () const { return sect_; }
      int sectorIndex() const { return isect_; }
      int debug() const { return debug_; }
      void reset() const { inter_ = Intersection(); sect_ = RecPtr(); isect_ = -1; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
    private:
      double maxDt_; // maximum extrapolation time
      double dptol_; // fractional momentum tolerance
      double intertol_; // intersection tolerance (mm)
      double minvnorm_; // minimum vel normal (outwards) to plane
      mutable Intersection inter_; // cache of most recent intersection
      mutable RecPtr sect_; // cache of most recent intersection
      mutable int isect_; // index to intersected sector
      int debug_; // debug level
      std::vector<RecPtr> sectors_; // rectangles at the (layer) middle of each sector
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
    if(debug_ > 4)std::cout << "CRV extrap tdir " << tdir << " stime " << stime << " etime " << etime << std::endl;
    // test from the start of the end poiece
    auto pos = ktraj.position3(stime);
    auto vel = ktraj.velocity(stime);
    if(debug_ > 4)std::cout << "CRV extrap testing time " << time << " vel " << vel << " pos " << pos << std::endl;
    for(size_t isect = 0; isect < sectors_.size(); ++isect ){
      auto const& sector = sectors_[isect];
      double normvel = vel.Dot(sector->normal())*timeDirSign(tdir); // sign by extrapolation direction
      double sdist = (sector->center()-pos).Dot(vel)*timeDirSign(tdir);
      if(debug_ > 4)std::cout << "CRV extrap normvel " << normvel << " time " << time << " sdist " << sdist << std::endl;
      // stop if horizontal or plane is past the current point.
      if(fabs(normvel) < minvnorm_ || sdist < 0 )continue;
      // try to intersect
      auto trange = tdir == TimeDir::forwards ? TimeRange(stime,etime) : TimeRange(etime,stime);
      auto newinter = KinKal::intersect(fittraj,*sector,trange,intertol_,tdir);
      if(debug_ > 3)std::cout << "CRV " << newinter  << std::endl;
      if(newinter.good()){
        inter_ = newinter;
        sect_ = sector;
        isect_ = isect;
        if(debug_ > 1)std::cout << "Good CRV " <<  newinter << std::endl;
        break;
      } else if ( newinter.onsurface_ && newinter.inbounds_) { // inbounds might be too strict for CentralHelix tracks, will need to check TODO
        // there's a potential intersection, but the trajectory hasn't gotten there yet. Tell the track to keep extending
        retval = true;
        if(debug_ > 2)std::cout << "Potential CRV " <<  newinter << std::endl;
      }
    }
    return retval;
  }
}
#endif
