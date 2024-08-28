// predicate to extrapolate to the next intersection with the IPA
// track changes direction.
#ifndef Mu2eKinKal_ExtrapolateIPA_hh
#define Mu2eKinKal_ExtrapolateIPA_hh
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/General/TimeDir.hh"
#include <limits>
#include "cetlib_except/exception.h"
namespace mu2e {
  class ExtrapolateIPA {
    public:
      using CylPtr = std::shared_ptr<KinKal::Cylinder>;
      using KinKal::TimeRange;
      ExtrapolateIPA() : maxDt_(-1.0), tol_(1e10) {}
      ExtrapolateIPA(double maxdt, double tol, double zmin, double zmax) :
        maxDt_(maxdt), tol_(tol), zmin_(zmin), zmax_(zmax) {
          if(zmin >= zmax) throw cet::exception("RECO")<<"Mu2eKinKal::ExtrapolateIPA: range configuration error, zmin "<< zmin << " zmax " << zmax <<endl;
        }
      // interface for extrapolation
      double maxDt() const { return maxDt_; } // maximum time to extend the track
      double tolerance() const { return tol_; } // intersection tolerance
      CylPtr const& IPACylinder() const { return ipa_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::Track<KTRAJ> const& ktrk, KinKal::TimeDir tdir, double time) const;
    private:
      double maxDt_; // maximum extrapolation time
      double tol_; // intersection tolerance (mm)
      CylPtr ipa_; // IPA cylinder
      Intersection inter_; // most recent intersection
  };

  template <class KTRAJ> bool ExtrapolateIPA::needsExtrapolation(KinKal::Track<KTRAJ> const& kktrk, KinKal::TimeDir tdir, double time) const {
    // sanity check on initial position
    double dz = kktrk.fitTraj().position3(time).Z() - ipa_->center().Z();
    double halflen = ipa_->halfLength();
    if(dz < -halflen || dz > halflen) return false; // need to qualify this in the case of reflection TODO
    // estimate time range based on local piece Z speed transit time
    double tz = std::min(2*halflen/std::max(fabs(kktrk.fitTraj().velocity(time).Z()),1.0),maxDT); // protect against reflection (zero z speed)
    TimeRange trange = tdir == TimeDirection::forwards TimeRange(time,time+tz) : TimeRange(time-tz,time)
    // update intersection
    inter_ = KinKal::intersect(kktrk.fitTraj(),*ipa_,trange,tol_,tdir);
    return inter_.onsurface_ && inter_.inbounds_ && inter_.time_ < time;
  }
}
#endif
