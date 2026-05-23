// predicate to extrapolate to CRV
#ifndef Mu2eKinKal_ExtrapolateCRV_hh
#define Mu2eKinKal_ExtrapolateCRV_hh
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
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
  using RecPtr = std::shared_ptr<KinKal::Rectangle>;
  class ExtrapolateCRV {
    public:
      ExtrapolateCRV() : maxDt_(-1.0), dptol_(1e10), intertol_(1e10), minv_(1e-5),
      step_(0),
      debug_(0){}

      ExtrapolateCRV(double maxdt, double dptol, double intertol, double minv, CRV const& crv, int debug=0) :
        maxDt_(maxdt), dptol_(dptol), intertol_(intertol), minv_(minv), step_(0), debug_(debug),sectors_(crv.sectors())
    {
      // interface for extrapolation
      double maxDt() const { return maxDt_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      double step() const { return step_; }
      auto const& modules () const { return sectors_; }
      auto const& intersection() const { return inter_; }
      auto const& module () const { return rec_; }
      int debug() const { return debug_; }
      void reset() const { inter_ = Intersection(); rec_ = RecPtr();}
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
    private:
      double maxDt_; // maximum extrapolation time
      double dptol_; // fractional momentum tolerance
      double intertol_; // intersection tolerance (mm)
      double minv_; // minimum vel.Y
      mutable double step_; // predicted step length
      mutable Intersection inter_; // cache of most recent intersection
      mutable RecPtr rec_; // cache of most recent intersection foil surface
      std::vector<RecPtr> sectors_;
      int debug_; // debug level
  };

  template <class KTRAJ> bool ExtrapolateCRV::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    bool retval(false);
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    // Make sure the range is positive definite
    auto stime = tdir == TimeDir::forwards ? ktraj.range().begin()+epsilon : ktraj.range().end()-epsilon;
    auto etime = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    auto valid = tdir == TimeDir::forwards ? stime < etime : etime < stime;
    if (valid){
    auto time = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    auto vel = ktraj.velocity(time);
    for(auto const& sector : sectors_){
      double normvel = vel.Dot(sector->normal())*timeDirSign(tdir); // sign by extrapolation direction
    if(debug_ > 2)std::cout << "CRV extrap tdir " << tdir << " normvel " << normvel << " time " << time << std::endl;
    // stop if horizontal
    if(fabs(normvel) < minv_)continue;

      auto trange = tdir == TimeDir::forwards ? TimeRange(stime,etime) : TimeRange(etime,stime);
      for (size_t i=0;i<sectors_.size();i++){
        auto newinter = KinKal::intersect(fittraj,*sectors_[i],trange,intertol_,tdir);
        if(debug_ > 2)std::cout << "CRV " << newinter  << std::endl;
        if(newinter.good()){
          inter_ = newinter;
          rec_ = sectors_[i];
          if(debug_ > 0)std::cout << "Good CRV " <<  newinter << std::endl;
          return false;
        }
      }
    }
    reset();

}
#endif
