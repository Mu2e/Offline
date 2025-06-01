// predicate to extrapolate to TestCRV
#ifndef Mu2eKinKal_ExtrapolateTCRV_hh
#define Mu2eKinKal_ExtrapolateTCRV_hh
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "Offline/KinKalGeom/inc/TestCRV.hh"
#include "cetlib_except/exception.h"
#include <memory>
#include <vector>
#include <limits>
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::TimeRange;
  using KinKalGeom::TestCRV;
  using KinKal::Rectangle;
  using KinKal::Intersection;
  using RecPtr = std::shared_ptr<KinKal::Rectangle>;
  class ExtrapolateTCRV {
    public:
      ExtrapolateTCRV() : maxDt_(-1.0), dptol_(1e10), intertol_(1e10), minv_(1e-5),
      step_(0),
      ymin_(std::numeric_limits<float>::max()),
      ymax_(-std::numeric_limits<float>::max()),
      debug_(0){}

      ExtrapolateTCRV(double maxdt, double dptol, double intertol, double minv, TestCRV const& tcrv, int debug=0) :
        maxDt_(maxdt), dptol_(dptol), intertol_(intertol), minv_(minv), step_(0), debug_(debug)
    {
      modules_.push_back(tcrv.ex1Ptr());
      modules_.push_back(tcrv.t1Ptr());
      modules_.push_back(tcrv.t2Ptr());
      ymin_ = std::min(std::min(tcrv.ex1Ptr()->center().Y(),tcrv.t1Ptr()->center().Y()),tcrv.t2Ptr()->center().Y());
      ymax_ = std::max(std::max(tcrv.ex1Ptr()->center().Y(),tcrv.t1Ptr()->center().Y()),tcrv.t2Ptr()->center().Y());
    }

      // interface for extrapolation
      double maxDt() const { return maxDt_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      double step() const { return step_; }
      double ymin() const { return ymin_; }
      double ymax() const { return ymax_; }
      auto const& modules () const { return modules_; }
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
      // cache of front and back Z positions
      double ymin_, ymax_; // z range of ST volume
      std::vector<RecPtr> modules_;
      int debug_; // debug level
  };

  template <class KTRAJ> bool ExtrapolateTCRV::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    auto time = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    auto vel = ktraj.velocity(time);
    auto pos = ktraj.position3(time);
    double yvel = vel.Y()*timeDirSign(tdir); // sign by extrapolation direction
    double yval = pos.Y();
    if(debug_ > 2)std::cout << "TCRV extrap tdir " << tdir << " start y " << pos.Y() << " yvel " << yvel << " time " << time << std::endl;
    // stop if horizontal
    if(fabs(yvel) < minv_){
      reset(); // clear any cache
      return false;
    }

    static const double epsilon(1e-6); // small step to avoid re-intersecting

    // if the particle is going in the right direction but hasn't yet reached the CRV yet keep going
    if( (yvel > 0 && pos.Y() < ymin_) || (yvel < 0 && pos.Y() > ymax_) ){
      if (yvel > 0){
        step_ = std::min(maxDt_,fabs((ymax_ - pos.Y())/yvel)+epsilon);
      }else{
        step_ = std::min(maxDt_,fabs((pos.Y()-ymin_)/yvel)+epsilon);
      }
      reset(); // clear any cache
      return true;
    }

    // if we get to here we need to test for an intersection with the actual cylinder. Make sure the range is positive definite
    auto stime = tdir == TimeDir::forwards ? ktraj.range().begin()+epsilon : ktraj.range().end()-epsilon;
    auto etime = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    auto valid = tdir == TimeDir::forwards ? stime < etime : etime < stime;
    if (valid){
      auto trange = tdir == TimeDir::forwards ? TimeRange(stime,etime) : TimeRange(etime,stime);
      for (size_t i=0;i<modules_.size();i++){
        auto newinter = KinKal::intersect(fittraj,*modules_[i],trange,intertol_,tdir);
        if(debug_ > 2)std::cout << "TCRV " << newinter  << std::endl;
        if(newinter.good()){
          inter_ = newinter;
          rec_ = modules_[i];
          if(debug_ > 0)std::cout << "Good TCRV " <<  newinter << std::endl;
          return false;
        }
      }
    }
    reset();

    // if no intersections left
    // stop if we're heading away from the target z
    if((yvel > 0 && yval > ymax_ ) || (yvel < 0 && yval < ymin_))return false;

    // no intersections yet: keep extending in Y till we clear the CRV
    if(yvel > 0.0){
      step_ = std::min(maxDt_,fabs((ymax_ - pos.Y())/yvel)+epsilon);
      return pos.Y() < ymax_;
    }
    else{
      step_ = std::min(maxDt_,fabs((pos.Y()-ymin_)/yvel)+epsilon);
      return pos.Y() > ymin_;
    }
  }
}
#endif
