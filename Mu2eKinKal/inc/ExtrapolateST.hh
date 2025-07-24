// predicate to extrapolate to foils in the ST
// track changes direction.
#ifndef Mu2eKinKal_ExtrapolateST_hh
#define Mu2eKinKal_ExtrapolateST_hh
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "Offline/KinKalGeom/inc/StoppingTarget.hh"
#include "cetlib_except/exception.h"
#include <memory>
#include <vector>
#include <limits>
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::TimeRange;
  using KinKalGeom::StoppingTarget;
  using KinKal::Annulus;
  using KinKal::Intersection;
  using AnnPtr = std::shared_ptr<KinKal::Annulus>;
  using FoilCol = std::vector<AnnPtr>;
  using CylPtr = std::shared_ptr<KinKal::Cylinder>;
  class ExtrapolateST {
    public:
      ExtrapolateST() : maxDt_(-1.0), dptol_(1e10), intertol_(1e10),
      zmin_(std::numeric_limits<float>::max()),
      zmax_(-std::numeric_limits<float>::max()),
      rmin_(std::numeric_limits<float>::max()),
      rmax_(-std::numeric_limits<float>::max()),
      debug_(0){}

      ExtrapolateST(double maxdt, double dptol, double intertol, StoppingTarget const& stoptarg, int debug=0) :
        maxDt_(maxdt), dptol_(dptol), intertol_(intertol),
        zmin_( (stoptarg.outer().center() - stoptarg.outer().axis()*stoptarg.outer().halfLength()).Z()),
        zmax_( (stoptarg.outer().center() + stoptarg.outer().axis()*stoptarg.outer().halfLength()).Z()),
        rmin_( stoptarg.inner().radius()), rmax_(stoptarg.outer().radius()),
        foils_(stoptarg.foils()),debug_(debug) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      double zmin() const { return zmin_; }
      double zmax() const { return zmax_; }
      double rmin() const { return rmin_; }
      double rmax() const { return rmax_; }
      auto const& foils () const { return foils_; }
      auto const& intersection() const { return inter_; }
      auto const& foil() const { return ann_; }
      auto const& foilId() const { return sid_; }
      int debug() const { return debug_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
      // reset between tracks
      void reset() const { inter_ = Intersection(); sid_ = SurfaceId(); ann_ = AnnPtr();}
      // find the nearest foil to a z positionin a given z direction
      size_t nearestFoil(double zpos, double zdir) const;
    private:
      double maxDt_; // maximum extrapolation time
      double dptol_; // fractional momentum tolerance
      double intertol_; // intersection tolerance (mm)
      mutable Intersection inter_; // cache of most recent intersection
      mutable SurfaceId sid_; // cache of most recent intersection foil SID
      mutable AnnPtr ann_; // cache of most recent intersection foil surface
      // cache of front and back Z positions
      double zmin_, zmax_; // z range of ST volume
      double rmin_, rmax_; // inner and outer radii of the anuli
      FoilCol foils_; // foils
      int debug_; // debug level
  };

  template <class KTRAJ> bool ExtrapolateST::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    // we are answering the question: did the segment last added to this extrapolated track hit a foil or not? If so we are done
    // extrapolating (for now) and we want to find all the intersections in that piece. If not, and if we're still inside or heading towards the
    // ST, keep going.
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    // add a small buffer to the test range to prevent re-intersection with the same piece
    static const double epsilon(1e-7); // small difference to avoid re-intersecting
    if(ktraj.range().range() <= epsilon) return true; // keep going if the step is very small
    auto stime = tdir == TimeDir::forwards ? ktraj.range().begin()+epsilon : ktraj.range().end()-epsilon;
    auto etime = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    auto vel = ktraj.velocity(stime); // physical velocity
    if(tdir == TimeDir::backwards) vel *= -1.0;
    auto spos = ktraj.position3(stime);
    auto epos = ktraj.position3(etime);
    if(debug_ > 2)std::cout << "ST extrap tdir " << tdir << " start z " << spos.Z() << " end z " << epos.Z() << " zvel " << vel.Z() << " rho " << spos.Rho() << std::endl;
    // stop if the particle is heading away from the ST
    if( (vel.Z() > 0 && spos.Z() > zmax_ ) || (vel.Z() < 0 && spos.Z() < zmin_)){
      reset(); // clear any cache
      if(debug_ > 1)std::cout << "Heading away from ST: done" << std::endl;
      return false;
    }
    // if the particle is going in the right direction but haven't yet reached the ST in Z just keep going
    if( (vel.Z() > 0 && epos.Z() < zmin_) || (vel.Z() < 0 && epos.Z() > zmax_) ){
      reset();
      if(debug_ > 2)std::cout << "Heading towards ST, z " << spos.Z()<< std::endl;
      return true;
    }
    // if we get to here we are in the correct Z range. Test foils.
    int ifoil = nearestFoil(spos.Z(),vel.Z());
    if(ifoil >= (int)foils_.size())return true;
    if(debug_ > 2)std::cout << "Looping on foils " << std::endl;
    int dfoil = vel.Z() > 0.0 ? 1 : -1; // iteration direction
    auto trange = tdir == TimeDir::forwards ? TimeRange(stime,ktraj.range().end()) : TimeRange(ktraj.range().begin(),stime);
    // loop over foils in the z range of this piece
    while(ifoil >= 0 && ifoil < (int)foils_.size() && (foils_[ifoil]->center().Z() - epos.Z())*dfoil < 0.0){
      auto foilptr = foils_[ifoil];
      if(debug_ > 2)std::cout << "foil " << ifoil << " z " << foilptr->center().Z() << std::endl;
      auto newinter = KinKal::intersect(ktraj,*foilptr,trange,intertol_,tdir);
      if(debug_ > 2)std::cout << "ST foil inter " << newinter  << std::endl;
      if(newinter.good()){
        // update the cache
        inter_ = newinter;
        ann_ = foils_[ifoil];
        sid_ = SurfaceId(SurfaceIdEnum::ST_Foils,ifoil);
        if(debug_ > 0)std::cout << "Good ST foil " << newinter << " sid " << sid_ << std::endl;
        return false;
      }
      ifoil += dfoil; // otherwise continue loopin on foils
    }
    // no more intersections: keep extending in Z till we clear the ST
    reset();
    if(debug_ > 1)std::cout << "Extrapolating to ST edge, z " << spos.Z() << std::endl;
    if(vel.Z() > 0.0)
      return spos.Z() < zmax_;
    else
      return spos.Z() > zmin_;
  }

  size_t ExtrapolateST::nearestFoil(double zpos, double zvel) const {
    size_t retval = foils_.size();
    if(zvel > 0.0){ // going forwards in z
      for(auto ifoil= foils_.begin(); ifoil != foils_.end(); ifoil++){
        auto const& foilptr = *ifoil;
        if(foilptr->center().Z() > zpos){
          retval = std::distance(foils_.begin(),ifoil);
          break;
        }
      }
    } else {
      for(auto ifoil= foils_.rbegin(); ifoil != foils_.rend(); ifoil++){
        auto const& foilptr = *ifoil;
        if(foilptr->center().Z() < zpos){
          auto jfoil = ifoil.base()-1; // points to the equivalent forwards object
          retval = std::distance(foils_.begin(),jfoil);
          break;
        }
      }
    }
    return retval;
  }

}
#endif
