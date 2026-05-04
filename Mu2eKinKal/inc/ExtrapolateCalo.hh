// predicate to extrapolate to calo
// Sophie Middleton (2025)
#ifndef Mu2eKinKal_ExtrapolateCalo_hh
#define Mu2eKinKal_ExtrapolateCalo_hh
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "Offline/KinKalGeom/inc/Calo.hh"
#include "cetlib_except/exception.h"
#include <memory>
#include <vector>
#include <limits>
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::TimeRange;
  using KKGeom::Calo;
  using KinKal::Annulus;
  using KinKal::Intersection;
  using AnnPtr = std::shared_ptr<KinKal::Annulus>;
  using CaloDisk = std::vector<AnnPtr>;
  using CylPtr = std::shared_ptr<KinKal::Cylinder>;
  class ExtrapolateCalo {
    public:
      ExtrapolateCalo() : maxDt_(-1.0), dptol_(1e10), intertol_(1e10),
      d0zmin_(std::numeric_limits<float>::max()),
      d0zmax_(-std::numeric_limits<float>::max()),
      d1zmin_(std::numeric_limits<float>::max()),
      d1zmax_(-std::numeric_limits<float>::max()),
      rmin_(std::numeric_limits<float>::max()),
      rmax_(-std::numeric_limits<float>::max()),
      debug_(0){}

      ExtrapolateCalo(double maxdt, double dptol, double intertol, Calo const& calo, int debug=0) :
        maxDt_(maxdt), dptol_(dptol), intertol_(intertol),
        d0zmin_( (calo.EMC_Disk_0_Outer().center() - calo.EMC_Disk_0_Outer().axis()*calo.EMC_Disk_0_Outer().halfLength()).Z()),
        d0zmax_( (calo.EMC_Disk_0_Outer().center() + calo.EMC_Disk_0_Outer().axis()*calo.EMC_Disk_0_Outer().halfLength()).Z()),
        d1zmin_( (calo.EMC_Disk_1_Outer().center() - calo.EMC_Disk_1_Outer().axis()*calo.EMC_Disk_1_Outer().halfLength()).Z()),
        d1zmax_( (calo.EMC_Disk_1_Outer().center() + calo.EMC_Disk_1_Outer().axis()*calo.EMC_Disk_1_Outer().halfLength()).Z()),
        rmin_( calo.EMC_Disk_0_Inner().radius()), rmax_(calo.EMC_Disk_0_Outer().radius()),
        disks_(2),debug_(debug) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      double d0zmin() const { return d0zmin_; }
      double d0zmax() const { return d0zmax_; }
      double d1zmin() const { return d1zmin_; }
      double d1zmax() const { return d1zmax_; }
      double rmin() const { return rmin_; }
      double rmax() const { return rmax_; }
      auto const& disks () const { return disks_; }
      auto const& intersection() const { return inter_; }

      int debug() const { return debug_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
      // reset between tracks
      void reset() const { inter_ = Intersection(); sid_ = SurfaceId(); ann_ = AnnPtr();}
      // find the nearest disk to a z positionin a given z direction
      size_t nearestDisk(double zpos, double zdir) const;
    private:
      double maxDt_; // maximum extrapolation time
      double dptol_; // fractional momentum tolerance
      double intertol_; // intersection tolerance (mm)
      mutable Intersection inter_; // cache of most recent intersection
      mutable SurfaceId sid_; // cache of most recent intersection disk SID
      mutable AnnPtr ann_; // cache of most recent intersection disk surface
      // cache of front and back Z positions
      double d0zmin_, d0zmax_; // z range of disk0 volume
      double d1zmin_, d1zmax_; // z range of disk1 volume
      double rmin_, rmax_; // inner and outer radii of the anuli
      CaloDisk disks_; // disks
      int debug_; // debug level
  };

  template <class KTRAJ> bool ExtrapolateCalo::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    // we are answering the question: did the segment last added to this extrapolated track hit a calo disk or not? If so we are done
    // extrapolating (for now) and we want to find all the intersections in that piece. If not, and if we're still inside or heading towards the
    // disks, keep going.
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
    if(debug_ > 2)std::cout << "calo extrap tdir " << tdir << " start z " << spos.Z() << " end z " << epos.Z() << " zvel " << vel.Z() << " rho " << spos.Rho() << std::endl;
    // stop if the particle is heading away from the calo
    if( (vel.Z() > 0 && spos.Z() > d1zmax_ ) || (vel.Z() < 0 && spos.Z() < d0zmin_)){
      reset(); // clear any cache
      if(debug_ > 1)std::cout << "Heading away from calo: done" << std::endl;
      return false;
    }
    // if the particle is going in the right direction but haven't yet reached the calo in Z just keep going
    if( (vel.Z() > 0 && epos.Z() < d0zmin_) || (vel.Z() < 0 && epos.Z() > d1zmax_) ){
      reset();
      if(debug_ > 2)std::cout << "Heading towards calo, z " << spos.Z()<< std::endl;
      return true;
    }
    // if we get to here we are in the correct Z range. Test disks.
    int idisk = nearestDisk(spos.Z(),vel.Z());
    if(idisk >= (int)disks_.size())return true;
    if(debug_ > 2)std::cout << "Looping on disks " << std::endl;
    int ddisk = vel.Z() > 0.0 ? 1 : -1; // iteration direction
    auto trange = tdir == TimeDir::forwards ? TimeRange(stime,ktraj.range().end()) : TimeRange(ktraj.range().begin(),stime);
    // loop over disks in the z range of this piece
    while(idisk >= 0 && idisk < (int)disks_.size() && (disks_[idisk]->center().Z() - epos.Z())*ddisk < 0.0){
      auto diskptr = disks_[idisk];
      if(debug_ > 2)std::cout << "disk " << idisk << " z " << diskptr->center().Z() << std::endl;
      auto newinter = KinKal::intersect(ktraj,*diskptr,trange,intertol_,tdir);
      if(debug_ > 2)std::cout << "calo disk inter " << newinter  << std::endl;
      if(newinter.good()){
        // update the cache
        inter_ = newinter;
        ann_ = disks_[idisk];
        //sid_ = SurfaceId(SurfaceIdEnum::calo_Foils,idisk); //FIXME
        if(debug_ > 0)std::cout << "Good calo disk " << newinter << " sid " << sid_ << std::endl;
        return false;
      }
      idisk += ddisk; // otherwise continue loopin on disks
    }
    // no more intersections: keep extending in Z till we clear the calo
    reset();
    if(debug_ > 1)std::cout << "Extrapolating to calo edge, z " << spos.Z() << std::endl;
    if(vel.Z() > 0.0)
      return spos.Z() < d1zmax_;
    else
      return spos.Z() > d0zmin_;
  }

  size_t ExtrapolateCalo::nearestDisk(double zpos, double zvel) const {
    size_t retval = disks_.size();
    if(zvel > 0.0){ // going forwards in z
      for(auto idisk= disks_.begin(); idisk != disks_.end(); idisk++){
        auto const& diskptr = *idisk;
        if(diskptr->center().Z() > zpos){
          retval = std::distance(disks_.begin(),idisk);
          break;
        }
      }
    } else {
      for(auto idisk= disks_.rbegin(); idisk != disks_.rend(); idisk++){
        auto const& diskptr = *idisk;
        if(diskptr->center().Z() < zpos){
          auto jdisk = idisk.base()-1; // points to the equivalent forwards object
          retval = std::distance(disks_.begin(),jdisk);
          break;
        }
      }
    }
    return retval;
  }

}
#endif
