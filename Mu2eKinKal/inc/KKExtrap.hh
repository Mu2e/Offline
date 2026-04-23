#ifndef Mu2eKinKal_KKExtrap_hh
#define Mu2eKinKal_KKExtrap_hh
//
// Helper class for extrapolating fits through materials and surfaces.
// original author: Dave Brown (LBN), 9/8/2025
//
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/KinKalGeom/inc/Tracker.hh"
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateToZ.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateIPA.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateST.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateCalo.hh"
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
#include "Offline/Mu2eKinKal/inc/KKMaterial.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
namespace mu2e {
  using KinKal::VEC3;
  using KinKal::TimeDir;
  using SurfacePtr = std::shared_ptr<KinKal::Surface>;
  using CylPtr = std::shared_ptr<KinKal::Cylinder>;
  using DiskPtr = std::shared_ptr<KinKal::Disk>;
  using FruPtr = std::shared_ptr<KinKal::Frustrum>;
  using AnnPtr = std::shared_ptr<KinKal::Annulus>;
  using Mu2eKinKal::KKExtrapConfig;

  class KKExtrap {
    public:
      explicit KKExtrap(KKExtrapConfig const& exconfig,KKMaterial const& kkmat);
      // extrapolation functions; these are templated on the type of trajectory
      template <class KTRAJ> void extrapolate(KKTrack<KTRAJ>& ktrk) const;
      template <class KTRAJ> void toTrackerEnds(KKTrack<KTRAJ>& ktrk) const;
      template <class KTRAJ> bool extrapolateIPA(KKTrack<KTRAJ>& ktrk,TimeDir trkdir) const;
      template <class KTRAJ> bool extrapolateST(KKTrack<KTRAJ>& ktrk,TimeDir trkdir) const;
      template <class KTRAJ> bool extrapolateTracker(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void toCaloD0(KKTrack<KTRAJ>& ktrk) const;
      template <class KTRAJ> void toCaloD1(KKTrack<KTRAJ>& ktrk) const;
      template <class KTRAJ> bool extrapolateTSDA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void toOPA(KKTrack<KTRAJ>& ktrk, double tstart, TimeDir tdir) const;

    private:
      int debug_;
      double btol_, intertol_, maxdt_;
      KKMaterial const& kkmat_;
      AnnPtr tsdaptr_;
      DiskPtr trkfrontptr_, trkmidptr_, trkbackptr_;
      FruPtr opaptr_;
      bool backToTracker_, toOPA_, toTrackerEnds_, upstream_;
      bool toCaloD0_, toCaloD1_;
      // calo surfaces and predicates
      DiskPtr calod0frontptr_, calod0backptr_, calod1frontptr_, calod1backptr_;
      mutable ExtrapolateToZ calod0Front_, calod0Back_, calod1Front_, calod1Back_;
      mutable ExtrapolateToZ TSDA_, trackerFront_, trackerBack_; // extrapolation predicate based on Z values
      mutable ExtrapolateIPA extrapIPA_; // extrapolation to intersections with the IPA
      mutable ExtrapolateST extrapST_; // extrapolation to intersections with the ST
      mutable bool geom_initialized_ = false;
      void initGeometry() const; // lazy initialization method
      double ipathick_ = 0.511; // ipa thickness: should come from geometry service TODO
      double stthick_ = 0.1056; // st foil thickness: should come from geometry service TODO
  };

  KKExtrap::KKExtrap(KKExtrapConfig const& extrapconfig,KKMaterial const& kkmat) :
    debug_(extrapconfig.Debug()),
    btol_(extrapconfig.btol()),
    intertol_(extrapconfig.interTol()),
    maxdt_(extrapconfig.MaxDt()),
    kkmat_(kkmat),
    tsdaptr_(nullptr),
    trkfrontptr_(nullptr),
    trkmidptr_(nullptr),
    trkbackptr_(nullptr),
    opaptr_(nullptr),
    backToTracker_(extrapconfig.BackToTracker()),
    toOPA_(extrapconfig.ToOPA()),
    toTrackerEnds_(extrapconfig.ToTrackerEnds()),
    upstream_(extrapconfig.Upstream()),
    toCaloD0_(extrapconfig.ToCaloD0()),
    toCaloD1_(extrapconfig.ToCaloD1()),
    geom_initialized_(false)
  {
    // geometry initialization is deferred to first use via initGeometry()
  }

  inline void KKExtrap::initGeometry() const {
    if(geom_initialized_) return; // already initialized
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;

    const_cast<AnnPtr&>(tsdaptr_) = kkg.DS().upstreamAbsorberPtr();
    const_cast<DiskPtr&>(trkfrontptr_) = kkg.tracker().frontPtr();
    const_cast<DiskPtr&>(trkmidptr_) = kkg.tracker().middlePtr();
    const_cast<DiskPtr&>(trkbackptr_) = kkg.tracker().backPtr();
    const_cast<FruPtr&>(opaptr_) = kkg.DS().outerProtonAbsorberPtr();

    // calo surfaces
    const_cast<DiskPtr&>(calod0frontptr_) = kkg.calo().EMC_Disk_0_FrontPtr();
    const_cast<DiskPtr&>(calod0backptr_) = kkg.calo().EMC_Disk_0_BackPtr();
    const_cast<DiskPtr&>(calod1frontptr_) = kkg.calo().EMC_Disk_1_FrontPtr();
    const_cast<DiskPtr&>(calod1backptr_) = kkg.calo().EMC_Disk_1_BackPtr();

    // initialize predicates
    const_cast<ExtrapolateToZ&>(TSDA_) = ExtrapolateToZ(maxdt_,btol_,kkg.DS().upstreamAbsorber().center().Z(),debug_);
    const_cast<ExtrapolateToZ&>(trackerFront_) = ExtrapolateToZ(maxdt_,btol_,kkg.tracker().front().center().Z(),debug_);
    const_cast<ExtrapolateToZ&>(trackerBack_) = ExtrapolateToZ(maxdt_,btol_,kkg.tracker().back().center().Z(),debug_);
    const_cast<ExtrapolateIPA&>(extrapIPA_) = ExtrapolateIPA(maxdt_,btol_,intertol_,kkg.DS().innerProtonAbsorberPtr(),debug_);
    const_cast<ExtrapolateST&>(extrapST_) = ExtrapolateST(maxdt_,btol_,intertol_,kkg.ST(),debug_);

    // calo predicates - use local Z values stored in Calo class
    const_cast<ExtrapolateToZ&>(calod0Front_) = ExtrapolateToZ(maxdt_,btol_,kkg.calo().EMC_Disk_0_Front_Z(),debug_);
    const_cast<ExtrapolateToZ&>(calod0Back_) = ExtrapolateToZ(maxdt_,btol_,kkg.calo().EMC_Disk_0_Back_Z(),debug_);
    const_cast<ExtrapolateToZ&>(calod1Front_) = ExtrapolateToZ(maxdt_,btol_,kkg.calo().EMC_Disk_1_Front_Z(),debug_);
    const_cast<ExtrapolateToZ&>(calod1Back_) = ExtrapolateToZ(maxdt_,btol_,kkg.calo().EMC_Disk_1_Back_Z(),debug_);

    if(debug_ > 0)std::cout << "IPA limits z " << extrapIPA_.zmin() << " " << extrapIPA_.zmax() << std::endl;
    if(debug_ > 0)std::cout << "ST limits z " << extrapST_.zmin() << " " << extrapST_.zmax() << " r " << extrapST_.rmin() << " " << extrapST_.rmax() << std::endl;
    const_cast<bool&>(geom_initialized_) = true;
  }

  template <class KTRAJ> void KKExtrap::extrapolate(KKTrack<KTRAJ>& ktrk) const {
    initGeometry(); // lazy initialization: initialize geometry on first use
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;
    ExtrapolateToZ trackerFront(maxdt_,btol_,kkg.tracker().front().center().Z(),debug_);

    // define the time direction according to the fit direction inside the tracker
    auto const& ftraj = ktrk.fitTraj();
    if(toTrackerEnds_)toTrackerEnds(ktrk);
    if(toCaloD0_)toCaloD0(ktrk);
    if(toCaloD1_)toCaloD1(ktrk);
    if(upstream_){
      auto dir0 = ftraj.direction(ftraj.t0());
      TimeDir tdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
      double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
      // extrapolate through the IPA in this direction.
      bool exitsIPA = extrapolateIPA(ktrk,tdir);

      if(exitsIPA){ // if it exits out the back, extrapolate through the target
        bool exitsST = extrapolateST(ktrk,tdir);
        if(exitsST) { // if it exits out the back, extrapolate to the TSDA (DS rear absorber)
          bool hitTSDA = extrapolateTSDA(ktrk,tdir);
          // if we hit the TSDA we are done. Otherwise if we reflected, go back through the ST
          if(!hitTSDA){ // reflection upstream of the target: go back through the target
            extrapolateST(ktrk,tdir);
            if(backToTracker_){ // optionally extrapolate back through the IPA, then to the tracker entrance
              extrapolateIPA(ktrk,tdir);
              extrapolateTracker(ktrk,tdir);
            }
          }
        } else { // reflection inside the ST; extrapolate back through the IPA, then to the tracker entrance
          if(backToTracker_){
            extrapolateIPA(ktrk,tdir);
            extrapolateTracker(ktrk,tdir);
          }
        }
      } else { // reflection inside the IPA; extrapolate back through the IPA, then to the tracker entrance
        if(backToTracker_)ktrk.extrapolate(tdir,trackerFront);
      }
      // optionally test for intersection with the OPA
      if(toOPA_)toOPA(ktrk,starttime,tdir);
    }
  }

  template <class KTRAJ> void KKExtrap::toTrackerEnds(KKTrack<KTRAJ>& ktrk) const {
    initGeometry();
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;

    ExtrapolateToZ trackerFront(maxdt_,btol_,kkg.tracker().front().center().Z(),debug_);
    ExtrapolateToZ trackerBack(maxdt_,btol_,kkg.tracker().back().center().Z(),debug_);
    // time direction to reach the bounding surfaces from the active region depends on the z momentum. This calculation assumes the particle doesn't
    // reflect inside the tracker volume
    auto const& ftraj = ktrk.fitTraj();
    auto dir0 = ftraj.direction(ftraj.t0());
    TimeDir fronttdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
    TimeDir backtdir = (dir0.Z() > 0) ? TimeDir::forwards : TimeDir::backwards;
    auto tofront = ktrk.extrapolate(fronttdir,trackerFront);
    auto toback = ktrk.extrapolate(backtdir,trackerBack);
    // record the standard tracker intersections
    static const SurfaceId tt_front("TT_Front");
    static const SurfaceId tt_mid("TT_Mid");
    static const SurfaceId tt_back("TT_Back");

    // start with the middle
    auto midinter = KinKal::intersect(ftraj,kkg.tracker().middle(),ftraj.range(),intertol_);
    if(midinter.good()) ktrk.addIntersection(tt_mid,midinter);
    if(tofront){
      // check the front piece first; that is usually correct
      // track extrapolation to the front succeeded, but the intersection failed. Use the last trajectory to force an intersection
      auto fhel = fronttdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
      auto frontinter = KinKal::intersect(fhel,kkg.tracker().front(),fhel.range(),intertol_,fronttdir);
      if(!frontinter.good()){
        // start from the middle
        TimeRange frange = ftraj.range();
        if(midinter.good())frange = fronttdir == TimeDir::forwards ? TimeRange(midinter.time_,ftraj.range().end()) : TimeRange(ftraj.range().begin(),midinter.time_);
        frontinter = KinKal::intersect(ftraj,kkg.tracker().front(),frange,intertol_,fronttdir);
      }
      if(frontinter.good()) ktrk.addIntersection(tt_front,frontinter);
    }
    if(toback){
      // start from the middle
      TimeRange brange = ftraj.range();
      if(midinter.good())brange = backtdir == TimeDir::forwards ? TimeRange(midinter.time_,ftraj.range().end()) : TimeRange(ftraj.range().begin(),midinter.time_);
      auto backinter = KinKal::intersect(ftraj,kkg.tracker().back(),brange,intertol_,backtdir);
      if(backinter.good())ktrk.addIntersection(tt_back,backinter);
    }
  }

    template <class KTRAJ> void KKExtrap::toCaloD0(KKTrack<KTRAJ>& ktrk) const {
    initGeometry();
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;
    auto const& ftraj = ktrk.fitTraj();
    auto dir0 = ftraj.direction(ftraj.t0());
    TimeDir fronttdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
    TimeDir backtdir = (dir0.Z() > 0) ? TimeDir::forwards : TimeDir::backwards;
    std::cout<<"calod0Front_ "<<kkg.calo().EMC_Disk_0_Front().center().Z()<<std::endl;
    std::cout<<"calod0Back_ "<<kkg.calo().EMC_Disk_0_Back().center().Z()<<std::endl;
    auto tocalofront = ktrk.extrapolate(fronttdir,calod0Front_);
    auto tocaloback = ktrk.extrapolate(backtdir,calod0Back_);
    // record the standard tracker intersections
    static const SurfaceId d0_front("EMC_Disk_0_Front");
    static const SurfaceId d0_back("EMC_Disk_0_Back");
    static const SurfaceId d0_inner("EMC_Disk_0_Inner");
    static const SurfaceId d0_outer("EMC_Disk_0_Outer");
    std::cout<<"extrapolating to calo d0"<<std::endl;
    if(tocalofront){
      // use full trajectory range for intersection calculation
      TimeRange frange = ftraj.range();
      auto frontinter = KinKal::intersect(ftraj,*calod0frontptr_,frange,intertol_,fronttdir);
      if(frontinter.good()) ktrk.addIntersection(d0_front,frontinter);
      std::cout<<"to front "<<std::endl;
    }
    if(tocaloback){
      // start from the middle
      TimeRange brange = ftraj.range();
      auto backinter = KinKal::intersect(ftraj,*calod0backptr_,brange,intertol_,backtdir);
      if(backinter.good())ktrk.addIntersection(d0_back,backinter);
      std::cout<<"to back "<<std::endl;
    }
    // get intersections with inner and outer cylinders
    auto innerinter = KinKal::intersect(ftraj,kkg.calo().EMC_Disk_0_Inner(),ftraj.range(),intertol_);
    if(innerinter.good()) ktrk.addIntersection(d0_inner,innerinter);
    auto outerinter = KinKal::intersect(ftraj,kkg.calo().EMC_Disk_0_Outer(),ftraj.range(),intertol_);
    if(outerinter.good()) ktrk.addIntersection(d0_outer,outerinter);
  }

   template <class KTRAJ> void KKExtrap::toCaloD1(KKTrack<KTRAJ>& ktrk) const {
    initGeometry();
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;
    auto const& ftraj = ktrk.fitTraj();
    auto dir0 = ftraj.direction(ftraj.t0());
    TimeDir fronttdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
    TimeDir backtdir = (dir0.Z() > 0) ? TimeDir::forwards : TimeDir::backwards;
    auto tocalofront = ktrk.extrapolate(fronttdir,calod1Front_);
    auto tocaloback = ktrk.extrapolate(backtdir,calod1Back_);
    // record the standard tracker intersections
    static const SurfaceId d1_front("EMC_Disk_1_Front");
    static const SurfaceId d1_back("EMC_Disk_1_Back");
    static const SurfaceId d1_inner("EMC_Disk_1_Inner");
    static const SurfaceId d1_outer("EMC_Disk_1_Outer");

    if(tocalofront){
      // use full trajectory range for intersection calculation
      TimeRange frange = ftraj.range();
      auto frontinter = KinKal::intersect(ftraj,*calod1frontptr_,frange,intertol_,fronttdir);
      if(frontinter.good()) ktrk.addIntersection(d1_front,frontinter);
    }
    if(tocaloback){
      // start from the middle
      TimeRange brange = ftraj.range();
      auto backinter = KinKal::intersect(ftraj,*calod1backptr_,brange,intertol_,backtdir);
      if(backinter.good())ktrk.addIntersection(d1_back,backinter);
    }
    // get intersections with inner and outer cylinders
    auto innerinter = KinKal::intersect(ftraj,kkg.calo().EMC_Disk_1_Inner(),ftraj.range(),intertol_);
    if(innerinter.good()) ktrk.addIntersection(d1_inner,innerinter);
    auto outerinter = KinKal::intersect(ftraj,kkg.calo().EMC_Disk_1_Outer(),ftraj.range(),intertol_);
    if(outerinter.good()) ktrk.addIntersection(d1_outer,outerinter);
  }

  template <class KTRAJ> bool KKExtrap::extrapolateIPA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    initGeometry();
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;
    using KKIPAXING = KKShellXing<KTRAJ,KinKal::Cylinder>;
    using KKIPAXINGPTR = std::shared_ptr<KKIPAXING>;
    // extraplate the fit through the IPA. This will add material effects for each intersection. It will continue till the
    // track exits the IPA
    ExtrapolateIPA extrapIPA(maxdt_,btol_,intertol_,kkg.DS().innerProtonAbsorberPtr(),debug_);
    if(extrapIPA.debug() > 0)std::cout << "extrapolating to IPA " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId IPASID("IPA");
    double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto startdir = ftraj.direction(starttime);
    do {
      ktrk.extrapolate(tdir,extrapIPA);
      if(extrapIPA.intersection().good()){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        auto const& IPA = kkg.DS().innerProtonAbsorberPtr();
        KKIPAXINGPTR ipaxingptr = std::make_shared<KKIPAXING>(IPA,IPASID,*kkmat_.IPAMaterial(),extrapIPA.intersection(),reftrajptr,ipathick_,extrapIPA.interTolerance());
        if(extrapIPA.debug() > 0){
          double dmom, paramomvar, perpmomvar;
          ipaxingptr->materialEffects(dmom,paramomvar,perpmomvar);
          std::cout << "IPA Xing dmom " << dmom << " para momsig " << sqrt(paramomvar) << " perp momsig " << sqrt(perpmomvar) << std::endl;
          std::cout << " before append mom = " << reftrajptr->momentum();
        }
        ktrk.addIPAXing(ipaxingptr,tdir);
        if(extrapIPA.debug() > 0){
          auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
          std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
        }
      }
    } while(extrapIPA.intersection().good());
    // check if the particle exited in the same physical direction or not (reflection)
    double endtime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto enddir = ftraj.direction(endtime);
    if(enddir.Z() * startdir.Z() > 0.0){
      return true;
    }
    return false;
  }

  template <class KTRAJ> bool KKExtrap::extrapolateST(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    using KKSTXING = KKShellXing<KTRAJ,KinKal::Annulus>;
    using KKSTXINGPTR = std::shared_ptr<KKSTXING>;
    initGeometry();
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;

    // extraplate the fit through the ST. This will add material effects for each foil intersection. It will continue till the
    // track exits the ST in Z
    ExtrapolateST extrapST(maxdt_,btol_,intertol_,kkg.ST(),debug_);
    auto const& ftraj = ktrk.fitTraj();
    double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto startdir = ftraj.direction(starttime);
    if(extrapST.debug() > 0)std::cout << "extrapolating to ST " << std::endl;
    do {
      ktrk.extrapolate(tdir,extrapST);
      if(extrapST.intersection().good()){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        KKSTXINGPTR stxingptr = std::make_shared<KKSTXING>(extrapST.foil(),extrapST.foilId(),*kkmat_.STMaterial(),extrapST.intersection(),reftrajptr,stthick_,extrapST.interTolerance());
        if(extrapST.debug() > 0){
          double dmom, paramomvar, perpmomvar;
          stxingptr->materialEffects(dmom,paramomvar,perpmomvar);
          std::cout << "ST Xing dmom " << dmom << " para momsig " << sqrt(paramomvar) << " perp momsig " << sqrt(perpmomvar) << std::endl;
          std::cout << " before append mom = " << reftrajptr->momentum();
        }
        // Add the xing. This truncates the fit
        ktrk.addSTXing(stxingptr,tdir);
        if(extrapST.debug() > 0){
          auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
          std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
        }
      }
    } while(extrapST.intersection().good());
    // check if the particle exited in the same physical direction or not (reflection)
    double endtime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto enddir = ftraj.direction(endtime);
    if(enddir.Z() * startdir.Z() > 0.0){
      return true;
    }
    return false;
  }

  template <class KTRAJ> bool KKExtrap::extrapolateTracker(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    initGeometry();
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;
    ExtrapolateToZ trackerFront(maxdt_,btol_,kkg.tracker().front().center().Z(),debug_);
    if(trackerFront.debug() > 0)std::cout << "extrapolating to Tracker " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId TrackerSID("TT_Front");
    ktrk.extrapolate(tdir,trackerFront);
    // the last piece appended should cover the necessary range
    auto const& ktraj = tdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
    auto trkfrontinter = KinKal::intersect(ftraj,kkg.tracker().front(),ktraj.range(),intertol_,tdir);
    if(trkfrontinter.onsurface_){ // dont worry about bounds here
      ktrk.addIntersection(TrackerSID,trkfrontinter);
      return true;
    }
    return false;
  }

  template <class KTRAJ> bool KKExtrap::extrapolateTSDA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    initGeometry();
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;
    ExtrapolateToZ TSDA(maxdt_,btol_,kkg.DS().upstreamAbsorber().center().Z(),debug_);
    if(TSDA.debug() > 0)std::cout << "extrapolating to TSDA " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId TSDASID("TSDA");
    ktrk.extrapolate(tdir,TSDA);
    // if we reflected we're done. Otherwize, save the TSDA intersection
    double tend = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto epos = ftraj.position3(tend);
    bool retval = epos.Z() < TSDA.zVal();
    if(retval){
      auto const& ktraj = tdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
      auto tsdainter = KinKal::intersect(ftraj,kkg.DS().upstreamAbsorber(),ktraj.range(),intertol_,tdir);
      if(tsdainter.onsurface_)ktrk.addIntersection(TSDASID,tsdainter);
    }
    return retval;
  }

  template <class KTRAJ> void KKExtrap::toOPA(KKTrack<KTRAJ>& ktrk, double tstart, TimeDir tdir) const {
    initGeometry();
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId OPASID("OPA");
    TimeRange trange = tdir == TimeDir::forwards ? TimeRange(tstart,ftraj.range().end()) : TimeRange(ftraj.range().begin(),tstart);
    auto opainter = KinKal::intersect(ftraj,kkg.DS().outerProtonAbsorber(),trange,intertol_,tdir);
    if(opainter.good()){
      ktrk.addIntersection(OPASID,opainter);
    }
  }
}
#endif
