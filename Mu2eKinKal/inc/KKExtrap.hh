#ifndef Mu2eKinKal_KKExtrap_hh
#define Mu2eKinKal_KKExtrap_hh
//
// Helper class for extrapolating fits through materials and surfaces.
// original author: Dave Brown (LBN), 9/8/2025
//
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/KinKalGeom/inc/SurfaceMap.hh"
#include "Offline/KinKalGeom/inc/Tracker.hh"
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateToZ.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateIPA.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateST.hh"
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
#include "Offline/Mu2eKinKal/inc/KKMaterial.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
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
      template <class KTRAJ> bool extrapolateTSDA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void toOPA(KKTrack<KTRAJ>& ktrk, double tstart, TimeDir tdir) const;

    private:
      int debug_;
      double btol_, intertol_, maxdt_;
      SurfaceMap smap_;
      KKMaterial const& kkmat_;
      AnnPtr tsdaptr_;
      DiskPtr trkfrontptr_, trkmidptr_, trkbackptr_;
      FruPtr opaptr_;
      bool backToTracker_, toOPA_, toTrackerEnds_, upstream_;
      ExtrapolateToZ TSDA_, trackerFront_, trackerBack_; // extrapolation predicate based on Z values
      ExtrapolateIPA extrapIPA_; // extrapolation to intersections with the IPA
      ExtrapolateST extrapST_; // extrapolation to intersections with the ST
      double ipathick_ = 0.511; // ipa thickness: should come from geometry service TODO
      double stthick_ = 0.1056; // st foil thickness: should come from geometry service TODO
  };

  KKExtrap::KKExtrap(KKExtrapConfig const& extrapconfig,KKMaterial const& kkmat) :
    debug_(extrapconfig.Debug()),
    btol_(extrapconfig.btol()),
    intertol_(extrapconfig.interTol()),
    maxdt_(extrapconfig.MaxDt()),
    kkmat_(kkmat),
    tsdaptr_(smap_.DS().upstreamAbsorberPtr()),
    trkfrontptr_(smap_.tracker().frontPtr()),
    trkmidptr_(smap_.tracker().middlePtr()),
    trkbackptr_(smap_.tracker().backPtr()),
    opaptr_(smap_.DS().outerProtonAbsorberPtr()),
    backToTracker_(extrapconfig.BackToTracker()),
    toOPA_(extrapconfig.ToOPA()),
    toTrackerEnds_(extrapconfig.ToTrackerEnds()),
    upstream_(extrapconfig.Upstream()),
    TSDA_(maxdt_,btol_,smap_.DS().upstreamAbsorber().center().Z(),debug_),
    trackerFront_(maxdt_,btol_,smap_.tracker().front().center().Z(),debug_),
    trackerBack_(maxdt_,btol_,smap_.tracker().back().center().Z(),debug_),
    extrapIPA_(maxdt_,btol_,intertol_,smap_.DS().innerProtonAbsorberPtr(),debug_),
    extrapST_(maxdt_,btol_,intertol_,smap_.ST(),debug_)
  {
    if(debug_ > 0)std::cout << "IPA limits z " << extrapIPA_.zmin() << " " << extrapIPA_.zmax() << std::endl;
    if(debug_ > 0)std::cout << "ST limits z " << extrapST_.zmin() << " " << extrapST_.zmax() << " r " << extrapST_.rmin() << " " << extrapST_.rmax() << std::endl;
  }

  template <class KTRAJ> void KKExtrap::extrapolate(KKTrack<KTRAJ>& ktrk) const {
    // define the time direction according to the fit direction inside the tracker
    auto const& ftraj = ktrk.fitTraj();
    if(toTrackerEnds_)toTrackerEnds(ktrk);
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
        if(backToTracker_)ktrk.extrapolate(tdir,trackerFront_);
      }
      // optionally test for intersection with the OPA
      if(toOPA_)toOPA(ktrk,starttime,tdir);
    }
  }

  template <class KTRAJ> void KKExtrap::toTrackerEnds(KKTrack<KTRAJ>& ktrk) const {
    // time direction to reach the bounding surfaces from the active region depends on the z momentum. This calculation assumes the particle doesn't
    // reflect inside the tracker volumei
    auto const& ftraj = ktrk.fitTraj();
    auto dir0 = ftraj.direction(ftraj.t0());
    TimeDir fronttdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
    TimeDir backtdir = (dir0.Z() > 0) ? TimeDir::forwards : TimeDir::backwards;
    auto tofront = ktrk.extrapolate(fronttdir,trackerFront_);
    auto toback = ktrk.extrapolate(backtdir,trackerBack_);
    // record the standard tracker intersections
    static const SurfaceId tt_front("TT_Front");
    static const SurfaceId tt_mid("TT_Mid");
    static const SurfaceId tt_back("TT_Back");

    // start with the middle
    auto midinter = KinKal::intersect(ftraj,*trkmidptr_,ftraj.range(),intertol_);
    if(midinter.good()) ktrk.addIntersection(tt_mid,midinter);
    if(tofront){
      // check the front piece first; that is usually correct
      // track extrapolation to the front succeeded, but the intersection failed. Use the last trajectory to force an intersection
      auto fhel = fronttdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
      auto frontinter = KinKal::intersect(fhel,*trkfrontptr_,fhel.range(),intertol_,fronttdir);
      if(!frontinter.good()){
        // start from the middle
        TimeRange frange = ftraj.range();
        if(midinter.good())frange = fronttdir == TimeDir::forwards ? TimeRange(midinter.time_,ftraj.range().end()) : TimeRange(ftraj.range().begin(),midinter.time_);
        frontinter = KinKal::intersect(ftraj,*trkfrontptr_,frange,intertol_,fronttdir);
      }
      if(frontinter.good()) ktrk.addIntersection(tt_front,frontinter);
    }
    if(toback){
      // start from the middle
      TimeRange brange = ftraj.range();
      if(midinter.good())brange = backtdir == TimeDir::forwards ? TimeRange(midinter.time_,ftraj.range().end()) : TimeRange(ftraj.range().begin(),midinter.time_);
      auto backinter = KinKal::intersect(ftraj,*trkbackptr_,brange,intertol_,backtdir);
      if(backinter.good())ktrk.addIntersection(tt_back,backinter);
    }
  }

  template <class KTRAJ> bool KKExtrap::extrapolateIPA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    using KKIPAXING = KKShellXing<KTRAJ,KinKal::Cylinder>;
    using KKIPAXINGPTR = std::shared_ptr<KKIPAXING>;

    if(extrapIPA_.debug() > 0)std::cout << "extrapolating to IPA " << std::endl;
    // extraplate the fit through the IPA. This will add material effects for each intersection. It will continue till the
    // track exits the IPA
    extrapIPA_.reset();
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId IPASID("IPA");
    double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto startdir = ftraj.direction(starttime);
    do {
      ktrk.extrapolate(tdir,extrapIPA_);
      if(extrapIPA_.intersection().good()){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        auto const& IPA = smap_.DS().innerProtonAbsorberPtr();
        KKIPAXINGPTR ipaxingptr = std::make_shared<KKIPAXING>(IPA,IPASID,*kkmat_.IPAMaterial(),extrapIPA_.intersection(),reftrajptr,ipathick_,extrapIPA_.interTolerance());
        if(extrapIPA_.debug() > 0){
          double dmom, paramomvar, perpmomvar;
          ipaxingptr->materialEffects(dmom,paramomvar,perpmomvar);
          std::cout << "IPA Xing dmom " << dmom << " para momsig " << sqrt(paramomvar) << " perp momsig " << sqrt(perpmomvar) << std::endl;
          std::cout << " before append mom = " << reftrajptr->momentum();
        }
        ktrk.addIPAXing(ipaxingptr,tdir);
        if(extrapIPA_.debug() > 0){
          auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
          std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
        }
      }
    } while(extrapIPA_.intersection().good());
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

    // extraplate the fit through the ST. This will add material effects for each foil intersection. It will continue till the
    // track exits the ST in Z
    extrapST_.reset();
    auto const& ftraj = ktrk.fitTraj();
    double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto startdir = ftraj.direction(starttime);
    if(extrapST_.debug() > 0)std::cout << "extrapolating to ST " << std::endl;
    do {
      ktrk.extrapolate(tdir,extrapST_);
      if(extrapST_.intersection().good()){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        KKSTXINGPTR stxingptr = std::make_shared<KKSTXING>(extrapST_.foil(),extrapST_.foilId(),*kkmat_.STMaterial(),extrapST_.intersection(),reftrajptr,stthick_,extrapST_.interTolerance());
        if(extrapST_.debug() > 0){
          double dmom, paramomvar, perpmomvar;
          stxingptr->materialEffects(dmom,paramomvar,perpmomvar);
          std::cout << "ST Xing dmom " << dmom << " para momsig " << sqrt(paramomvar) << " perp momsig " << sqrt(perpmomvar) << std::endl;
          std::cout << " before append mom = " << reftrajptr->momentum();
        }
        // Add the xing. This truncates the fit
        ktrk.addSTXing(stxingptr,tdir);
        if(extrapST_.debug() > 0){
          auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
          std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
        }
      }
    } while(extrapST_.intersection().good());
    // check if the particle exited in the same physical direction or not (reflection)
    double endtime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto enddir = ftraj.direction(endtime);
    if(enddir.Z() * startdir.Z() > 0.0){
      return true;
    }
    return false;
  }

  template <class KTRAJ> bool KKExtrap::extrapolateTracker(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    if(trackerFront_.debug() > 0)std::cout << "extrapolating to Tracker " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId TrackerSID("TT_Front");
    ktrk.extrapolate(tdir,trackerFront_);
    // the last piece appended should cover the necessary range
    auto const& ktraj = tdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
    auto trkfrontinter = KinKal::intersect(ftraj,*trkfrontptr_,ktraj.range(),intertol_,tdir);
    if(trkfrontinter.onsurface_){ // dont worry about bounds here
      ktrk.addIntersection(TrackerSID,trkfrontinter);
      return true;
    }
    return false;
  }

  template <class KTRAJ> bool KKExtrap::extrapolateTSDA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    if(TSDA_.debug() > 0)std::cout << "extrapolating to TSDA " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId TSDASID("TSDA");
    ktrk.extrapolate(tdir,TSDA_);
    // if we reflected we're done. Otherwize, save the TSDA intersection
    double tend = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto epos = ftraj.position3(tend);
    bool retval = epos.Z() < TSDA_.zVal();
    if(retval){
      auto const& ktraj = tdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
      auto tsdainter = KinKal::intersect(ftraj,*tsdaptr_,ktraj.range(),intertol_,tdir);
      if(tsdainter.onsurface_)ktrk.addIntersection(TSDASID,tsdainter);
    }
    return retval;
  }

  template <class KTRAJ> void KKExtrap::toOPA(KKTrack<KTRAJ>& ktrk, double tstart, TimeDir tdir) const {
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId OPASID("OPA");
    TimeRange trange = tdir == TimeDir::forwards ? TimeRange(tstart,ftraj.range().end()) : TimeRange(ftraj.range().begin(),tstart);
    auto opainter = KinKal::intersect(ftraj,*opaptr_,trange,intertol_,tdir);
    if(opainter.good()){
      ktrk.addIntersection(OPASID,opainter);
    }
  }
}
#endif
