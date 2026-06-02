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
#include "Offline/Mu2eKinKal/inc/ExtrapolateCRV.hh"
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
#include "Offline/KinKalGeom/inc/KKMaterial.hh"
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
      explicit KKExtrap(KKExtrapConfig const& exconfig);
      // extrapolation functions; these are templated on the type of trajectory
      template <class KTRAJ> void extrapolate(KKTrack<KTRAJ>& ktrk) const;
      template <class KTRAJ> bool extrapolateIPA(KKTrack<KTRAJ>& ktrk,TimeDir trkdir) const;
      template <class KTRAJ> bool extrapolateST(KKTrack<KTRAJ>& ktrk,TimeDir trkdir) const;
      template <class KTRAJ> bool extrapolateTracker(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> bool extrapolateTSDA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void extrapolateCRV(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void extrapolateOPA(KKTrack<KTRAJ>& ktrk, double tstart, TimeDir tdir) const;
      template <class KTRAJ> void toTrackerEnds(KKTrack<KTRAJ>& ktrk) const;
        // TODO add DS and shielding material

    private:
      int debug_;
      double btol_, intertol_, maxdt_, maxdtstep_, minv_;
      bool backToTracker_, extrapolateOPA_, toTrackerEnds_, upstream_, toCRV_;
      double ipathick_ = 0.511; // ipa thickness: should come from geometry service TODO
      double stthick_ = 0.1056; // st foil thickness: should come from geometry service TODO
  };

  KKExtrap::KKExtrap(KKExtrapConfig const& extrapconfig) :
    debug_(extrapconfig.Debug()),
    btol_(extrapconfig.btol()),
    intertol_(extrapconfig.interTol()),
    maxdt_(extrapconfig.MaxDt()),
    maxdtstep_(extrapconfig.MaxDtStep()),
    minv_(extrapconfig.MinV()),
    backToTracker_(extrapconfig.BackToTracker()),
    extrapolateOPA_(extrapconfig.ToOPA()),
    toTrackerEnds_(extrapconfig.ToTrackerEnds()),
    upstream_(extrapconfig.Upstream()),
    toCRV_(extrapconfig.ToCRV())
  {}

  template <class KTRAJ> void KKExtrap::extrapolate(KKTrack<KTRAJ>& ktrk) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    // define the time direction according to the fit direction inside the tracker
    auto const& ftraj = ktrk.fitTraj();
    auto dir0 = ftraj.direction(ftraj.t0());
    TimeDir tdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
    if(toTrackerEnds_)toTrackerEnds(ktrk);
    if(upstream_){
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
        ExtrapolateToZ trackerFront(maxdt_,maxdtstep_,btol_,kkg_h->tracker()->front().center().Z(),debug_);
        if(backToTracker_)ktrk.extrapolate(tdir,trackerFront);
      }
      // optionally test for intersection with the OPA
      if(extrapolateOPA_)extrapolateOPA(ktrk,starttime,tdir);
    }
    // try CRV in both directions
    if(toCRV_){
      extrapolateCRV(ktrk,TimeDir::backwards);
      extrapolateCRV(ktrk,TimeDir::forwards);
    }
  }

  template <class KTRAJ> void KKExtrap::toTrackerEnds(KKTrack<KTRAJ>& ktrk) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;

    ExtrapolateToZ trackerFront(maxdt_,maxdtstep_,btol_,kkg_h->tracker()->front().center().Z(),debug_);
    ExtrapolateToZ trackerBack(maxdt_,maxdtstep_,btol_,kkg_h->tracker()->back().center().Z(),debug_);
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
    auto midinter = KinKal::intersect(ftraj,kkg_h->tracker()->middle(),ftraj.range(),intertol_);
    if(midinter.good()) ktrk.addIntersection(tt_mid,midinter);
    if(tofront){
      // check the front piece first; that is usually correct
      // track extrapolation to the front succeeded, but the intersection failed. Use the last trajectory to force an intersection
      auto fhel = fronttdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
      auto frontinter = KinKal::intersect(fhel,kkg_h->tracker()->front(),fhel.range(),intertol_,fronttdir);
      if(!frontinter.good()){
        // start from the middle
        TimeRange frange = ftraj.range();
        if(midinter.good())frange = fronttdir == TimeDir::forwards ? TimeRange(midinter.time_,ftraj.range().end()) : TimeRange(ftraj.range().begin(),midinter.time_);
        frontinter = KinKal::intersect(ftraj,kkg_h->tracker()->front(),frange,intertol_,fronttdir);
      }
      if(frontinter.good()) ktrk.addIntersection(tt_front,frontinter);
    }
    if(toback){
      // start from the middle
      TimeRange brange = ftraj.range();
      if(midinter.good())brange = backtdir == TimeDir::forwards ? TimeRange(midinter.time_,ftraj.range().end()) : TimeRange(ftraj.range().begin(),midinter.time_);
      auto backinter = KinKal::intersect(ftraj,kkg_h->tracker()->back(),brange,intertol_,backtdir);
      if(backinter.good())ktrk.addIntersection(tt_back,backinter);
    }
  }

  template <class KTRAJ> bool KKExtrap::extrapolateIPA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    GeomHandle<mu2e::KKMaterial> kkmat_h;
    using KKIPAXING = KKShellXing<KTRAJ,KinKal::Cylinder>;
    using KKIPAXINGPTR = std::shared_ptr<KKIPAXING>;
    // extraplate the fit through the IPA. This will add material effects for each intersection. It will continue till the
    // track exits the IPA
    ExtrapolateIPA extrapIPA(maxdt_,maxdtstep_,btol_,intertol_,kkg_h->DS()->innerProtonAbsorberPtr(),debug_);
    if(extrapIPA.debug() > 2)std::cout << "extrapolating to IPA " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId IPASID("IPA");
    double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto startdir = ftraj.direction(starttime);
    do {
      ktrk.extrapolate(tdir,extrapIPA);
      if(extrapIPA.intersection().good()){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        auto const& IPA = kkg_h->DS()->innerProtonAbsorberPtr();
        KKIPAXINGPTR ipaxingptr = std::make_shared<KKIPAXING>(IPA,IPASID,*kkmat_h->IPAMaterial(),extrapIPA.intersection(),reftrajptr,ipathick_,extrapIPA.interTolerance());
        if(extrapIPA.debug() > 2){
          double dmom, paramomvar, perpmomvar;
          ipaxingptr->materialEffects(dmom,paramomvar,perpmomvar);
          std::cout << "IPA Xing dmom " << dmom << " para momsig " << sqrt(paramomvar) << " perp momsig " << sqrt(perpmomvar) << std::endl;
          std::cout << " before append mom = " << reftrajptr->momentum();
        }
        ktrk.addIPAXing(ipaxingptr,tdir);
        if(extrapIPA.debug() > 2){
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
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    GeomHandle<mu2e::KKMaterial> kkmat_h;

    // extraplate the fit through the ST. This will add material effects for each foil intersection. It will continue till the
    // track exits the ST in Z
    ExtrapolateST extrapST(maxdt_,maxdtstep_,btol_,intertol_,*kkg_h->ST(),debug_);
    auto const& ftraj = ktrk.fitTraj();
    double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto startdir = ftraj.direction(starttime);
    if(extrapST.debug() > 2)std::cout << "extrapolating to ST " << std::endl;
    do {
      ktrk.extrapolate(tdir,extrapST);
      if(extrapST.intersection().good()){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        KKSTXINGPTR stxingptr = std::make_shared<KKSTXING>(extrapST.foil(),extrapST.foilId(),*kkmat_h->STMaterial(),extrapST.intersection(),reftrajptr,stthick_,extrapST.interTolerance());
        if(extrapST.debug() > 2){
          double dmom, paramomvar, perpmomvar;
          stxingptr->materialEffects(dmom,paramomvar,perpmomvar);
          std::cout << "ST Xing dmom " << dmom << " para momsig " << sqrt(paramomvar) << " perp momsig " << sqrt(perpmomvar) << std::endl;
          std::cout << " before append mom = " << reftrajptr->momentum();
        }
        // Add the xing. This truncates the fit
        ktrk.addSTXing(stxingptr,tdir);
        if(extrapST.debug() > 2){
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
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    ExtrapolateToZ trackerFront(maxdt_,maxdtstep_,btol_,kkg_h->tracker()->front().center().Z(),debug_);
    if(trackerFront.debug() > 2)std::cout << "extrapolating to Tracker " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId TrackerSID("TT_Front");
    ktrk.extrapolate(tdir,trackerFront);
    // the last piece appended should cover the necessary range
    auto const& ktraj = tdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
    auto trkfrontinter = KinKal::intersect(ftraj,kkg_h->tracker()->front(),ktraj.range(),intertol_,tdir);
    if(trkfrontinter.onsurface_){ // dont worry about bounds here
      ktrk.addIntersection(TrackerSID,trkfrontinter);
      return true;
    }
    return false;
  }

  template <class KTRAJ> bool KKExtrap::extrapolateTSDA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    ExtrapolateToZ TSDA(maxdt_,maxdtstep_,btol_,kkg_h->DS()->upstreamAbsorber().center().Z(),debug_);
    if(TSDA.debug() > 2)std::cout << "extrapolating to TSDA " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId TSDASID("TSDA");
    ktrk.extrapolate(tdir,TSDA);
    // if we reflected we're done. Otherwize, save the TSDA intersection
    double tend = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto epos = ftraj.position3(tend);
    bool retval = epos.Z() < TSDA.zVal();
    if(retval){
      auto const& ktraj = tdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
      auto tsdainter = KinKal::intersect(ftraj,kkg_h->DS()->upstreamAbsorber(),ktraj.range(),intertol_,tdir);
      if(tsdainter.onsurface_)ktrk.addIntersection(TSDASID,tsdainter);
    }
    return retval;
  }

  //extrapolate to the OPA (DS entrance). only 1 of these is possible
  template <class KTRAJ> void KKExtrap::extrapolateOPA(KKTrack<KTRAJ>& ktrk, double tstart, TimeDir tdir) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId OPASID("OPA");
    TimeRange trange = tdir == TimeDir::forwards ? TimeRange(tstart,ftraj.range().end()) : TimeRange(ftraj.range().begin(),tstart);
    auto opainter = KinKal::intersect(ftraj,kkg_h->DS()->outerProtonAbsorber(),trange,intertol_,tdir);
    if(opainter.good()){
      ktrk.addIntersection(OPASID,opainter);
    }
  }

  //extrapolate to the CRV sectors. This only makes sense for KKLine or CentralHelix
  template <class KTRAJ> void KKExtrap::extrapolateCRV(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    using KKCRVXING = KKShellXing<KTRAJ,KinKal::Rectangle>;
    using KKCRVXINGPTR = std::shared_ptr<KKCRVXING>;
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    GeomHandle<mu2e::KKMaterial> kkmat_h;
    // extrapolate to the extracted CRV. Loop to cover multiple intersections
    auto extrapCRV = ExtrapolateCRV(maxdt_,maxdtstep_,btol_,intertol_,minv_,*kkg_h->CRV(),debug_);
    if(debug_ > 5){std::cout << "Extrapolating to CRV with " << extrapCRV.sectors().size() << " sectors" << std::endl;
      for(auto const& sector : kkg_h->CRV()->sectors()) {
        std::cout << sector.sname_ << " position " << sector.sector_->center() << " halfwidth " << sector.whw_ << std::endl;
      }
    }
    auto const& ftraj = ktrk.fitTraj();
    do {
      // iterate until we no longer hit CRV modules
      ktrk.extrapolate(tdir,extrapCRV);
      if(debug_ > 0) std::cout << "Found " << extrapCRV.intersections().size() << " CRV intersections " << std::endl;
      for(auto const& inter : extrapCRV.intersections()){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        auto crvxingptr = std::make_shared<KKCRVXING>(extrapCRV.sector((size_t)inter.isect_).sector_,
            SurfaceId(kkg_h->CRV()->sectorName(inter.isect_)),*kkmat_h->CRVMaterial(),inter.inter_,reftrajptr,
            2*inter.whw_, extrapCRV.interTolerance());
        ktrk.addCRVXing(crvxingptr,tdir);
        if(debug_ > 1) std::cout << "Good CRV " << inter.inter_ << ftraj.range() << std::endl;
      }
    } while(extrapCRV.intersections().size()>0);
  }

}
#endif
