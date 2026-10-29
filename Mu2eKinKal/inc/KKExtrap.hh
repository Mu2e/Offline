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
#include "Offline/Mu2eKinKal/inc/ExtrapolateToTrackerPerimeter.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateIPA.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateST.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateCRVRegion.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateCylinders.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolatePlanes.hh"
#include "Offline/Mu2eKinKal/inc/CaloMaterialLayers.hh"
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
#include "Offline/KinKalGeom/inc/KKMaterial.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "cetlib_except/exception.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include <unordered_set>
#include <map>
#include <optional>
#include <utility>
#include <algorithm>
#include <limits>
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
      // throw at configuration if the fit config's minimum domain step can't fit this extrapolation's maximum step
      void checkConfig(KinKal::Config const& config) const;
      // extrapolation functions; these are templated on the type of trajectory
      template <class KTRAJ> void extrapolate(KKTrack<KTRAJ>& ktrk) const;
      template <class KTRAJ> bool extrapolateIPA(KKTrack<KTRAJ>& ktrk,TimeDir trkdir) const;
      template <class KTRAJ> bool extrapolateST(KKTrack<KTRAJ>& ktrk,TimeDir trkdir) const;
      template <class KTRAJ> bool extrapolateTracker(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> bool extrapolateTSDA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void extrapolateDSMaterial(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void extrapolatePassiveMaterial(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void extrapolateCaloMaterial(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> bool extrapolateCaloSupport(KKTrack<KTRAJ>& ktrk,KKGeom::CaloInnerSupport const& support,double& tstart) const;
      template <class KTRAJ> void extrapolateCRV(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void extrapolateOPA(KKTrack<KTRAJ>& ktrk, double tstart, TimeDir tdir) const;
      template <class KTRAJ> void toTrackerPerimeter(KKTrack<KTRAJ>& ktrk) const;
      template <class KTRAJ> void toCaloDisk(KKTrack<KTRAJ>& ktrk, int idisk, double tref, TimeDir tdir) const;

    private:
      int debug_;
      double btol_, intertol_, maxdt_, maxdtstep_, minv_;
      bool backToTracker_, extrapolateOPA_, toTrackerPerimeter_, upstream_, toCRV_;
      bool toCaloD0_, toCaloD1_, toCaloMaterial_;
      // IPA/ST thicknesses are nominal literals, consistent with the rest of KinKalGeom, whose IPA
      // (single Cylinder) and ST (uniform Annulus foils) surfaces are themselves hard-coded nominal
      // placeholders (see KinKalGeomMaker::makeDS/makeTarget). Sourcing these from the real geometry
      // service is part of the broader "build KinKalGeom from the geometry service" work, not a standalone
      // fix -- doing it for the thicknesses alone while the surfaces stay nominal would be inconsistent. TODO
      double ipathick_ = 0.511; // ipa thickness (mm)
      double stthick_ = 0.1056; // st foil thickness (mm)
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
    toTrackerPerimeter_(extrapconfig.ToTrackerPerimeter()),
    upstream_(extrapconfig.Upstream()),
    toCRV_(extrapconfig.ToCRV()),
    toCaloD0_(extrapconfig.ToCaloD0()),
    toCaloD1_(extrapconfig.ToCaloD1()),
    toCaloMaterial_(extrapconfig.ToCaloMaterial())
  {}

  inline void KKExtrap::checkConfig(KinKal::Config const& config) const {
    if(config.bfcorr_ && config.mindtstep_ > maxdtstep_)
      throw cet::exception("RECO") << "mu2e::KKExtrap: MinDomainTimeStep " << config.mindtstep_
        << " ns exceeds ExtrapolationSettings.MaxDtStep " << maxdtstep_ << " ns" << std::endl;
  }

  template <class KTRAJ> void KKExtrap::extrapolate(KKTrack<KTRAJ>& ktrk) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    // define the time direction according to the fit direction inside the tracker
    auto const& ftraj = ktrk.fitTraj();
    // capture t0 now: once the fit is extrapolated upstream, ftraj.t0() can resolve onto the reflected leg
    double t0 = ftraj.t0();
    auto dir0 = ftraj.direction(t0);
    TimeDir tdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
    if(toTrackerPerimeter_)toTrackerPerimeter(ktrk);
    if(upstream_){
      double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
      // extrapolate through the IPA in this direction.
      bool exitsIPA = extrapolateIPA(ktrk,tdir);
      if(exitsIPA){ // if it exits out the back, extrapolate through the target
        bool exitsST = extrapolateST(ktrk,tdir);
        if(exitsST) { // if it exits out the back, extrapolate to the TSDA (DS rear absorber)
          bool hitTSDA = extrapolateTSDA(ktrk,tdir);
          // if we hit the TSDA we are done. Otherwise if we reflected, go back through the ST
          if( (!hitTSDA) ){ // reflection upstream of the target: go back through the target
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
        if(backToTracker_)extrapolateTracker(ktrk,tdir);
      }
      // optionally test for intersection with the OPA
      if(extrapolateOPA_)extrapolateOPA(ktrk,starttime,tdir);
    }
    // try CRV in both directions
    if(toCRV_){
      extrapolateDSMaterial(ktrk,TimeDir::backwards);
      extrapolatePassiveMaterial(ktrk,TimeDir::backwards);
      extrapolateCRV(ktrk,TimeDir::backwards);
      extrapolateDSMaterial(ktrk,TimeDir::forwards);
      extrapolatePassiveMaterial(ktrk,TimeDir::forwards);
      extrapolateCRV(ktrk,TimeDir::forwards);
    }
    // calorimeter disks MUST stay after all material passes: the calo extension adds no material, and running it
    // first moves DS shell crossings onto interior pieces where they are dropped (see PR #2059)
    // the calo is downstream of the tracker, opposite to the upstream direction tdir
    TimeDir calotdir = tdir == TimeDir::forwards ? TimeDir::backwards : TimeDir::forwards;
    // the front panels sit upstream of the disk faces, so their material must be added before the disk extension
    if(toCaloMaterial_)extrapolateCaloMaterial(ktrk,calotdir);
    if(toCaloD0_)toCaloDisk(ktrk,0,t0,calotdir);
    if(toCaloD1_)toCaloDisk(ktrk,1,t0,calotdir);
  }

  // Extrapolate a fit to the tracker perimeter (the closed cylinder bounding the active volume) and record
  // the bounding face it reaches: a z-end plane (TT_Front/TT_Back) for a downstream track, the outer wall
  // (TT_Outer) for a cosmic. R and the z-ends come from KinKalGeom::tracker(), i.e. from GeometryService.
  template <class KTRAJ> void KKExtrap::toTrackerPerimeter(KKTrack<KTRAJ>& ktrk) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& trk = *kkg_h->tracker();
    double const R  = trk.outer().radius();
    double const zf = trk.front().center().Z();
    double const zb = trk.back().center().Z();
    double const zmin = std::min(zf,zb); // upstream z-end plane (do not assume front<back)
    double const zmax = std::max(zf,zb); // downstream z-end plane
    ExtrapolateToTrackerPerimeter toPerim(maxdt_,maxdtstep_,btol_,R,zmin,zmax,debug_);
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId tt_outer("TT_Outer");
    static const SurfaceId tt_front("TT_Front");
    static const SurfaceId tt_back("TT_Back");
    static const SurfaceId tt_mid("TT_Mid");

    // middle reference, recorded once if the track crosses it
    auto midinter = KinKal::intersect(ftraj,trk.middle(),ftraj.range(),intertol_);
    if(midinter.good()) ktrk.addIntersection(tt_mid,midinter);

    // The track meets the perimeter at both ends; extend each time direction and keep the face crossed
    // first. Anchor the range and time origin on the tracker interior, and intersect the whole trajectory:
    // a calo-cluster track is already extrapolated to the calo, so its end piece lies beyond TT_Back.
    double tref = midinter.good() ? midinter.time_ : ftraj.t0();
    ftraj.range().forceRange(tref);
    auto recordFace = [&](TimeDir tdir){
      if(!ktrk.extrapolate(tdir,toPerim)) return;
      double const sgn = KinKal::timeDirSign(tdir);
      TimeRange const rng = tdir == TimeDir::forwards ?
        TimeRange(tref,ftraj.range().end()) : TimeRange(ftraj.range().begin(),tref);
      std::optional<std::pair<SurfaceId,KinKal::Intersection>> best;
      double bestdt = std::numeric_limits<double>::max();
      auto consider = [&](SurfaceId const& sid, auto const& surf){
        auto inter = KinKal::intersect(ftraj,surf,rng,intertol_,tdir);
        if(inter.good()){
          double const dt = sgn*(inter.time_ - tref); // time from the tracker interior out to this face
          if(dt >= 0.0 && dt < bestdt){ bestdt = dt; best.emplace(sid,inter); }
        }
      };
      consider(tt_outer, trk.outer());
      consider(tt_front, trk.front());
      consider(tt_back,  trk.back());
      if(best) ktrk.addIntersection(best->first,best->second);
    };
    recordFace(TimeDir::backwards);
    recordFace(TimeDir::forwards);
  }

  // Extrapolate a fit to a calorimeter disk (0 or 1) and record the intersections with its front/back faces
  // and inner/outer edges. tref is a time inside the tracker and tdir the time direction towards the calo, both
  // taken before any extrapolation
  template <class KTRAJ> void KKExtrap::toCaloDisk(KKTrack<KTRAJ>& ktrk, int idisk, double tref, TimeDir tdir) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    if(!kkg_h->calo()) return; // no calorimeter in this geometry
    auto const& calo = *kkg_h->calo();
    bool d0 = idisk == 0;
    auto const& front = d0 ? calo.EMC_Disk_0_Front() : calo.EMC_Disk_1_Front();
    auto const& back  = d0 ? calo.EMC_Disk_0_Back()  : calo.EMC_Disk_1_Back();
    auto const& inner = d0 ? calo.EMC_Disk_0_Inner() : calo.EMC_Disk_1_Inner();
    auto const& outer = d0 ? calo.EMC_Disk_0_Outer() : calo.EMC_Disk_1_Outer();
    static const SurfaceId frontsid[2] = {SurfaceId("EMC_Disk_0_Front"), SurfaceId("EMC_Disk_1_Front")};
    static const SurfaceId backsid[2]  = {SurfaceId("EMC_Disk_0_Back"),  SurfaceId("EMC_Disk_1_Back")};
    static const SurfaceId innersid[2] = {SurfaceId("EMC_Disk_0_Inner"), SurfaceId("EMC_Disk_1_Inner")};
    static const SurfaceId outersid[2] = {SurfaceId("EMC_Disk_0_Outer"), SurfaceId("EMC_Disk_1_Outer")};
    auto const& ftraj = ktrk.fitTraj();
    // extend to the back face; the front face and edges are then covered by the trajectory too
    ExtrapolateToZ caloBack(maxdt_,maxdtstep_,btol_,intertol_,back.center().Z(),debug_);
    if(debug_ > 0)std::cout << "toCaloDisk " << idisk << " front z " << front.center().Z() << " back z " << back.center().Z() << std::endl;
    if(!ktrk.extrapolate(tdir,caloBack)) return;
    // record the calo disk intersections, searching outward from the tracker interior only: after upstream
    // extrapolation the far end piece is a nominal-field helix whose range can reach back to the calo z
    TimeRange const rng = tdir == TimeDir::forwards ?
      TimeRange(tref,ftraj.range().end()) : TimeRange(ftraj.range().begin(),tref);
    auto frontinter = KinKal::intersect(ftraj,front,rng,intertol_,tdir);
    if(frontinter.good()) ktrk.addIntersection(frontsid[idisk],frontinter);
    auto backinter = KinKal::intersect(ftraj,back,rng,intertol_,tdir);
    if(backinter.good())ktrk.addIntersection(backsid[idisk],backinter);
    // get intersections with inner and outer cylinders
    auto innerinter = KinKal::intersect(ftraj,inner,rng,intertol_,tdir);
    if(innerinter.good()) ktrk.addIntersection(innersid[idisk],innerinter);
    auto outerinter = KinKal::intersect(ftraj,outer,rng,intertol_,tdir);
    if(outerinter.good()) ktrk.addIntersection(outersid[idisk],outerinter);
  }

  template <class KTRAJ> bool KKExtrap::extrapolateIPA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    GeomHandle<mu2e::KKMaterial> kkmat_h;
    using KKIPAXING = KKShellXing<KTRAJ,KinKal::Cylinder>;
    using KKIPAXINGPTR = std::shared_ptr<KKIPAXING>;
    // extraplate the fit through the IPA. This will add material effects for each intersection. It will continue till the
    // track exits the IPA
    ExtrapolateIPA extrapIPA(maxdt_,maxdtstep_,btol_,intertol_,kkg_h->DS()->innerProtonAbsorberPtr(),debug_);
    if(debug_ > 2)std::cout << "extrapolating to IPA " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId IPASID("IPA");
    double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    double startz = ftraj.position3(starttime).Z() - extrapIPA.zmid(); // z relative to the IPA center
    do {
      ktrk.extrapolate(tdir,extrapIPA);
      if(extrapIPA.intersection().good()){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        auto const& IPA = kkg_h->DS()->innerProtonAbsorberPtr();
        KKIPAXINGPTR ipaxingptr = std::make_shared<KKIPAXING>(IPA,IPASID,*kkmat_h->IPAMaterial(),extrapIPA.intersection(),reftrajptr,ipathick_,extrapIPA.interTolerance());
        if(debug_ > 2){
          double dmom, paramomvar, perpmomvar;
          ipaxingptr->materialEffects(dmom,paramomvar,perpmomvar);
          std::cout << "IPA Xing dmom " << dmom << " para momsig " << sqrt(paramomvar) << " perp momsig " << sqrt(perpmomvar) << std::endl;
          std::cout << " before append mom = " << reftrajptr->momentum();
        }
        ktrk.addIPAXing(ipaxingptr,tdir);
        if(debug_ > 2){
          auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
          std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
        }
      }
    } while(extrapIPA.intersection().good());
    // check if the particle exited the same side it entered (ie it reflected inside the IPA)
    double endtime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    double endz = ftraj.position3(endtime).Z() - extrapIPA.zmid(); // z relative to the IPA center
    if(endz * startz < 0.0){
      if(debug_ > 1)std::cout <<"Exiting opposite side IPA" << std::endl;
      return true;
    } else {
      if(debug_ > 1)std::cout <<"Exiting same side IPA" << std::endl;
      return false;
    }
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
    double startz = ftraj.position3(starttime).Z() - extrapST.zmid();
    if(debug_ > 2)std::cout << "extrapolating to ST " << std::endl;
    do {
      ktrk.extrapolate(tdir,extrapST);
      if(extrapST.intersection().good()){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        KKSTXINGPTR stxingptr = std::make_shared<KKSTXING>(extrapST.foil(),extrapST.foilId(),*kkmat_h->STMaterial(),extrapST.intersection(),reftrajptr,stthick_,extrapST.interTolerance());
        if(debug_ > 2){
          double dmom, paramomvar, perpmomvar;
          stxingptr->materialEffects(dmom,paramomvar,perpmomvar);
          std::cout << "ST Xing dmom " << dmom << " para momsig " << sqrt(paramomvar) << " perp momsig " << sqrt(perpmomvar) << std::endl;
          std::cout << " before append mom = " << reftrajptr->momentum();
        }
        // Add the xing. This truncates the fit
        ktrk.addSTXing(stxingptr,tdir);
        if(debug_ > 2){
          auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
          std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
        }
      }
    } while(extrapST.intersection().good());
    // check if the particle exited in the same physical direction or not (reflection)
    double endtime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    double endz = ftraj.position3(endtime).Z() - extrapST.zmid(); // z relative to the ST center
    if(endz * startz < 0.0){
      if(debug_ > 1)std::cout <<"Exiting opposite side ST" << std::endl;
      return true;
    } else {
      if(debug_ > 1)std::cout <<"Exiting same side ST" << std::endl;
      return false;
    }
  }

  template <class KTRAJ> bool KKExtrap::extrapolateTracker(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    ExtrapolateToZ trackerFront(maxdt_,maxdtstep_,btol_,intertol_,kkg_h->tracker()->front().center().Z(),debug_);
    if(debug_ > 2)std::cout << "extrapolating to Tracker " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    double tstart = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    static const SurfaceId TrackerSID("TT_Front");
    ktrk.extrapolate(tdir,trackerFront);
    TimeRange testrange = tdir == TimeDir::forwards ? TimeRange(tstart, ftraj.range().end())
      : TimeRange(ftraj.range().begin(),tstart);
    auto trkfrontinter = KinKal::intersect(ftraj,kkg_h->tracker()->front(),testrange,intertol_,tdir);
    if(trkfrontinter.onsurface_){ // dont worry about bounds here
      ktrk.addIntersection(TrackerSID,trkfrontinter);
      return true;
    }
    return false;
  }

  template <class KTRAJ> bool KKExtrap::extrapolateTSDA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    ExtrapolateToZ TSDA(maxdt_,maxdtstep_,btol_,intertol_,kkg_h->DS()->upstreamAbsorber().center().Z(),debug_);
    if(debug_ > 2)std::cout << "extrapolating to TSDA " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    double tstart = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    static const SurfaceId TSDASID("TSDA");
    ktrk.extrapolate(tdir,TSDA);
    TimeRange testrange = tdir == TimeDir::forwards ? TimeRange(tstart, ftraj.range().end())
      : TimeRange(ftraj.range().begin(),tstart);
    // test for an intersection with TSDA downstream surface
    auto tsdainter = KinKal::intersect(ftraj,kkg_h->DS()->upstreamAbsorber(),testrange,intertol_,tdir);
    if(tsdainter.onsurface_){
      ktrk.addIntersection(TSDASID,tsdainter);
      return true;
    }
    return false;
  }

  //extrapolate to the OPA (DS entrance). only 1 of these is currently possible
  // muons can intersect multiple times: TODO
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
  template <class KTRAJ> void KKExtrap::extrapolateDSMaterial(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    using KKMATCYLXING = KKShellXing<KTRAJ,KinKal::Cylinder>;
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    GeomHandle<mu2e::KKMaterial> kkmat_h;
    auto const& cylinders = kkg_h->DS()->materialCylinders();
    if(cylinders.empty()) return;
    auto extrapCylinders = ExtrapolateCylinders(maxdt_,maxdtstep_,btol_,intertol_,minv_,cylinders,debug_);
    auto const& ftraj = ktrk.fitTraj();
    using MaterialCylinder = ExtrapolateCylinders::MaterialCylinder;
    // NB: this bank-all-sorted idiom differs deliberately from the one-crossing-at-a-time idiom in
    // extrapolatePassiveMaterial / extrapolateCRV below. It is safe HERE because the DS shells are
    // CONCENTRIC and closely spaced: one step overshoots the whole cluster and the resulting trajectory
    // piece spans all of them, so all crossings can be added in a single sorted pass. The planar surfaces
    // are STACKED and spatially separated, so adding more than the nearest crossing before re-extrapolating
    // advances the front past the farther planes and strands them -- hence they must go one at a time.
    std::unordered_set<MaterialCylinder const*> banked; // bank each cylinder once (avoid re-found duplicates)
    do {
      ktrk.extrapolate(tdir,extrapCylinders);
      if(debug_ > 0) std::cout << "Found " << extrapCylinders.intersections().size() << " DS material intersections " << std::endl;
      // Bank every crossing found in this pass, in strict tracker-inward time order (descending time
      // for backward, ascending for forward). The DS material shells are concentric and closely spaced,
      // so a single extrapolation step overshoots the whole cluster and the trajectory piece spans all
      // of them; adding the crossings in this order makes each prepend/append shrink the piece down to
      // exactly the next crossing, so every one stays within the addable front/back range and none is
      // stranded as interior. (intersections() is not reliably time-sorted, so sort a local copy.) The
      // 'banked' set skips a crossing re-found on a later step; once a pass adds nothing new, we are done.
      auto inters = extrapCylinders.intersections();
      std::sort(inters.begin(),inters.end(),[tdir](auto const& a,auto const& b){
        return tdir == TimeDir::forwards ? a.inter_.time_ < b.inter_.time_ : a.inter_.time_ > b.inter_.time_; });
      bool newbank = false;
      for(auto const& inter : inters) {
        if(!banked.insert(inter.cylinder_).second) continue;
        newbank = true;
        auto const& cylinder = *inter.cylinder_;
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        auto matxingptr = std::make_shared<KKMATCYLXING>(cylinder.surface_,cylinder.sid_,*kkmat_h->material(cylinder.material_),
            inter.inter_,reftrajptr,cylinder.thickness_,extrapCylinders.interTolerance());
        ktrk.addMaterialCylXing(matxingptr,tdir);
      }
      if(!newbank) break;
    } while(extrapCylinders.intersections().size()>0);
  }

  template <class KTRAJ> void KKExtrap::extrapolatePassiveMaterial(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    using KKMATRECXING = KKShellXing<KTRAJ,KinKal::Plane>;  // passive plane (Rectangle or Annulus)
    using MaterialPlane = ExtrapolatePlanes::MaterialPlane;
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    GeomHandle<mu2e::KKMaterial> kkmat_h;
    auto const& planes = kkg_h->passiveMaterialPlanes();
    if(planes.empty()) return;
    auto extrapPlanes = ExtrapolatePlanes(maxdt_,maxdtstep_,btol_,intertol_,minv_,planes,debug_);
    auto const& ftraj = ktrk.fitTraj();

    // Normal path: process each plane at most once, ONE crossing per extrapolate pass (add nearest,
    // re-extrapolate, repeat). Unlike extrapolateDSMaterial's bank-all-sorted pass, the passive planes
    // are stacked and spatially separated, so adding more than the nearest crossing before re-extrapolating
    // would advance the front past the farther planes and strand them (see the cross-reference note there).
    // Stop before calling ktrk.extrapolate when all planes have been processed so
    // we don't advance the trajectory past the last plane (which would overshoot the
    // CRV and prevent extrapolateCRV from finding it).
    std::unordered_set<MaterialPlane const*> processed;
    // The CRV strongback planes are crossed together with the CRV sectors in extrapolateCRV
    // (interleaved, time-ordered): processing all strongbacks HERE advances the front past the
    // lower sectors of a stacked geometry and strands them. Exclude the strongbacks from this
    // (concrete-only for the extracted geometry) pass by pre-marking them processed.
    SurfaceId const sbid(SurfaceIdEnum::CRV_StrongBack);
    for(auto const& plane : planes) if(plane.sid_.id() == sbid.id()) processed.insert(&plane);
    do {
      if(processed.size() >= planes.size()) break;
      ktrk.extrapolate(tdir,extrapPlanes);
      if(debug_ > 0) std::cout << "Found " << extrapPlanes.intersections().size() << " passive material plane intersections " << std::endl;
      bool found = false;
      for(auto const& inter : extrapPlanes.intersections()) {
        if(processed.count(inter.plane_) > 0) continue;
        auto const& plane = *inter.plane_;
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        auto matxingptr = std::make_shared<KKMATRECXING>(plane.surface_,plane.sid_,*kkmat_h->material(plane.material_),
            inter.inter_,reftrajptr,plane.thickness_,extrapPlanes.interTolerance());
        ktrk.addMaterialPlaneXing(matxingptr,tdir);
        processed.insert(inter.plane_);
        found = true;
        break;
      }
      if(!found) break;
    } while(extrapPlanes.intersections().size()>0);
  }

  // Extrapolate through the calorimeter front panel layers and their cooling pipes (both disks), adding their material
  // effects. All are described by planes normal to z, so visit them in z order along the extrapolation direction: extend
  // the fit to each plane's z, then add a crossing if the track intersects the plane within the extension (ie it doesn't
  // pass through the hole) and, for the pipe planes, if it crosses a pipe wall there.
  // The search starts at the last measurement, not at the end of the fit: the fit extends past its last measurement,
  // possibly into the front panel, and its end can move on refit, so crossings between the two would otherwise be lost.
  // Crossings can only be added in the last trajectory piece, so one before it is skipped: crossings inside the fit
  // are added by KKFit AddCaloMaterial, which can refit.
  // Layers already in the fit are skipped, and if an earlier pass (ie ToCRV) took the fit past the calorimeter no
  // front panel material is added.
  // Tracks moving downstream then cross the crystal case inner support rings of a disk after its front panel, see
  // extrapolateCaloSupport.
  template <class KTRAJ> void KKExtrap::extrapolateCaloMaterial(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    using KKMATRECXING = KKShellXing<KTRAJ,KinKal::Plane>;
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    GeomHandle<mu2e::KKMaterial> kkmat_h;
    auto layers = caloMaterialLayers(*kkg_h);
    if(layers.empty() && kkg_h->caloInnerSupports().empty()) return;
    auto const& ftraj = ktrk.fitTraj();
    auto arange = ktrk.activeRange();
    // a measurement (ie the CaloHit) can be beyond the end of the fit
    double tstart = std::clamp(tdir == TimeDir::forwards ? arange.end() : arange.begin(),ftraj.range().begin(),ftraj.range().end());
    // order the layers along the direction of motion in z at the start of the extrapolation
    double zdir = ftraj.direction(tstart).Z()*KinKal::timeDirSign(tdir);
    std::sort(layers.begin(),layers.end(),[zdir](auto const& a,auto const& b){
        return zdir > 0 ? a.surface_->center().Z() < b.surface_->center().Z() : a.surface_->center().Z() > b.surface_->center().Z(); });
    std::vector<KKGeom::CaloInnerSupport const*> supports;
    if(tdir == TimeDir::forwards && zdir > 0.0)
      for(auto const& support : kkg_h->caloInnerSupports())supports.push_back(&support);
    std::sort(supports.begin(),supports.end(),[](auto a,auto b){ return a->front()->center().Z() < b->front()->center().Z(); });
    auto nextsupport = supports.begin();
    bool usable(true);
    for(auto const& layer : layers) {
      while(usable && nextsupport != supports.end() && (*nextsupport)->front()->center().Z() < layer.surface_->center().Z())
        usable = extrapolateCaloSupport(ktrk,**nextsupport++,tstart);
      if(!usable) break;
      bool infit(false);
      for(auto const& matxing : ktrk.materialPlaneXings())
        if(matxing->surfaceId() == layer.sid_){ infit = true; break; }
      if(infit){
        if(debug_ > 0)std::cout << "Calo front panel crossing " << layer.sid_ << " already in the fit" << std::endl;
        continue;
      }
      ExtrapolateToZ toLayer(maxdt_,maxdtstep_,btol_,intertol_,layer.surface_->center().Z(),debug_);
      if(!ktrk.extrapolate(tdir,toLayer)){ usable = false; break; } // the fit is no longer usable
      TimeRange testrange = tdir == TimeDir::forwards ? TimeRange(tstart,ftraj.range().end()) : TimeRange(ftraj.range().begin(),tstart);
      auto inter = KinKal::intersect(ftraj,*layer.surface_,testrange,intertol_,tdir);
      if(debug_ > 1)std::cout << "Calo front panel " << layer.sid_ << " " << inter << std::endl;
      if(!inter.good()) continue;
      double thickness = layer.thickness(inter.pos_);
      if(thickness <= 0.0) continue; // between the pipes
      auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
      // a crossing can only be added at the end of the fit; one before it is skipped
      if(!reftrajptr->range().inRange(inter.time_)){
        if(debug_ > 0)std::cout << "Skipping calo front panel crossing " << layer.sid_ << " not at the end of the fit" << std::endl;
        continue;
      }
      auto matxingptr = std::make_shared<KKMATRECXING>(layer.surface_,layer.sid_,*kkmat_h->material(*layer.material_),
          inter,reftrajptr,thickness,intertol_);
      if(debug_ > 2){
        double dmom, paramomvar, perpmomvar;
        matxingptr->materialEffects(dmom,paramomvar,perpmomvar);
        std::cout << layer.sid_ << " thickness " << thickness << " Xing dmom " << dmom << " para momsig " << sqrt(paramomvar)
          << " perp momsig " << sqrt(perpmomvar) << " before append mom = " << reftrajptr->momentum();
      }
      ktrk.addMaterialPlaneXing(matxingptr,tdir);
      if(debug_ > 2){
        auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
      }
      // the next layer must be crossed after this one
      tstart = inter.time_;
    }
    while(usable && nextsupport != supports.end())
      usable = extrapolateCaloSupport(ktrk,**nextsupport++,tstart);
  }

  // Extrapolate forwards through the inner support rings of a calorimeter crystal case, see caloSupportXings: extend the
  // fit to the case front and, if the track is then in the disk hole, to the back of the case, and add the crossings
  // found after tstart. As for the front panel, crossings can only be added in the last trajectory piece, and rings
  // already in the fit are skipped. Returns false if the fit is no longer usable.
  template <class KTRAJ> bool KKExtrap::extrapolateCaloSupport(KKTrack<KTRAJ>& ktrk,KKGeom::CaloInnerSupport const& support,double& tstart) const {
    using Support = KKGeom::CaloInnerSupport;
    using KKMATRECXING = KKShellXing<KTRAJ,KinKal::Plane>;
    using KKMATCYLXING = KKShellXing<KTRAJ,KinKal::Cylinder>;
    static const TimeDir tdir(TimeDir::forwards);
    GeomHandle<mu2e::KKMaterial> kkmat_h;
    auto ringid = support.surfaceId(Support::alFront).id();
    auto isring = [&ringid](auto const& matxing){ return matxing->surfaceId().id() == ringid; };
    if(std::any_of(ktrk.materialPlaneXings().begin(),ktrk.materialPlaneXings().end(),isring) ||
        std::any_of(ktrk.materialCylXings().begin(),ktrk.materialCylXings().end(),isring)){
      if(debug_ > 0)std::cout << "Calo inner ring crossing " << support.surfaceId(Support::alFront) << " already in the fit" << std::endl;
      return true;
    }
    auto const& ftraj = ktrk.fitTraj();
    double zfront = support.front()->center().Z();
    ExtrapolateToZ toFront(maxdt_,maxdtstep_,btol_,intertol_,zfront,debug_);
    if(!ktrk.extrapolate(tdir,toFront)) return false;
    if(support.transverseRadius(ftraj.position3(ftraj.range().end())) < support.rin()){
      ExtrapolateToZ toBack(maxdt_,maxdtstep_,btol_,intertol_,zfront + 2.0*support.inner()->halfLength(),debug_);
      if(!ktrk.extrapolate(tdir,toBack)) return false;
    }
    for(auto const& sxing : caloSupportXings(ftraj,support,TimeRange(tstart,ftraj.range().end()),intertol_)){
      auto sid = support.surfaceId(sxing.xing_);
      auto const& reftrajptr = ftraj.backPtr();
      if(!reftrajptr->range().inRange(sxing.inter_.time_)){
        if(debug_ > 0)std::cout << "Skipping calo inner ring crossing " << sid << " not at the end of the fit" << std::endl;
        continue;
      }
      auto const& dmat = *kkmat_h->material(support.materialName(Support::material(sxing.xing_)));
      if(debug_ > 1)std::cout << "Calo inner ring crossing " << sid << " path " << sxing.path_ << " " << sxing.inter_ << std::endl;
      if(sxing.xing_ == Support::alFront || sxing.xing_ == Support::cfFront)
        ktrk.addMaterialPlaneXing(std::make_shared<KKMATRECXING>(support.front(),sid,dmat,sxing.inter_,reftrajptr,sxing.thickness(),intertol_),tdir);
      else
        ktrk.addMaterialCylXing(std::make_shared<KKMATCYLXING>(sxing.xing_ == Support::alHole ? support.inner() : support.cfInner(),
              sid,dmat,sxing.inter_,reftrajptr,sxing.thickness(),intertol_),tdir);
      tstart = sxing.inter_.time_;
    }
    return true;
  }

  template <class KTRAJ> void KKExtrap::extrapolateCRV(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const {
    using KKCRVXING = KKShellXing<KTRAJ,KinKal::Rectangle>;
    using KKMATRECXING = KKShellXing<KTRAJ,KinKal::Plane>;  // strongback (passive) plane
    using MaterialPlane = ExtrapolateCRVRegion::MaterialPlane;
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    GeomHandle<mu2e::KKMaterial> kkmat_h;
    // The CRV sector planes and their aluminium strongbacks are spatially interleaved (the sector plane
    // is modeled at the scintillator-stack center; its tracker-facing strongback plane sits ~60 mm --
    // a stack half-thickness -- toward the tracker), so they MUST be crossed in ONE outward,
    // time-ordered pass -- a xing can only be added at the trajectory front/back, never the interior.
    // Crossing all strongbacks first (the old extrapolatePassiveMaterial pass) advanced the front past
    // the lower sectors of a stacked geometry (extracted EX/T1/T2 ~150 mm apart in y) and stranded them
    // (CRV_T1 lost on every track). Hand both surface types to ExtrapolateCRVRegion, which stops at the
    // next crossing of either; add them one at a time, innermost first, re-extrapolating between each
    // so every crossing is added at the current front (mirrors extrapolatePassiveMaterial). Run-2 sees
    // <=1 sector per direction, so the order (strongback then sector) is unchanged.
    ExtrapolateCRVRegion::MaterialPlaneCollection strongbacks;
    SurfaceId const sbid(SurfaceIdEnum::CRV_StrongBack);
    for(auto const& plane : kkg_h->passiveMaterialPlanes())
      if(plane.sid_.id() == sbid.id()) strongbacks.push_back(plane);
    auto extrapCRV = ExtrapolateCRVRegion(maxdt_,maxdtstep_,btol_,intertol_,minv_,*kkg_h->CRV(),strongbacks,debug_);
    if(debug_ > 5){std::cout << "Extrapolating to CRV with " << extrapCRV.sectors().size() << " sectors + "
        << strongbacks.size() << " strongbacks" << std::endl;
      for(auto const& sector : kkg_h->CRV()->sectors())
        std::cout << sector.sname_ << " position " << sector.sector_->center() << " halfwidth " << sector.whw_ << std::endl;
    }
    auto const& ftraj = ktrk.fitTraj();
    std::unordered_set<int> processedSectors;
    std::unordered_set<MaterialPlane const*> processedPlanes;
    do {
      ktrk.extrapolate(tdir,extrapCRV);
      if(debug_ > 0) std::cout << "Found " << extrapCRV.intersections().size() << " CRV-region intersections " << std::endl;
      bool added = false;
      // intersections are time-sorted (innermost first for this tdir): add the nearest un-processed
      // one, then re-extrapolate so the next is again at the front.
      for(auto const& xing : extrapCRV.intersections()){
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        if(xing.isSector_){
          if(processedSectors.count(xing.sectorIdx_) > 0) continue;
          auto crvxingptr = std::make_shared<KKCRVXING>(extrapCRV.sector((size_t)xing.sectorIdx_).sector_,
              SurfaceId(kkg_h->CRV()->sectorName(xing.sectorIdx_)),*kkmat_h->CRVMaterial(),xing.inter_,reftrajptr,
              2*xing.whw_, extrapCRV.interTolerance());
          ktrk.addCRVXing(crvxingptr,tdir);
          processedSectors.insert(xing.sectorIdx_);
          if(debug_ > 1) std::cout << "Good CRV sector " << xing.sectorIdx_ << " " << xing.inter_ << std::endl;
        } else {
          if(processedPlanes.count(xing.plane_) > 0) continue;
          auto const& plane = *xing.plane_;
          auto matxingptr = std::make_shared<KKMATRECXING>(plane.surface_,plane.sid_,*kkmat_h->material(plane.material_),
              xing.inter_,reftrajptr,plane.thickness_,extrapCRV.interTolerance());
          ktrk.addMaterialPlaneXing(matxingptr,tdir);
          processedPlanes.insert(xing.plane_);
          if(debug_ > 1) std::cout << "Good CRV strongback " << plane.sid_ << " " << xing.inter_ << std::endl;
        }
        added = true;
        break;
      }
      if(!added) break;
    } while(extrapCRV.intersections().size()>0);
  }

}
#endif
