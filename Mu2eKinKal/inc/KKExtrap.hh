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
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
#include "Offline/KinKalGeom/inc/KKMaterial.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
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
      // extrapolation functions; these are templated on the type of trajectory
      template <class KTRAJ> void extrapolate(KKTrack<KTRAJ>& ktrk) const;
      template <class KTRAJ> bool extrapolateIPA(KKTrack<KTRAJ>& ktrk,TimeDir trkdir) const;
      template <class KTRAJ> bool extrapolateST(KKTrack<KTRAJ>& ktrk,TimeDir trkdir) const;
      template <class KTRAJ> bool extrapolateTracker(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> bool extrapolateTSDA(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void extrapolateDSMaterial(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void extrapolatePassiveMaterial(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void extrapolateCRV(KKTrack<KTRAJ>& ktrk,TimeDir tdir) const;
      template <class KTRAJ> void extrapolateOPA(KKTrack<KTRAJ>& ktrk, double tstart, TimeDir tdir) const;
      template <class KTRAJ> void toTrackerPerimeter(KKTrack<KTRAJ>& ktrk) const;

    private:
      int debug_;
      double btol_, intertol_, maxdt_, maxdtstep_, minv_;
      bool backToTracker_, extrapolateOPA_, toTrackerPerimeter_, upstream_, toCRV_;
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
    toCRV_(extrapconfig.ToCRV())
  {}

  template <class KTRAJ> void KKExtrap::extrapolate(KKTrack<KTRAJ>& ktrk) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    // define the time direction according to the fit direction inside the tracker
    auto const& ftraj = ktrk.fitTraj();
    auto dir0 = ftraj.direction(ftraj.t0());
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
