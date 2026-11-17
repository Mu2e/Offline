#ifndef Mu2eKinKal_KKExtrap_hh
#define Mu2eKinKal_KKExtrap_hh
// Helper class for time-ordered extrapolation through requested materials and surfaces.
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/KinKalGeom/inc/Tracker.hh"
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateSweep.hh"
#include "Offline/Mu2eKinKal/inc/CaloMaterialLayers.hh"
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
#include "Offline/KinKalGeom/inc/KKMaterial.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "cetlib_except/exception.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include <utility>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
namespace mu2e {
  using KinKal::VEC3;
  using KinKal::TimeDir;
  using Mu2eKinKal::KKExtrapConfig;

  class KKExtrap {
    public:
      explicit KKExtrap(KKExtrapConfig const& exconfig);
      // throw at configuration if the fit config's minimum domain step can't fit this extrapolation's maximum step
      void checkConfig(KinKal::Config const& config) const;
      // extrapolate the fit in both time directions through every requested surface
      template <class KTRAJ> void extrapolate(KKTrack<KTRAJ>& ktrk) const;
      // the surfaces requested in one time direction; updir is the direction that leaves the tracker upstream
      // To add passive material: give it a SurfaceId and a KinKalGeom surface, add a SweepTarget here (effect_, surf_ and the typed surface, mat_, thick_, reach_, mincos_),
      // and, for a new Effect, a case in addCrossing that builds its KinKal material xing; a surface that is neither a Plane nor a Cylinder is searched with KinKal::intersect
      // and is never cut by withinReach; add a case to ExtrapolateSweepTest and, if optional, a KKExtrapConfig switch.
      SweepTargetCollection targets(bool updir) const;
    private:
      // where a track entered a calorimeter crystal case inner support: the material still to cross at the carbon fiber ring inner radius
      struct SupportEntry { double path_ = 0.0; KKGeom::CaloInnerSupport::Crossing xing_ = KKGeom::CaloInnerSupport::cfRing; double time_ = 0.0; };
      using SupportEntries = std::map<KKGeom::CaloInnerSupport const*,SupportEntry>;
      template <class KTRAJ> void sweep(KKTrack<KTRAJ>& ktrk, SweepTargetCollection const& tgts, TimeDir tdir, double t0, bool updir) const;
      template <class KTRAJ> void addFailure(KKTrack<KTRAJ>& ktrk, TimeDir tdir, std::string const& reason) const;
      template <class KTRAJ> bool rangeOut(KKTrack<KTRAJ>& ktrk, SweepTarget const& tg, KinKal::Intersection const& inter) const;
      template <class KTRAJ> bool addCrossing(KKTrack<KTRAJ>& ktrk, SweepTarget const& tg, KinKal::Intersection const& inter, TimeDir tdir) const;
      template <class KTRAJ> bool supportCrossing(KKTrack<KTRAJ>& ktrk, ExtrapolateSweep& xs, size_t itg, KinKal::Intersection const& inter,
          TimeDir tdir, SupportEntries& entries) const;
      template <class KTRAJ> bool reflect(ExtrapolateSweep& xs, KinKal::ParticleTrajectory<KTRAJ> const& ftraj, TimeDir tdir) const;
      int debug_;
      double btol_, intertol_, maxdt_, maxdtstep_, minv_;
      bool backToTracker_, extrapolateOPA_, toTrackerPerimeter_, upstream_, toCRV_;
      bool toCaloD0_, toCaloD1_, toCaloMaterial_;
      // nominal literal like the hard-coded KinKalGeom IPA surface (KinKalGeomMaker::makeDS); TODO: take it from the geometry service.
      // The ST foil thicknesses come from the geometry service (KKGeom::StoppingTarget::foilThickness)
      double ipathick_ = 0.511; // ipa thickness (mm)
      static constexpr unsigned maxCrossings_ = 1000; // per sweep, against endless reflection loops
      static constexpr double dsShellMinCos_ = 0.1; // thin-shell approximation: skip near-tangent DS shell crossings
      static constexpr double ipaMinCos_ = 0.1; // as the DS shells; no truth IPA crossing on CE has |cos| < 0.2
      static constexpr double stFoilMinCos_ = 0.02; // no truth foil crossing on CE has |cos| < 0.02; the thin-foil path is then 50x the thickness
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
    auto const& ftraj = ktrk.fitTraj();
    if(ftraj.range().null()) return;
    // Clamp the anchor because the helix time at z=0 need not lie inside the fitted range.
    double t0 = ftraj.t0();
    ftraj.range().forceRange(t0);
    TimeDir updir = (ftraj.direction(t0).Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
    if(toTrackerPerimeter_){
      GeomHandle<mu2e::KinKalGeom> kkg_h;
      static const SurfaceId tt_mid("TT_Mid");
      auto midinter = KinKal::intersect(ftraj,kkg_h->tracker()->middle(),ftraj.range(),intertol_);
      if(midinter.good()) ktrk.addIntersection(tt_mid,midinter);
    }
    for(auto tdir : {TimeDir::backwards, TimeDir::forwards}){
      auto tgts = targets(updir == tdir); // geometry or service failures stop the job
      try {
        sweep(ktrk,tgts,tdir,t0,updir == tdir);
      } catch(std::bad_alloc const&) {
        throw;
      } catch(std::exception const& error) {
        // a problem with this track ends this direction only; the crossings already added are kept
        addFailure(ktrk,tdir,error.what());
      }
    }
  }

  template <class KTRAJ> void KKExtrap::addFailure(KKTrack<KTRAJ>& ktrk, TimeDir tdir, std::string const& reason) const {
    std::string failure = std::string(tdir == TimeDir::forwards ? "forwards: " : "backwards: ") + reason;
    if(debug_ > 0) std::cout << "KKExtrap: extrapolation failed " << failure << std::endl;
    ktrk.addExtrapolationFailure(failure);
  }

  // The particle stops in this material: the forwards direction ends before it, keeping the crossings already added.
  template <class KTRAJ> bool KKExtrap::rangeOut(KKTrack<KTRAJ>& ktrk, SweepTarget const& tg, KinKal::Intersection const& inter) const {
    if(debug_ > 0) std::cout << "KKExtrap: ranged out in " << tg.sid_ << " at " << inter.time_ << std::endl;
    ktrk.setRangedOut();
    return false;
  }

  inline SweepTargetCollection KKExtrap::targets(bool updir) const {
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    GeomHandle<mu2e::KKMaterial> kkmat_h;
    SweepTargetCollection tgts;
    auto const& ds = *kkg_h->DS();
    double const bore = ds.inner().radius(); // surfaces inside the DS bore cannot be reached once a track leaves it
    auto zrange = [](KinKal::Cylinder const& cyl, double& zlo, double& zhi){
      double z0 = (cyl.center() - cyl.axis()*cyl.halfLength()).Z();
      double z1 = (cyl.center() + cyl.axis()*cyl.halfLength()).Z();
      zlo = std::min(z0,z1); zhi = std::max(z0,z1);
    };
    // tracker perimeter: the first face crossed, in each direction
    if(toTrackerPerimeter_){
      auto const& trk = *kkg_h->tracker();
      auto face = [&](char const* name, KinKal::Surface const* surf){
        SweepTarget tg;
        tg.sid_ = SurfaceId(name);
        tg.surf_ = surf;
        tg.reach_ = SweepTarget::inside;
        tg.group_ = 0;
        tg.rreg_ = trk.outer().radius();
        tg.zreglo_ = std::min(trk.front().center().Z(),trk.back().center().Z());
        tg.zreghi_ = std::max(trk.front().center().Z(),trk.back().center().Z());
        return tg;
      };
      auto outer = face("TT_Outer",trk.outerPtr().get());
      tgts.push_back(outer);
      auto front = face("TT_Front",trk.frontPtr().get());
      tgts.push_back(front);
      auto back = face("TT_Back",trk.backPtr().get());
      tgts.push_back(back);
    }
    // upstream: IPA, ST foils, TSDA (ends the sweep), OPA, and optionally the tracker entrance after a reflection
    if(updir && upstream_){
      auto const& ipa = ds.innerProtonAbsorberPtr();
      SweepTarget ipatg;
      ipatg.effect_ = SweepTarget::ipashell;
      ipatg.reach_ = SweepTarget::zrange;
      ipatg.sid_ = SurfaceId("IPA");
      ipatg.surf_ = ipa.get();
      ipatg.cyl_ = ipa;
      ipatg.mat_ = &*kkmat_h->IPAMaterial();
      ipatg.thick_ = ipathick_;
      ipatg.maxcross_ = 0;
      ipatg.mincos_ = std::max(minv_,ipaMinCos_);
      zrange(*ipa,ipatg.zreglo_,ipatg.zreghi_);
      tgts.push_back(ipatg);
      auto const& st = *kkg_h->ST();
      double stzlo, stzhi;
      zrange(st.outer(),stzlo,stzhi);
      for(size_t ifoil = 0; ifoil < st.foils().size(); ++ifoil){
        auto const& foil = st.foils()[ifoil];
        SweepTarget tg;
        tg.effect_ = SweepTarget::stfoil;
        tg.reach_ = SweepTarget::zrange;
        tg.sid_ = SurfaceId(SurfaceIdEnum::ST_Foils,ifoil);
        tg.surf_ = foil.get();
        tg.pl_ = foil.get(); // a grazing foil is crossed again only after the track leaves its tolerance
        tg.ann_ = foil;
        tg.mat_ = &*kkmat_h->STMaterial();
        tg.thick_ = st.foilThickness(ifoil);
        tg.maxcross_ = 0;
        tg.mincos_ = std::max(minv_,stFoilMinCos_);
        tg.zreglo_ = stzlo; tg.zreghi_ = stzhi;
        tgts.push_back(tg);
      }
      SweepTarget tsda;
      tsda.reach_ = SweepTarget::zplane;
      tsda.sid_ = SurfaceId("TSDA");
      tsda.surf_ = ds.upstreamAbsorberPtr().get();
      tsda.zval_ = ds.upstreamAbsorber().center().Z();
      tsda.terminal_ = true;
      tsda.inbounds_ = false;
      tgts.push_back(tsda);
      if(extrapolateOPA_){
        SweepTarget opa;
        opa.sid_ = SurfaceId("OPA");
        opa.surf_ = ds.outerProtonAbsorberPtr().get();
        opa.drives_ = false;
        tgts.push_back(opa);
      }
      if(backToTracker_){
        auto const& trk = *kkg_h->tracker();
        SweepTarget entr;
        entr.reach_ = SweepTarget::zplane;
        entr.sid_ = SurfaceId("TT_Front");
        entr.surf_ = trk.frontPtr().get();
        entr.zval_ = trk.front().center().Z();
        entr.terminal_ = true;
        entr.inbounds_ = false;
        entr.active_ = false;
        tgts.push_back(entr);
      }
    }
    // DS material shells, passive planes (concrete, DS end walls, CRV strongbacks) and CRV sectors, both directions
    if(toCRV_){
      for(auto const& cyl : ds.materialCylinders()){
        SweepTarget tg;
        tg.effect_ = SweepTarget::cylshell;
        tg.reach_ = SweepTarget::radial;
        tg.sid_ = cyl.sid_;
        tg.surf_ = cyl.surface_.get();
        tg.cyl_ = cyl.surface_;
        tg.mat_ = &*kkmat_h->material(cyl.material_);
        tg.thick_ = cyl.thickness_;
        tg.mincos_ = std::max(minv_,dsShellMinCos_);
        tg.rmax_ = cyl.surface_->radius() + 0.5*cyl.thickness_;
        tgts.push_back(tg);
      }
      auto plane = [&](SweepTarget& tg, KinKal::Plane const& pl){
        tg.reach_ = SweepTarget::plane;
        tg.pl_ = &pl;
      };
      for(auto const& pmp : kkg_h->passiveMaterialPlanes()){
        SweepTarget tg;
        tg.effect_ = SweepTarget::planeshell;
        tg.sid_ = pmp.sid_;
        tg.surf_ = pmp.surface_.get();
        tg.plane_ = pmp.surface_;
        tg.mat_ = &*kkmat_h->material(pmp.material_);
        tg.thick_ = pmp.thickness_;
        plane(tg,*pmp.surface_);
        tgts.push_back(tg);
      }
      if(kkg_h->CRV()){
        for(auto const& sector : kkg_h->CRV()->sectors()){
          SweepTarget tg;
          tg.effect_ = SweepTarget::crvsector;
          tg.sid_ = SurfaceId(sector.sname_);
          tg.surf_ = sector.sector_.get();
          tg.rect_ = sector.sector_;
          tg.mat_ = &*kkmat_h->CRVMaterial();
          tg.thick_ = 2*sector.whw_;
          tg.crvsector_ = &sector;
          plane(tg,*sector.sector_);
          tgts.push_back(tg);
        }
      }
    }
    // calorimeter material, downstream only: front panel layers and cooling pipes, and the crystal case inner supports
    if(!updir && toCaloMaterial_){
      for(auto const& layer : caloMaterialLayers(*kkg_h)){
        SweepTarget tg;
        tg.effect_ = SweepTarget::calopanel;
        tg.reach_ = SweepTarget::plane;
        tg.sid_ = layer.sid_;
        tg.surf_ = layer.surface_.get();
        tg.pl_ = layer.surface_.get();
        tg.plane_ = layer.surface_;
        tg.mat_ = &*kkmat_h->material(*layer.material_);
        tg.thick_ = layer.thickness_;
        tg.pipes_ = layer.pipes_;
        tg.rmax_ = bore;
        tgts.push_back(tg);
      }
      // the case front finds where the case is entered; the hole and carbon fiber ring parts are switched on by supportCrossing
      for(auto const& support : kkg_h->caloInnerSupports()){
        double zback = support.front()->center().Z() + 2.0*support.inner()->halfLength();
        for(auto part : {SweepTarget::casefront, SweepTarget::casehole, SweepTarget::cfring}){
          SweepTarget tg;
          tg.effect_ = SweepTarget::calosupport;
          tg.support_ = &support;
          tg.part_ = part;
          tg.sid_ = support.surfaceId(KKGeom::CaloInnerSupport::alFront);
          tg.rmax_ = bore;
          if(part == SweepTarget::casefront){
            tg.reach_ = SweepTarget::plane;
            tg.surf_ = support.frontPlane().get();
            tg.pl_ = support.frontPlane().get();
          } else {
            tg.reach_ = SweepTarget::zplane;
            tg.zval_ = zback;
            tg.surf_ = part == SweepTarget::casehole ? support.inner().get() : support.cfInner().get();
            tg.active_ = false;
          }
          tgts.push_back(tg);
        }
      }
    }
    // calorimeter disks, downstream only: faces and edges, extending to the back face while the track is in the bore
    if(!updir && kkg_h->calo()){
      auto const& calo = *kkg_h->calo();
      for(int idisk = 0; idisk < 2; ++idisk){
        if(!(idisk == 0 ? toCaloD0_ : toCaloD1_)) continue;
        bool d0 = idisk == 0;
        auto const& front = d0 ? calo.EMC_Disk_0_FrontPtr() : calo.EMC_Disk_1_FrontPtr();
        auto const& back = d0 ? calo.EMC_Disk_0_BackPtr() : calo.EMC_Disk_1_BackPtr();
        auto const& inner = d0 ? calo.EMC_Disk_0_InnerPtr() : calo.EMC_Disk_1_InnerPtr();
        auto const& outer = d0 ? calo.EMC_Disk_0_OuterPtr() : calo.EMC_Disk_1_OuterPtr();
        auto face = [&](char const* name, KinKal::Surface const* surf){
          SweepTarget tg;
          tg.sid_ = SurfaceId(name);
          tg.surf_ = surf;
          tg.drives_ = false;
          tg.rmax_ = bore;
          return tg;
        };
        auto ftg = face(d0 ? "EMC_Disk_0_Front" : "EMC_Disk_1_Front",front.get());
        tgts.push_back(ftg);
        auto btg = face(d0 ? "EMC_Disk_0_Back" : "EMC_Disk_1_Back",back.get());
        btg.drives_ = true;
        btg.reach_ = SweepTarget::zplane;
        btg.zval_ = back->center().Z();
        tgts.push_back(btg);
        for(auto const& [name,cyl] : {std::make_pair(d0 ? "EMC_Disk_0_Inner" : "EMC_Disk_1_Inner",inner),
            std::make_pair(d0 ? "EMC_Disk_0_Outer" : "EMC_Disk_1_Outer",outer)}){
          auto ctg = face(name,cyl.get());
          tgts.push_back(ctg);
        }
      }
    }
    return tgts;
  }

  template <class KTRAJ> void KKExtrap::sweep(KKTrack<KTRAJ>& ktrk, SweepTargetCollection const& tgts, TimeDir tdir, double t0, bool updir) const {
    if(tgts.empty() || !ktrk.fitStatus().usable()) return;
    auto const& ftraj = ktrk.fitTraj();
    ExtrapolateSweep xs(tgts,tdir,maxdt_,maxdtstep_,btol_,intertol_,minv_,debug_);
    // material can only be added on the end piece: the sweep starts there, never before t0
    double tend = tdir == TimeDir::forwards ? ftraj.back().range().begin() : ftraj.front().range().end();
    double tstart = tdir == TimeDir::forwards ? std::max(t0,tend) : std::min(t0,tend);
    // calorimeter material starts at the last measurement (the fit can extend past it into the front panel) and is not
    // added again where the fit has it (KKFit AddCaloMaterial); the inner support is only entered forwards in time
    double tmeas = !updir && toCaloMaterial_ ? std::clamp(tdir == TimeDir::forwards ? ktrk.activeRange().end() : ktrk.activeRange().begin(),
        ftraj.range().begin(),ftraj.range().end()) : tstart;
    for(size_t itg = 0; itg < tgts.size(); ++itg){
      auto const& tg = tgts[itg];
      if(tg.effect_ != SweepTarget::calopanel && tg.effect_ != SweepTarget::calosupport) continue;
      xs.notBefore(itg,tmeas);
      bool skip = tg.effect_ == SweepTarget::calosupport && tdir == TimeDir::backwards;
      for(auto const& matxing : ktrk.caloMaterialPlaneXings())
        skip |= tg.effect_ == SweepTarget::calopanel ? matxing->surfaceId() == tg.sid_ : matxing->surfaceId().id() == tg.sid_.id();
      if(tg.effect_ == SweepTarget::calosupport)
        for(auto const& matxing : ktrk.caloMaterialCylXings()) skip |= matxing->surfaceId().id() == tg.sid_.id();
      if(skip) xs.setActive(itg,false);
    }
    SupportEntries entries;
    // record-only surfaces the fitted trajectory already crosses between t0 and the sweep start
    xs.start(t0);
    xs.limit(tstart,true);
    while(!xs.needsExtrapolation(ftraj,tdir) && xs.hasCrossing()){
      addCrossing(ktrk,xs.targets()[xs.crossing().itarget_],xs.crossing().inter_,tdir);
      xs.accept();
    }
    if(xs.failed()){
      addFailure(ktrk,tdir,xs.error());
      return;
    }
    // then extend the end piece from the sweep start, one crossing at a time
    xs.unlimit();
    xs.start(tstart);
    if(updir && upstream_){
      auto const& ktraj = tdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
      if(ktraj.bnom().R() > 0.0){
        double vb = ktraj.velocity(tstart).Dot(ktraj.bnom().Unit())*KinKal::timeDirSign(tdir);
        xs.watchReflection(vb > 0.0 ? 1.0 : -1.0);
      }
    }
    unsigned ncross(0);
    while(ncross < maxCrossings_){
      bool extended = !xs.needsExtrapolation(ftraj,tdir) || ktrk.extrapolate(tdir,xs);
      if(xs.failed()){
        addFailure(ktrk,tdir,xs.error()); // the search failed, before or during KinKal's walk; the fit stays usable
        break;
      }
      if(!extended){
        addFailure(ktrk,tdir,ktrk.fitStatus().comment_); // KinKal marked the fit unusable while extending it
        break;
      }
      if(xs.hasCrossing()){
        size_t itg = xs.crossing().itarget_;
        auto const& tg = xs.targets()[itg];
        bool added = tg.effect_ == SweepTarget::calosupport ? supportCrossing(ktrk,xs,itg,xs.crossing().inter_,tdir,entries)
          : addCrossing(ktrk,tg,xs.crossing().inter_,tdir);
        if(!added) break;
        xs.accept();
        ++ncross;
        if(tg.terminal_ || !ktrk.fitStatus().usable()) break;
      } else if(xs.reflected()) {
        if(!reflect(xs,ftraj,tdir)) break;
      } else {
        break; // nothing left to reach, or MaxDt since the last crossing
      }
    }
    if(ncross >= maxCrossings_) addFailure(ktrk,tdir,"sweep exceeded " + std::to_string(maxCrossings_) + " crossings");
  }

  // A reflected track returns to the tracker or finishes crossing its reflecting element.
  template <class KTRAJ> bool KKExtrap::reflect(ExtrapolateSweep& xs, KinKal::ParticleTrajectory<KTRAJ> const& ftraj, TimeDir tdir) const {
    xs.watchReflection(0.0);
    auto const& tgts = xs.targets();
    if(backToTracker_){
      for(size_t itg = 0; itg < tgts.size(); ++itg){
        if(tgts[itg].sid_ == SurfaceId("TSDA")) xs.setActive(itg,false);
        if(tgts[itg].terminal_ && tgts[itg].sid_ == SurfaceId("TT_Front")) xs.setActive(itg,true);
      }
      return true;
    }
    double z = ftraj.position3(tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin()).Z();
    SweepTarget::Effect keep(SweepTarget::record);
    for(auto const& tg : tgts){
      if(tg.effect_ == SweepTarget::ipashell && z > tg.zreglo_ && z < tg.zreghi_) keep = SweepTarget::ipashell;
      if(tg.effect_ == SweepTarget::stfoil && z < tg.zreghi_) keep = SweepTarget::stfoil;
    }
    if(keep == SweepTarget::record) return false;
    for(size_t itg = 0; itg < tgts.size(); ++itg)
      if(tgts[itg].drives_ && tgts[itg].effect_ != keep) xs.setActive(itg,false);
    return true;
  }

  template <class KTRAJ> bool KKExtrap::addCrossing(KKTrack<KTRAJ>& ktrk, SweepTarget const& tg, KinKal::Intersection const& inter, TimeDir tdir) const {
    if(debug_ > 0) std::cout << "KKExtrap: " << (tdir == TimeDir::forwards ? "forwards " : "backwards ") << tg.sid_ << " " << inter << std::endl;
    if(tg.effect_ == SweepTarget::record){
      ktrk.addIntersection(tg.sid_,inter);
      return true;
    }
    double thick = tg.pipes_ ? tg.pipes_->wallDepth(inter.pos_) : tg.thick_;
    if(tg.crvsector_) thick *= tg.crvsector_->layerFraction(inter.pos_,inter.pdir_); // near a seam, the staggered layers are shared with the next sector
    if(tg.effect_ == SweepTarget::calopanel && thick <= 0.0) return true; // the pipe plane between the pipes
    if(tg.effect_ == SweepTarget::crvsector && thick <= 0.0){ ktrk.addIntersection(tg.sid_,inter); return true; } // inside the envelope but between layers
    // KinKal records why a fit became unusable (e.g. a failed domain walk); the direction just ends
    if(!ktrk.fitStatus().usable()) return false;
    auto const& ftraj = ktrk.fitTraj();
    auto const& endpiece = tdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
    auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
    // material can only be added on the end piece, and the shell must reference that same piece
    if(!endpiece.range().inRange(inter.time_) || ftraj.nearestTraj(inter.time_) != reftrajptr){
      std::ostringstream reason;
      reason << "material crossing " << tg.sid_ << " at " << inter.time_ << " is not owned by the end piece " << endpiece.range();
      addFailure(ktrk,tdir,reason.str());
      return false;
    }
    switch(tg.effect_) {
      case SweepTarget::cylshell:
      case SweepTarget::ipashell: {
        auto shell = std::make_shared<KKShellXing<KTRAJ,KinKal::Cylinder>>(tg.cyl_,tg.sid_,*tg.mat_,inter,reftrajptr,tg.thick_,intertol_);
        if(rangesOut(*shell,tdir)) return rangeOut(ktrk,tg,inter);
        if(tg.effect_ == SweepTarget::cylshell) ktrk.addMaterialCylXing(shell,tdir);
        else ktrk.addIPAXing(shell,tdir);
        return true;
      }
      case SweepTarget::stfoil: {
        auto shell = std::make_shared<KKShellXing<KTRAJ,KinKal::Annulus>>(tg.ann_,tg.sid_,*tg.mat_,inter,reftrajptr,tg.thick_,intertol_);
        if(rangesOut(*shell,tdir)) return rangeOut(ktrk,tg,inter);
        ktrk.addSTXing(shell,tdir);
        return true;
      }
      case SweepTarget::planeshell:
      case SweepTarget::calopanel: {
        auto shell = std::make_shared<KKShellXing<KTRAJ,KinKal::Plane>>(tg.plane_,tg.sid_,*tg.mat_,inter,reftrajptr,thick,intertol_);
        if(rangesOut(*shell,tdir)) return rangeOut(ktrk,tg,inter);
        ktrk.addMaterialPlaneXing(shell,tdir);
        return true;
      }
      case SweepTarget::crvsector: {
        auto shell = std::make_shared<KKShellXing<KTRAJ,KinKal::Rectangle>>(tg.rect_,tg.sid_,*tg.mat_,inter,reftrajptr,thick,intertol_);
        if(rangesOut(*shell,tdir)) return rangeOut(ktrk,tg,inter);
        ktrk.addCRVXing(shell,tdir);
        return true;
      }
      case SweepTarget::record:
        return true;
      case SweepTarget::calosupport:
        break; // supportCrossing places the support material
    }
    throw cet::exception("RECO") << "mu2e::KKExtrap: unexpected sweep effect for " << tg.sid_;
  }

  // A track moving downstream onto a crystal case front enters the case there, or, from the disk hole, through the case
  // inner cylinder (see caloSupportXings). The first material is added at the entry and the other where the track then
  // reaches the carbon fiber ring inner radius, each as a thin shell with the straight-line path through it from the entry.
  template <class KTRAJ> bool KKExtrap::supportCrossing(KKTrack<KTRAJ>& ktrk, ExtrapolateSweep& xs, size_t itg,
      KinKal::Intersection const& inter, TimeDir tdir, SupportEntries& entries) const {
    using Support = KKGeom::CaloInnerSupport;
    auto const& tgts = xs.targets();
    auto const& tg = tgts[itg];
    auto const& support = *tg.support_;
    auto part = [&](SweepTarget::SupportPart sp){
      for(size_t jtg = 0; jtg < tgts.size(); ++jtg)
        if(tgts[jtg].support_ == tg.support_ && tgts[jtg].part_ == sp) return jtg;
      throw cet::exception("RECO") << "mu2e::KKExtrap: missing calorimeter support target for " << tg.sid_;
    };
    if(debug_ > 0) std::cout << "KKExtrap: calo support " << tg.sid_ << " part " << tg.part_ << " " << inter << std::endl;
    auto shell = [&](Support::Crossing xing, double path){
      GeomHandle<mu2e::KKMaterial> kkmat_h;
      SweepTarget mt;
      mt.sid_ = support.surfaceId(xing);
      mt.mat_ = &*kkmat_h->material(support.materialName(Support::material(xing)));
      mt.thick_ = path*std::fabs(inter.norm_.Dot(inter.pdir_));
      if(xing == Support::alFront || xing == Support::cfFront){
        mt.effect_ = SweepTarget::planeshell;
        mt.plane_ = support.front();
      } else {
        mt.effect_ = SweepTarget::cylshell;
        mt.cyl_ = xing == Support::alHole ? support.inner() : support.cfInner();
      }
      return addCrossing(ktrk,mt,inter,tdir);
    };
    auto enter = [&](Support::Crossing first, Support::Crossing second){
      auto paths = support.pathLengths(inter.pos_,inter.pdir_);
      if(paths[Support::material(first)] > 0.0 && !shell(first,paths[Support::material(first)])) return false;
      if(paths[Support::material(second)] > 0.0){
        entries[tg.support_] = SupportEntry{paths[Support::material(second)],second,inter.time_};
        xs.setActive(part(SweepTarget::cfring),true);
      }
      return true;
    };
    switch(tg.part_){
      case SweepTarget::casefront: {
        double rho = support.transverseRadius(inter.pos_);
        if(inter.pdir_.Z() <= 0.0 || rho >= support.rcrystal()) return true; // moving upstream, or onto the crystals, where it showers
        if(rho < support.rin()){
          xs.setActive(part(SweepTarget::casehole),true); // in the disk hole: the case is entered through its inner cylinder
          return true;
        }
        bool oncf = rho >= support.rcf();
        return enter(oncf ? Support::cfFront : Support::alFront, oncf ? Support::alRing : Support::cfRing);
      }
      case SweepTarget::casehole:
        return enter(Support::alHole,Support::cfRing);
      case SweepTarget::cfring: {
        auto const& entry = entries.at(tg.support_);
        if(inter.time_ == entry.time_) return true;
        return shell(entry.xing_,entry.path_);
      }
    }
    return true;
  }
}
#endif
