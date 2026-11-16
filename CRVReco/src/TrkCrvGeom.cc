#include "Offline/CRVReco/inc/TrkCrvGeom.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <map>

namespace mu2e {

  std::vector<TrkCrvSector> trkCrvSectors() {
    GeomHandle<CosmicRayShield> crs;
    GeomHandle<KinKalGeom> kkg;
    auto const& shields = crs->getCRSScintillatorShields();
    std::vector<TrkCrvSector> sectors(shields.size());
    std::map<std::string,int> nameToShield;
    for(size_t i = 0; i < shields.size(); ++i) {
      sectors[i].name = shields[i].getName();
      nameToShield[sectors[i].name] = static_cast<int>(i);
    }
    for(auto const& sec : kkg->CRV()->sectors()) {
      auto it = nameToShield.find(sec.sname_);
      if(it == nameToShield.end()) continue;
      auto& g = sectors[it->second];
      g.valid = true;
      g.rect = sec.sector_;
      g.center = sec.sector_->center();
      g.u = sec.sector_->uDirection();
      g.v = sec.sector_->vDirection();
      g.w = sec.sector_->normal();
      g.uHalf = sec.sector_->uHalfLength();
      g.vHalf = sec.sector_->vHalfLength();
      try { g.sid = SurfaceId(sec.sname_); }
      catch(std::exception const&) { mf::LogWarning("TrkCrvGeom") << "no SurfaceId for CRV sector " << sec.sname_; }
    }
    return sectors;
  }

  int trkCrvDominantShield(CrvCoincidenceCluster const& cluster, CosmicRayShield const& crs, bool* multiSector) {
    std::map<int,double> pes;
    for(auto const& rp : cluster.GetCrvRecoPulses()) pes[crs.getBar(rp->GetScintillatorBarIndex()).id().getShieldNumber()] += rp->GetPEs();
    int shield = -1;
    double best = -1;
    for(auto const& [s, pe] : pes) if(pe > best) { best = pe; shield = s; }
    if(multiSector) *multiSector = pes.size() > 1;
    return shield;
  }
}
