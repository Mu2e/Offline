#include "Offline/KinKalGeom/inc/CaloInnerSupport.hh"
namespace mu2e {
  namespace KKGeom {
    using KinKal::VEC3;

    CaloInnerSupport::CaloInnerSupport(SurfaceIdEnum::enum_type sid, VEC3 const& center, double rin, double rcf, double rcrystal,
        double zfront, double zback, std::array<std::string,nmaterials> const& materials) :
      sid_(sid), center_(center.X(),center.Y(),zfront), rcrystal_(rcrystal),
      frontplane_(std::make_shared<KinKal::Plane>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),center_)),
      front_(std::make_shared<KinKal::Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),center_,rin,rcrystal)),
      inner_(std::make_shared<KinKal::Cylinder>(VEC3(0.0,0.0,1.0),VEC3(center.X(),center.Y(),0.5*(zfront+zback)),rin,0.5*(zback-zfront))),
      cfinner_(std::make_shared<KinKal::Cylinder>(VEC3(0.0,0.0,1.0),VEC3(center.X(),center.Y(),0.5*(zfront+zback)),rcf,0.5*(zback-zfront))),
      materials_(materials) {}

    std::array<double,CaloInnerSupport::nmaterials> CaloInnerSupport::pathLengths(VEC3 const& pos, VEC3 const& dir) const {
      // step along the line: the rings are a few mm to cm thick, and the path is at most a few cm
      static const double step(0.1); // mm
      static const double maxpath(300.0); // mm
      double const rin = inner_->radius();
      double const zmin = center_.Z();
      double const zmax = zmin + 2.0*inner_->halfLength();
      std::array<double,nmaterials> paths{};
      auto udir = dir.Unit();
      // start half a step inside, so a point on the entry surface isn't counted as being in the hole
      for(double s = 0.5*step; s < maxpath; s += step){
        auto pt = pos + s*udir;
        double r = transverseRadius(pt);
        if(r < rin || r >= rcrystal_ || pt.Z() < zmin || pt.Z() > zmax) break;
        for(auto const& ring : rings_){
          if(r >= ring.rin_ && r < ring.rout_ && pt.Z() >= ring.zmin_ && pt.Z() < ring.zmax_){
            paths[ring.mat_] += step;
            break;
          }
        }
      }
      return paths;
    }
  }
}
