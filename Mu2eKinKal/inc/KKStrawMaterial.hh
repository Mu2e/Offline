#ifndef Mu2eKinKal_KKStrawMaterial_hh
#define Mu2eKinKal_KKStrawMaterial_hh
//
//  description of a local segment of a straw, including a
//  mixture for the straw, the gas, and the wire,
//
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/TrackerGeom/inc/StrawProperties.hh"
#include "Offline/Mu2eKinKal/inc/StrawXingUpdater.hh"

namespace MatEnv {
  class MatDBInfo;
  class DetMaterial;
}
namespace KinKal {
  class MaterialXing;
}
namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::MaterialXing;
  using MatEnv::DetMaterial;
  using MatEnv::MatDBInfo;

  class KKStrawMaterial {
    public:
      enum PathCalc {exactdoca=0,truncdoca,average,unknown};
      // construct from geometry and explicit materials
      KKStrawMaterial(StrawProperties const& sprops,
          const std::shared_ptr<DetMaterial> wallmat_,
          const std::shared_ptr<DetMaterial> gasmat_,
          const std::shared_ptr<DetMaterial> wiremat_);
      // construct using materials by name
      KKStrawMaterial(MatDBInfo const& matdbinfo,StrawProperties const& sprops,
          const std::string& wallmat="straw-wall", const std::string& gasmat="straw-gas", const std::string& wiremat="straw-wire");
      // pathlength through straw components, given closest approach
      void pathLengths(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig, double& wallpath, double& gaspath, double& wirepath) const;
      // transit length given closest approach
      double transitLength(ClosestApproachData const& cadata) const;
      // find the material crossings given doca and error on doca.  Should allow for straw and wire to have different axes TODO
      void findXings(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig, std::vector<MaterialXing>& mxings) const;
      double outerRadius() const { return sprops_.strawOuterRadius(); }
      double innerRadius() const { return sprops_.strawInnerRadius(); }
      double wallThickness() const { return sprops_.strawWallThickness(); }
      double wireRadius() const { return sprops_.wireRadius(); }
      DetMaterial const& wallMaterial() const { return *wallmat_; }
      DetMaterial const& gasMaterial() const { return *gasmat_; }
      DetMaterial const& wireMaterial() const { return *wiremat_; }
      auto pathCalculation() const& { return pathcalc_; }
    private:
      StrawProperties const& sprops_;
      double orad2_; // outer radius of the straw squared
      double irad2_; // inner radius of the straw squared
      double tpath_; // outer radius * thickness = average path at the tangent point
      mutable PathCalc pathcalc_ = unknown; // how were path lenghts calculated?
      const std::shared_ptr<DetMaterial> wallmat_; // material of the straw wall
      const std::shared_ptr<DetMaterial> gasmat_; // material of the straw gas
      const std::shared_ptr<DetMaterial> wiremat_; // material of the wire
      // utility to calculate material factor given the cosine of the angle of the particle WRT the straw
      double angleFactor(double dirdot) const;
  };

}
#endif
