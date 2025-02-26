#ifndef Mu2eKinKal_KKStrawMaterial_hh
#define Mu2eKinKal_KKStrawMaterial_hh
//
//  description of a local segment of a straw, including a
//  mixture for the straw, the gas, and the wire,
//
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/Detector/MaterialXing.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/TrackerGeom/inc/StrawProperties.hh"
#include "Offline/Mu2eKinKal/inc/StrawXingUpdater.hh"

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::MaterialXing;
  class KKStrawMaterial {
    public:
      // construct from geometry and explicit materials
      KKStrawMaterial(StrawProperties const& sprops,
          const std::shared_ptr<MatEnv::DetMaterial> wallmat_,
          const std::shared_ptr<MatEnv::DetMaterial> gasmat_,
          const std::shared_ptr<MatEnv::DetMaterial> wiremat_);
      // construct using materials by name
      KKStrawMaterial(MatEnv::MatDBInfo const& matdbinfo,StrawProperties const& sprops,
        const std::string& wallmat="straw-wall", const std::string& gasmat="straw-gas", const std::string& wiremat="straw-wire");
      // pathlength through straw components, given closest approach
      void pathLengths(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig, double& wallpath, double& gaspath, double& wirepath) const;
      // transit length given closest approach
      double transitLength(ClosestApproachData const& cadata) const;
      // find the material crossings given doca and error on doca.  Should allow for straw and wire to have different axes TODO
      void findXings(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig, std::vector<MaterialXing>& mxings) const;
      double gasRadius() const { return sprops_.strawInnerRadius(); }
      double strawRadius() const { return sprops_.strawOuterRadius(); }
      double wallThickness() const { return sprops_.strawWallThickness(); }
      double wireRadius() const { return sprops_.wireRadius(); }
      MatEnv::DetMaterial const& wallMaterial() const { return *wallmat_; }
      MatEnv::DetMaterial const& gasMaterial() const { return *gasmat_; }
      MatEnv::DetMaterial const& wireMaterial() const { return *wiremat_; }
    private:
      StrawProperties const& sprops_;
      double srad2_; // average outer transverse radius of the straw squared
      double grad2_; // effective gas volume radius squared
      const std::shared_ptr<MatEnv::DetMaterial> wallmat_; // material of the straw wall
      const std::shared_ptr<MatEnv::DetMaterial> gasmat_; // material of the straw gas
      const std::shared_ptr<MatEnv::DetMaterial> wiremat_; // material of the wire
      // utility to calculate material factor given the cosine of the angle of the particle WRT the straw
      double angleFactor(double dirdot) const;
      // maximum DOCA given straw irregularities
  };

}
#endif
