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
      enum PathCalc {range,average,unknown};
      // construct from geometry and explicit materials
      KKStrawMaterial(StrawProperties const& sprops,
          const std::shared_ptr<DetMaterial> wallmat_,
          const std::shared_ptr<DetMaterial> gasmat_,
          const std::shared_ptr<DetMaterial> wiremat_);
      // construct using materials by name
      KKStrawMaterial(MatDBInfo const& matdbinfo,StrawProperties const& sprops,
          const std::string& wallmat="straw-wall", const std::string& gasmat="straw-gas", const std::string& wiremat="straw-wire");
      // pathlength through straw components, given closest approach. Return the method used to compute the paths
      PathCalc pathLengths(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig, double& wallpath, double& gaspath, double& wirepath) const;
      PathCalc averagePathLengths(double& wallpath, double& gaspath, double& wirepath) const;
      // transit length given closest approach
      double transitLength(ClosestApproachData const& cadata) const;
      // find the material crossings given doca and error on doca.  Should allow for straw and wire to have different axes TODO
      PathCalc findXings(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig, std::vector<MaterialXing>& mxings) const;
      DetMaterial const& wallMaterial() const { return *wallmat_; }
      DetMaterial const& gasMaterial() const { return *gasmat_; }
      DetMaterial const& wireMaterial() const { return *wiremat_; }
      double wireRadius() const { return wrad_; }
    private:
      double irad_, orad_; // inner and outer radii
      double avggaspath_, avgwallpath_; // average wall and gas paths
      double wrad_; // wire radius
      const std::shared_ptr<DetMaterial> wallmat_; // material of the straw wall
      const std::shared_ptr<DetMaterial> gasmat_; // material of the straw gas
      const std::shared_ptr<DetMaterial> wiremat_; // material of the wire
      // utility to calculate material factor given the cosine of the angle of the particle WRT the straw
      static double angleFactor(double dirdot);
      static double segmentArea(double doca, double r); // area of a circular segment defined by the cord DOCA to the center
  };

}
#endif
