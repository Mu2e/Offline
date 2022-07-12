#include "Offline/Mu2eKinKal/inc/KKStrawMaterial.hh"
#include "Offline/Mu2eKinKal/inc/StrawXingUpdater.hh"
#include <cmath>
#include <algorithm>
#include "cetlib_except/exception.h"

namespace mu2e {

  KKStrawMaterial::KKStrawMaterial(StrawProperties const& sprops,
      const MatEnv::DetMaterial *wallmat,
      const MatEnv::DetMaterial *gasmat,
      const MatEnv::DetMaterial *wiremat) :
    sprops_(sprops),
    wallmat_(wallmat), gasmat_(gasmat), wiremat_(wiremat) {
      // compute some caches
      srad2_ = strawRadius()*strawRadius();
      grad2_ = gasRadius()*gasRadius();
    }

  KKStrawMaterial::KKStrawMaterial(MatEnv::MatDBInfo const& matdbinfo,StrawProperties const& sprops,
      const char* wallmat, const char* gasmat, const char* wiremat) :
    KKStrawMaterial(sprops,
        matdbinfo.findDetMaterial(wallmat),
        matdbinfo.findDetMaterial(gasmat),
        matdbinfo.findDetMaterial(wiremat)) {}

  void KKStrawMaterial::pathLengths(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig,
      double& wallpath, double& gaspath, double& wirepath) const {
    wallpath = gaspath = wirepath = 0.0;
    double adoca = fabs(cadata.doca());
    double docasig = sqrt(cadata.docaVar());
    if(adoca < caconfig.maxdoca_ || caconfig.hitstate_.active()){
      if (docasig < caconfig.maxdocasig_ && adoca < caconfig.maxddoca_) {  // use exact calculation based on DOCA
        double doca = std::min(adoca, gasRadius()); // truncate
        double ddoca = doca*doca;
        gaspath = 2.0*sqrt(grad2_- ddoca);
        wallpath = 2.0*wallThickness()*strawRadius()/sqrt(srad2_-ddoca);
      } else {
        // errors are large WRT the size of the straw, or DOCA is very far from the wire: just take the average over all impact parameters
        gaspath = M_PI_2*strawRadius();
        wallpath = M_PI*wallThickness();
      }
      if(isnan(wallpath) || isnan(gaspath))throw cet::exception("RECO")<<"mu2e::KKStrawMaterial: Invalid pathlength" << std::endl;
      // Model the wire as a diffuse gas, density constrained by DOCA TODO
      // correct for the angle WRT the axis
      double afac = angleFactor(cadata.dirDot());
      wallpath *= afac;
      gaspath *= afac;
    }
  }

  double KKStrawMaterial::transitLength(ClosestApproachData const& cadata) const {
    double doca = std::min(fabs(cadata.doca()),gasRadius());
    double tlen = 2.0*sqrt(srad2_-doca*doca);
    // correct for the angle WRT the axis
    double afac = angleFactor(cadata.dirDot());
    tlen *= afac;
    return tlen;
  }

  double KKStrawMaterial::angleFactor(double dirdot) const {
    // protect against nearly parallel angles
    static const double maxfac(10.0); // beyond this the straw irregularities dominate and the estimate is unreliable
    static const double minsin2(1.0/(maxfac*maxfac)); // minimum sin^2(theta) this implies
    double sin2 = std::max(1.0 -dirdot*dirdot,minsin2);
    return 1.0/sqrt(sin2);
  }

  void KKStrawMaterial::findXings(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig, std::vector<MaterialXing>& mxings) const {
    mxings.clear();
    double wallpath, gaspath, wirepath;
    pathLengths(cadata,caconfig,wallpath, gaspath, wirepath);
    if(wallpath > 0.0) mxings.push_back(MaterialXing(*wallmat_,wallpath));
    if(gaspath > 0.0) mxings.push_back(MaterialXing(*gasmat_,gaspath));
    if(wirepath > 0.0) mxings.push_back(MaterialXing(*wiremat_,wirepath));
  }

}
