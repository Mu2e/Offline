#include "Offline/Mu2eKinKal/inc/KKStrawMaterial.hh"
#include "Offline/Mu2eKinKal/inc/StrawXingUpdater.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/Detector/MaterialXing.hh"
#include <cmath>
#include <algorithm>
#include "cetlib_except/exception.h"

namespace mu2e {

  KKStrawMaterial::KKStrawMaterial(StrawProperties const& sprops,
      const std::shared_ptr<MatEnv::DetMaterial> wallmat,
      const std::shared_ptr<MatEnv::DetMaterial>gasmat,
      const std::shared_ptr<MatEnv::DetMaterial> wiremat) :
    orad2_(pow(sprops.strawOuterRadius(),2)),
    irad_(sprops.strawInnerRadius()),
    irad2_(pow(sprops.strawInnerRadius(),2)),
    wallonlypath_(sqrt(2.0*sprops.strawOuterRadius()*sprops.strawWallThickness())),
    avggaspath_(0.5*M_PI*sprops.strawInnerRadius()),
    avgwallpath_(0.5*M_PI*(orad2_-irad2_)/sprops.strawOuterRadius()),
    wrad_(sprops.wireRadius()),
    wallmat_(wallmat), gasmat_(gasmat), wiremat_(wiremat)
  {}

  KKStrawMaterial::KKStrawMaterial(MatEnv::MatDBInfo const& matdbinfo,StrawProperties const& sprops,
      const std::string& wallmat, const std::string& gasmat, const std::string& wiremat) :
    KKStrawMaterial(sprops,
        matdbinfo.findDetMaterial(wallmat),
        matdbinfo.findDetMaterial(gasmat),
        matdbinfo.findDetMaterial(wiremat)) {}

  KKStrawMaterial::PathCalc KKStrawMaterial::averagePathLengths(double& wallpath, double& gaspath, double& wirepath) const {
    gaspath = avggaspath_;
    wallpath = avgwallpath_;
    wallpath = 0.0;
    return KKStrawMaterial::average;
  }

  KKStrawMaterial::PathCalc KKStrawMaterial::pathLengths(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig,
      double& wallpath, double& gaspath, double& wirepath) const {
    wallpath = gaspath = wirepath = 0.0;
    PathCalc retval = KKStrawMaterial::unknown;
    double adoca = fabs(cadata.doca());
    double maxvar = pow(caconfig.maxdocasig_,2);
    if (cadata.docaVar() < maxvar ) {  // if measured doca is accurate, use it to compute the exact paths
      if(adoca < caconfig.maxddoca_ && adoca < irad_ ){
        double adoca2 = pow(adoca,2);
        gaspath = 2.0*sqrt(irad2_ - adoca2);
        wallpath = 2.0*sqrt(orad2_ - adoca2) - gaspath;
        retval = exactdoca;
      } else {
        // near the straw boundary, take the average over the outer part of the straw for the wall thickness. Note gas path is 0 here
        wallpath = wallonlypath_;
        retval = truncdoca;
      }
    } else {
      // errors are large WRT the size of the straw, or DOCA is very far from the wire: just take the average over all impact parameters
      gaspath = avggaspath_;
      wallpath = avgwallpath_;
      retval = average;
    }
    // Model the wire as a diffuse gas, density constrained by DOCA TODO
    // correct for the angle WRT the axis
    double afac = angleFactor(cadata.dirDot());
    wallpath *= afac;
    gaspath *= afac;
    return retval;
  }

  double KKStrawMaterial::transitLength(ClosestApproachData const& cadata) const {
    double doca = std::min(fabs(cadata.doca()),irad_);
    double tlen = 2.0*sqrt(orad2_-doca*doca);
    // correct for the angle WRT the axis
    double afac = angleFactor(cadata.dirDot());
    tlen *= afac;
    return tlen;
  }

  double KKStrawMaterial::angleFactor(double dirdot) const {
    // protect against nearly parallel angles
    static const double maxfac2(100.0); // beyond this the straw irregularities dominate and the estimate is unreliable
    static const double minsin2(1.0/(maxfac2)); // minimum sin^2(theta) this implies
    double sin2 = std::max(1.0 -dirdot*dirdot,minsin2);
    return 1.0/sqrt(sin2);
  }

  KKStrawMaterial::PathCalc KKStrawMaterial::findXings(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig,
      std::vector<MaterialXing>& mxings) const {
    mxings.clear();
    double wallpath, gaspath, wirepath;
    auto retval = pathLengths(cadata,caconfig,wallpath, gaspath, wirepath);
    if(wallpath > 0.0) mxings.push_back(MaterialXing(*wallmat_,wallpath));
    if(gaspath > 0.0) mxings.push_back(MaterialXing(*gasmat_,gaspath));
    if(wirepath > 0.0) mxings.push_back(MaterialXing(*wiremat_,wirepath));
    return retval;
  }

}
