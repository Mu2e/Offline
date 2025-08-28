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
    irad_(sprops.strawInnerRadius()),
    orad_(sprops.strawOuterRadius()),
    avggaspath_(0.5*M_PI*irad_),
    avgwallpath_(M_PI*sprops.strawWallThickness()),
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
    double docarange = 0.5*caconfig.nsig_*sqrt(std::max(0.0,cadata.docaVar()));
    // if the doca range covers the straw, use averages
    double adoca = fabs(cadata.doca());
    double mindoca = std::min(std::max(0.0,adoca-docarange),irad_);
    double imaxdoca = std::min(irad_,adoca+docarange);
    double omaxdoca = std::min(orad_,adoca+docarange);
    if(2*docarange > irad_){
      gaspath = avggaspath_;
      wallpath = avgwallpath_;
      retval = average;
    } else {
      // path average is the area covered by the DOCA ranged divided by the range
      double agas = segmentArea(mindoca,irad_) - segmentArea(imaxdoca,irad_);
      if(agas>0.0)gaspath = agas/(imaxdoca-mindoca);
      // wall area is the outer radius area - gas area
      double awall = segmentArea(mindoca,orad_) - segmentArea(omaxdoca,orad_) - agas;
      if(awall>0.0)wallpath = awall/(omaxdoca-mindoca);
      retval = range;
    }
    if(caconfig.diag_>0){
      std::cout << "KKStrawMaterial: DOCA " << fabs(cadata.doca()) << " DOCA Var " << cadata.docaVar() << " range [" << mindoca << "," << imaxdoca << "] wall path (mm)" << wallpath << " gas path " << gaspath << " retval " << retval << std::endl;
    }

    // Model the wire as a diffuse gas, density constrained by DOCA TODO
    // 3D pathlength includes projection along straw
    double afac = angleFactor(cadata.dirDot());
    wallpath *= afac;
    gaspath *= afac;
    return retval;
  }

  double KKStrawMaterial::transitLength(ClosestApproachData const& cadata) const {
    double doca = std::min(fabs(cadata.doca()),irad_);
    double tlen = 2.0*sqrt(orad_*orad_-doca*doca);
    // correct for the angle WRT the axis
    double afac = angleFactor(cadata.dirDot());
    tlen *= afac;
    return tlen;
  }

  KKStrawMaterial::PathCalc KKStrawMaterial::findXings(ClosestApproachData const& cadata,StrawXingUpdater const& caconfig,
      std::vector<MaterialXing>& mxings) const {
    mxings.clear();
    double wallpath, gaspath, wirepath;
    auto retval = pathLengths(cadata,caconfig,wallpath, gaspath, wirepath);
    if(wallpath > 0.0) mxings.push_back(MaterialXing(*wallmat_,wallpath));
    if(gaspath > 0.0) mxings.push_back(MaterialXing(*gasmat_,gaspath));
    if(wirepath > 0.0) mxings.push_back(MaterialXing(*wiremat_,wirepath));
    if(caconfig.diag_ > 1){
      std::cout << "Ce wall dE/dx " << wallmat_->energyLoss(105,1.0,0.511) << std::endl;
      std::cout << "Ce gas dE/dx " << gasmat_->energyLoss(105,1.0,0.511) << std::endl;
    }
    return retval;
  }

  double KKStrawMaterial::angleFactor(double dirdot){
    // protect against nearly parallel angles
    static const double maxfac2(100.0); // beyond this the straw irregularities dominate and the estimate is unreliable
    static const double minsin2(1.0/(maxfac2)); // minimum sin^2(theta) this implies
    double sin2 = std::max(1.0 -dirdot*dirdot,minsin2);
    return 1.0/sqrt(sin2);
  }

  double KKStrawMaterial::segmentArea(double doca, double r){
    if(doca<r){
      double theta = 2*acos(doca/r); // angle of segment
      return 0.5*r*r*(theta-sin(theta));
    } else
      return 0.0;
    }

}
