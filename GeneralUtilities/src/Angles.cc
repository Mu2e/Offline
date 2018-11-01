#include "GeneralUtilities/inc/Angles.hh"
namespace mu2e {
  namespace Angles {
    double deltaPhi(double& phi, double refphi) {
      double dphi = phi - refphi;
      static const double twopi = 2*M_PI;
      while(dphi > M_PI){
	dphi -= twopi;
	phi -= twopi;
      }
      while(dphi <= -M_PI){
	dphi += twopi;
	phi += twopi;
      }
      return dphi;
    }
    float deltaPhi(float& phi, float refphi) {
      double dophi(phi);
      double dorefphi(refphi);
      return deltaPhi(dophi,dorefphi);
    }
  }
}
