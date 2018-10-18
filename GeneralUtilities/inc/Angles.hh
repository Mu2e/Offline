//
//  Angular manipulation utilities.  These mostly deal with angle variables obey different arithemtic
//  (x=0; y=2pi; x==y is true) than assumed by standard math libraries.
//  Radians are assumed unless otherwise explicitly stated
//  Original author: Dave Brown (LBNL) 8/9/2016
//
#ifndef GeneralUtilities_Angles_hh
#define GeneralUtilities_Angles_hh
#include <math.h>
namespace mu2e {
  namespace Angles {
   // find difference between 2 angles in the range (-pi,pi]
   // Note that this also updates the input angle phase so that dphi = phi - refphi
    double deltaPhi(double& phi, double refphi=0.0);
    float deltaPhi(float& phi, float refphi=0.0);
  }
}
#endif
