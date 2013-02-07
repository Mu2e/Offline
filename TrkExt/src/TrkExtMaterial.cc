//
//  $Id: TrkExtMaterial.cc,v 1.2 2013/02/07 02:09:47 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2013/02/07 02:09:47 $
//
//  Original author MyeongJae Lee
//
// Note : here we define various parameters which are used in material effect calculation.
// The algorithm should be updated to a more rubust algorithm: TODO
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "TrkExt/inc/TrkExtMaterial.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"

using namespace CLHEP;

using namespace std;

namespace mu2e {




 
  TrkExtMaterial::TrkExtMaterial( ) :
    _name("undefined")
  {
    _matid = Undefined;
  }
    
    
  TrkExtMaterial::TrkExtMaterial( string n) :
    _name(n)
  {
    if (n.find("PE")>=0) {
      _dpmp[0] = -0.719604;
      _dpmp[1] = 0.904384;
      _dpmp[2] = -0.0136965;
      _dpmp[3] = 1.18018;
      _dpmp[4] = 1.66496e-5;
      _dpmp[5] = 2.56602e-6;
      _dpsp[0] = 2.01486;
      _edpsp[0] = 0.0098517;
      _dpsp[1] = 0.0220202;
      _edpsp[1] = 0.000184769;
      _dpsp[2] = -1.19048e-06;
      _edpsp[2] = 8.58465e-07;
      _thpar[0] = 0.6063;
      _thpar[1] = 0.7636;
      _thpar[2] = 0.0875;
      _rho =  0.89; //g cm-3
      _matid = PE;
    }
    else if (n.find("Al")>=0) {
      _dpmp[0] = -0.435264;
      _dpmp[1] = 1.02322;
      _dpmp[2] = -0.00790427;
      _dpmp[3] = 1.05155;
      _dpmp[4] = -1.45739e-5;
      _dpmp[5] = -1.64281e-6;
      _dpsp[0] = 1.54196;
      _edpsp[0] = 0.00855373;
      _dpsp[1] = 0.040025;
      _edpsp[1] = 0.000160425;
      _dpsp[2] = 1.66667e-06;
      _edpsp[2] = 7.45356e-07;
      _thpar[0] = 1.4420;
      _thpar[1] = 0.8294;
      _thpar[2] = 0.0875;
      _rho = 2.699;
      _matid = Al;
    }
    else if (n.find("Vac")>=0) {
      _dpmp[0] = 0;
      _dpmp[1] = 0;
      _dpmp[2] = 0;
      _dpmp[3] = 0;
      _dpmp[4] = 0;
      _dpmp[5] = 0;
      _dpsp[0] = 0;
      _dpsp[1] = 0;
      _dpsp[2] = 0;
      _thpar[0] = 0;
      _thpar[1] = 0;
      _thpar[2] = 0;
      _rho = 0;
      _matid = Vac;
    }
    else {
      _dpmp[0] = 0;
      _dpmp[1] = 0;
      _dpmp[2] = 0;
      _dpmp[3] = 0;
      _dpmp[4] = 0;
      _dpmp[5] = 0;
      _dpsp[0] = 0;
      _dpsp[1] = 0;
      _dpsp[2] = 0;
      _thpar[0] = 0;
      _thpar[1] = 0;
      _thpar[2] = 0;
      _rho = 0;
      _matid = Undefined;
    }


  }


  double TrkExtMaterial::meanEnergyLoss (const Hep3Vector& p, double ds) {
    double ke = safeSqrt(p.mag2() + _mec22) - _mec2;
    return _rho * (_dpsp[0] + _dpsp[1]*ke + _dpsp[2]*ke*ke) * 0.1 * ds;
    return 0;
  }

  double TrkExtMaterial::mostProbableEnergyLoss (const Hep3Vector & p, double ds) {
    if (_matid == Vac || _matid == Undefined) return 0;
    double mom = p.mag() - 104.;
    double thick = log10(fabs(ds));
    double de = pow(10., (_dpmp[0]+_dpmp[1]*thick+_dpmp[2]*thick*thick)*(_dpmp[3]+_dpmp[4]*mom+_dpmp[5]*mom*mom));
    return de;
  }

  double TrkExtMaterial::scatteringAngle (const CLHEP::Hep3Vector& p, double ds) {
    double fabsds = fabs(ds);
    if (fabsds <=0) return 0;
    return _thpar[0] * safeSqrt(fabsds) / p.mag() * (_thpar[1] + _thpar[2] * log10(fabsds)) ; 
  }



} // end namespace mu2e

