//-----------------------------------------------------------------------------
// LsqSums2
//-----------------------------------------------------------------------------
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include <cmath>


LsqSums2::LsqSums2(double X0, double Y0) {
  clear();
  fX0 = X0;
  fY0 = Y0;
}


void LsqSums2::clear() {
  _qn   = 0;
  sw    = 0;
  sx    = 0;
  sy    = 0;
  sx2   = 0;
  sxy   = 0;
  sy2   = 0;

  fX0   = 0;
  fY0   = 0;
}

void LsqSums2::addPoint(double XX, double YY, double W) {
  double X, Y;
  // move to COG
  X = XX-fX0;
  Y = YY-fY0;

  _qn   += 1;
  sw    += W;
  sx    += X*W;
  sy    += Y*W;
  sx2   += X*X*W;
  sxy   += X*Y*W;
  sy2   += Y*Y*W;
}

void LsqSums2::removePoint(double XX, double YY, double W) {
  // move to COG
  double X, Y;
  X = XX-fX0;
  Y = YY-fY0;

  _qn   -= 1;
  sw    -= 1*W;
  sx    -= X*W;
  sy    -= Y*W;
  sx2   -= X*X*W;
  sxy   -= X*Y*W;
  sy2   -= Y*Y*W;
}

void LsqSums2::addSum(LsqSums2& Lsq) {

  _qn   += Lsq._qn;
  sw    += Lsq.sw;
  sx    += Lsq.sx;
  sy    += Lsq.sy;
  sx2   += Lsq.sx2;
  sxy   += Lsq.sxy;
  sy2   += Lsq.sy2;
}

void LsqSums2::removeSum(LsqSums2& Lsq) {

  _qn   -= Lsq._qn;
  sw    -= Lsq.sw;
  sx    -= Lsq.sx;
  sy    -= Lsq.sy;
  sx2   -= Lsq.sx2;
  sxy   -= Lsq.sxy;
  sy2   -= Lsq.sy2;
}

double LsqSums2::dydx() {

  double dfdz(0), D;
  if (_qn > 1) {
    D = sw*sx2 - sx*sx;

    dfdz = sw*sxy - sy*sx;
    dfdz /= D;
  }

  return dfdz;
}


double LsqSums2::dydxErr(){
  double slopeErr(0), N, D;
  if (_qn > 2) {
    N = sw;
    D = sw*sx2 - sx*sx;
    slopeErr = std::sqrt(N/D);
  }

  return slopeErr;
}

double LsqSums2::y0(){
  double y0, D;
  D = sw*sx2 - sx*sx;

  y0 = sx2*sy - sx*sxy;
  y0 /= D;

  return y0;
}

double LsqSums2::y0Err(){
  double interceptErr(0), N, D;
  if (_qn > 2) {
  N = sx2;
  D = sw*sx2 - sx*sx;

  interceptErr = std::sqrt(N/D);
  }

  return interceptErr;
}

double LsqSums2::chi2Dof() {

  double chi2(0);
  if (_qn > 2) {
    chi2  = sigYY()*sigXX() - sigXY()*sigXY();
    chi2 /= sigXX();
    chi2 *= sw/(_qn-2);
  }
  //  double chi2_new = sigYY() - dfdz()*sigXY();
  //  printf("[LsqSum4::chi2rphiDofCircle] chi2 = %5.3e chi2_new = %5.3e\n", chi2 , chi2_new);

  return chi2;
}
