//-----------------------------------------------------------------------------
// LsqSums4
//-----------------------------------------------------------------------------
#include <cmath>

#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"

LsqSums4::LsqSums4(double X0, double Y0) {
  clear();
  fX0 = X0;
  fY0 = Y0;
}

void LsqSums4::clear() {
  _qn   = 0;
  sw    = 0;
  sx    = 0;
  sy    = 0;
  sx2   = 0;
  sxy   = 0;
  sy2   = 0;
  sx3   = 0;
  sx2y  = 0;
  sxy2  = 0;
  sy3   = 0;
  sx4   = 0;
  sx3y  = 0;
  sx2y2 = 0;
  sxy3  = 0;
  sy4   = 0;

  fX0   = 0;
  fY0   = 0;
}

//-----------------------------------------------------------------------------
void LsqSums4::addPoint(double XX, double YY, double W) {
  // move to COG to improve numerical accuracy
  const double X = XX-fX0;
  const double Y = YY-fY0;

  // Only compute these factors once
  const double X2   = X * X;
  const double Y2   = Y * Y;
  const double XY   = X * Y;
  const double W_X  = W * X;
  const double W_Y  = W * Y;
  const double W_X2 = W * X2;
  const double W_Y2 = W * Y2;

  _qn   += 1;
  sw    += W;          // W
  sx    += W_X;        // W * X
  sy    += W_Y;        // W * Y
  sx2   += W_X2;       // W * X * X
  sxy   += W*XY;       // W * X * Y
  sy2   += W_Y2;       // W * Y * Y
  sx3   += W_X2 * X;   // W * X * X * X
  sx2y  += W_X2 * Y;   // W * X * X * Y
  sxy2  += W_Y2 * X;   // W * Y * Y * X
  sy3   += W_Y2 * Y;   // W * Y * Y * Y
  sx4   += W_X2 * X2;  // W * X * X * X * X
  sx3y  += W_X2 * XY;  // W * X * X * X * Y
  sx2y2 += W_X2 * Y2;  // W * X * X * Y * Y
  sxy3  += W_Y2 * XY;  // W * Y * Y * Y * X
  sy4   += W_Y2 * Y2;  // W * Y * Y * Y * Y

  // _qn   += 1;
  // sw    += W;
  // sx    += X*W;
  // sy    += Y*W;
  // sx2   += X*X*W;
  // sxy   += X*Y*W;
  // sy2   += Y*Y*W;
  // sx3   += X*X*X*W;
  // sx2y  += X*X*Y*W;
  // sxy2  += X*Y*Y*W;
  // sy3   += Y*Y*Y*W;
  // sx4   += X*X*X*X*W;
  // sx3y  += X*X*X*Y*W;
  // sx2y2 += X*X*Y*Y*W;
  // sxy3  += X*Y*Y*Y*W;
  // sy4   += Y*Y*Y*Y*W;
}

void LsqSums4::removePoint(double XX, double YY, double W) {
  // move to COG to improve numerical accuracy
  const double X = XX-fX0;
  const double Y = YY-fY0;

  // Only compute these factors once
  const double X2   = X * X;
  const double Y2   = Y * Y;
  const double XY   = X * Y;
  const double W_X  = W * X;
  const double W_Y  = W * Y;
  const double W_X2 = W * X2;
  const double W_Y2 = W * Y2;

  _qn   -= 1;
  sw    -= W;          // W
  sx    -= W_X;        // W * X
  sy    -= W_Y;        // W * Y
  sx2   -= W_X2;       // W * X * X
  sxy   -= W*XY;       // W * X * Y
  sy2   -= W_Y2;       // W * Y * Y
  sx3   -= W_X2 * X;   // W * X * X * X
  sx2y  -= W_X2 * Y;   // W * X * X * Y
  sxy2  -= W_Y2 * X;   // W * Y * Y * X
  sy3   -= W_Y2 * Y;   // W * Y * Y * Y
  sx4   -= W_X2 * X2;  // W * X * X * X * X
  sx3y  -= W_X2 * XY;  // W * X * X * X * Y
  sx2y2 -= W_X2 * Y2;  // W * X * X * Y * Y
  sxy3  -= W_Y2 * XY;  // W * Y * Y * Y * X
  sy4   -= W_Y2 * Y2;  // W * Y * Y * Y * Y
}

//-----------------------------------------------------------------------------
double LsqSums4::x0() {
  double x;
  x  = (sigYY()*(sigX2X()+sigXY2())-sigXY()*(sigX2Y()+sigYY2()))/2/det();
  return x + fX0;
}

double LsqSums4::y0() {
  double y;
  y = (sigXX()*(sigX2Y()+sigYY2())-sigXY()*(sigX2X()+sigXY2()))/2/det();
  return y + fY0;
}

double LsqSums4::radius () {
  double r, x_0, y_0, dx, dy;

  x_0 = x0()-fX0;
  y_0 = y0()-fY0;

  dx = xMean()-x_0;
  dy = yMean()-y_0;

  r  = sqrt(sigXX()+sigYY()+dx*dx+dy*dy);

  return r;
}

double LsqSums4::phi0(){
  double phi0, D;
  D = sw*sx2 - sx*sx;

  phi0 = sx2*sy - sx*sxy;
  phi0 /= D;

  return phi0;
}

//-----------------------------------------------------------------------------
// the straight line slope
//-----------------------------------------------------------------------------
double LsqSums4::dfdz(){
  double dfdz, D;
  D     = sw*sx2 - sx*sx;
  dfdz  = sw*sxy - sy*sx;
  dfdz /= D;

  return dfdz;
}

double LsqSums4::chi2DofCircle() {

  double chi2, /*x_0, y_0,*/ r, sx2, sy2;

//   x_0 = x0();
//   y_0 = y0();
  r   = radius();

  sx2 = sigX2X()+sigXY2();
  sy2 = sigX2Y()+sigYY2();

  chi2 = sigX2X2()+2.*sigX2Y2()+sigY2Y2()-(sigYY()*sx2*sx2+sigXX()*sy2*sy2-2*sigXY()*sx2*sy2)/det();
  //-----------------------------------------------------------------------------
  // normalization, assuming spread of points << R
  //-----------------------------------------------------------------------------
  chi2 = chi2/(4*r*r);
  chi2 *= sw/_qn;

  return chi2;
}

double LsqSums4::chi2DofLine() {

  double chi2;
  chi2  = sigYY()*sigXX() - sigXY()*sigXY();
  chi2 /= sigXX();
  chi2 *= sw/_qn;

 //  double chi2_new = sigYY() - dfdz()*sigXY();
//   printf("[LsqSum4::chi2rphiDofCircle] chi2 = %5.3e chi2_new = %5.3e\n", chi2 , chi2_new);

  return chi2;
}
