//-----------------------------------------------------------------------------
// LsqSums2
//-----------------------------------------------------------------------------
#include "Mu2eUtilities/inc/LsqSums2.hh"


LsqSums2::LsqSums2() {
  clear();
}

LsqSums2::LsqSums2(const LsqSums2& S) {
  init(S);
  //  printf(" LsqSums2::LsqSums2 : dont use me!\n");
}

LsqSums2::~LsqSums2() {
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

void LsqSums2::init(const LsqSums2& S) {
  _qn   = S._qn  ;
  sw    = S.sw   ;
  sx    = S.sx   ;
  sy    = S.sy   ;
  sx2   = S.sx2  ; 
  sxy   = S.sxy  ;
  sy2   = S.sy2  ;

  fX0   = S.fX0;
  fY0   = S.fY0;
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

double LsqSums2::dydx() {

  double dfdz(0), D;
  if (_qn > 1) {
    D = sw*sx2 - sx*sx;				
  
    dfdz = sw*sxy - sy*sx;
    dfdz /= D;
  }
  
  return dfdz;
}

double LsqSums2::chi2Dof() {

  double chi2;
  chi2  = sigYY()*sigXX() - sigXY()*sigXY();
  chi2 /= sigXX();
  chi2 *= sw/_qn;

 //  double chi2_new = sigYY() - dfdz()*sigXY();
//   printf("[LsqSum4::chi2rphiDofCircle] chi2 = %5.3e chi2_new = %5.3e\n", chi2 , chi2_new);
  
  return chi2;
}
