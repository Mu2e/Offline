#ifndef _MU2E_UTILITIES_BUILDLINEARFITMATRIXSUMS_HH
#define _MU2E_UTILITIES_BUILDLINEARFITMATRIXSUMS_HH

#include "Math/VectorUtil.h"
#include "TMatrixD.h"
class BuildLinearFitMatrixSums {

protected:
  double betaX00, betaX10, gammaX00, gammaX01, gammaX11;
  double betaY00, betaY10, gammaY00, gammaY01, gammaY11;
  double deltaX, deltaY;
  double chi2;
public:
  BuildLinearFitMatrixSums();
  BuildLinearFitMatrixSums(const BuildLinearFitMatrixSums& Sums);
  ~BuildLinearFitMatrixSums();

  void   clear();
  void   init(const BuildLinearFitMatrixSums& S);
  void   addPoint(XYZVec point_i, XYZVec XPrime, XYZVec YPrime,XYZVec ZPrime, double errX, double errY);
 

  double Get2DParameter(int i, TMatrixD Alpha);
  void SetChi2(double chi2);
  TMatrixD GetAlphaX();
  TMatrixD GetGammaX();
  TMatrixD GetCovX();
  TMatrixD GetBetaX();
  TMatrixD GetAlphaY();
  TMatrixD GetGammaY();
  TMatrixD GetBetaY();
  TMatrixD GetCovY();
  
  double GetChi2X();
  double GetChi2Y();
  double GetTotalChi2();


};

#endif
