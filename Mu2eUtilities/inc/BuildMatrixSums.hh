#ifndef _MU2E_UTILITIES_BUILDMATRIXSUMS_HH
#define _MU2E_UTILITIES_BUILDMATRIXSUMS_HH

#include "Math/VectorUtil.h"
#include "TMatrixD.h"
class BuildMatrixSums {

protected:
  double betaX00, betaX10, gammaX00, gammaX01, gammaX11;
  double betaY00, betaY10, gammaY00, gammaY01, gammaY11;

public:
  BuildMatrixSums();
  BuildMatrixSums(const BuildMatrixSums& Sums);
  ~BuildMatrixSums();

  void   clear();
  void   init(const BuildMatrixSums& S);
  void addPoint(int i, XYZVec point_i, XYZVec track_direction, double errX, double errY);
  void removePoint(XYZVec point_i, XYZVec track_direction, double errX, double errY);

  double Get2DParameter(int i, TMatrixD Alpha);

  TMatrixD GetAlphaX();
  TMatrixD GetGammaX();
  TMatrixD GetBetaX();
  TMatrixD GetAlphaY();
  TMatrixD GetGammaY();
  TMatrixD GetBetaY();

};

#endif
