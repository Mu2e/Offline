#ifndef murat_LsqSums2
#define murat_LsqSums2
///////////////////////////////////////////////////////////////////////////////
// circle fit
// Author: P. Murat
// Date:
///////////////////////////////////////////////////////////////////////////////

class LsqSums2 {

protected:
  double _qn, sw, sx, sy, sx2, sxy, sy2;
  double fX0, fY0; // offsets, need to be defined in the very beginning, by default - 0


public:
  LsqSums2(double X0=0,double Y0=0);

  void   clear();

  void   addPoint   (double X, double Y, double W = 1.);
  void   removePoint(double X, double Y, double W = 1.);
  void   addSum     (LsqSums2& Lsq);
  void   removeSum  (LsqSums2& Lsq);

  double qn      () { return _qn; }
  double xMean   () { return sx/sw; }
  double yMean   () { return sy/sw; }
  double x2Mean  () { return sx2/sw; }
  double xyMean  () { return sxy/sw; }
  double y2Mean  () { return sy2/sw; }

  double sigXX   () { return x2Mean() - xMean()*xMean(); }
  double sigXY   () { return xyMean() - xMean()*yMean(); }
  double sigYY   () { return y2Mean() - yMean()*yMean(); }

  double det     () { return sigXX()*sigYY() -sigXY()*sigXY(); }

  // reconstructed parameters of the line
  // note: error computations assume that weights on data points (x,y) with error on y were set using weight = 1/error^2
  double dydx();
  double dydxErr();
  double y0();
  double y0Err();
  double chi2Dof();
  //  ClassDef(LsqSums2,0)

};

#endif
