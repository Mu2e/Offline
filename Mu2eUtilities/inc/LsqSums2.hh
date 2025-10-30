#ifndef murat_LsqSums2
#define murat_LsqSums2
///////////////////////////////////////////////////////////////////////////////
// accumulation of the LSQ sums for the straight line fit
///////////////////////////////////////////////////////////////////////////////
class LsqSums2 {
protected:
  double _qn, sw, sx, sy, sx2, sxy, sy2;
  double fX0, fY0; // offsets, need to be defined in the very beginning, by default - 0
public:
  LsqSums2(double X0 = 0, double Y0 = 0);

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
  double dydx();
  double y0(double X=0);
                                        // parameter errors are calculated assuming
                                        // the data point weights defined as weight=1./err^2
  double dydxErr();
                                        // to take into account correlations
  double y0Err  (double X=0);
  double chi2Dof();
};

#endif
