#ifndef murat_LsqSums4
#define murat_LsqSums4
///////////////////////////////////////////////////////////////////////////////
// PM: Kasa fit of the circle parameters
// Kasa, I. (1976). A curve fitting procedure and its error analysis.
// IEEE Trans. Inst. Meas. 25 8–14.
// the original Kasa paper is difficult to find in open access,
// so see https://arxiv.org/pdf/0907.0421 by N.Chernov instead
///////////////////////////////////////////////////////////////////////////////
class LsqSums4 {

protected:
  double _qn, sw, sx, sy, sx2, sxy, sy2, sx3, sx2y, sxy2, sy3, sx4, sx3y, sx2y2, sxy3, sy4;

  double fX0, fY0; // offsets, need to be defined in the very beginning, by default - 0


public:
  LsqSums4(double X0 = 0, double Y0 = 0);

  void   clear();

  void   addPoint   (double X, double Y, double W = 1.);
  void   removePoint(double X, double Y, double W = 1.);

  double qn      () { return _qn; }
  double xMean   () { return sx/sw; }
  double yMean   () { return sy/sw; }
  double x2Mean  () { return sx2/sw; }
  double xyMean  () { return sxy/sw; }
  double y2Mean  () { return sy2/sw; }
  double x3Mean  () { return sx3/sw; }
  double x2yMean () { return sx2y/sw; }
  double xy2Mean () { return sxy2/sw; }
  double y3Mean  () { return sy3/sw; }
  double x4Mean  () { return sx4/sw; }
  double x3yMean () { return sx3y/sw; }
  double x2y2Mean() { return sx2y2/sw; }
  double xy3Mean () { return sxy3/sw; }
  double y4Mean  () { return sy4/sw; }

  double sigXX   () { return x2Mean() - xMean()*xMean(); }
  double sigXY   () { return xyMean() - xMean()*yMean(); }
  double sigYY   () { return y2Mean() - yMean()*yMean(); }

  double sigX2X  () { return x3Mean () - xMean()*x2Mean(); }
  double sigX2Y  () { return x2yMean() - yMean()*x2Mean(); }
  double sigXY2  () { return xy2Mean() - xMean()*y2Mean(); }
  double sigYY2  () { return y3Mean () - yMean()*y2Mean(); }

  double sigXX3  () { return x4Mean  ()-xMean ()*x3Mean();}
  double sigX2X2 () { return x4Mean  ()-x2Mean()*x2Mean();}
  double sigX3Y  () { return x3yMean ()-x3Mean()*yMean ();}
  double sigX2Y2 () { return x2y2Mean()-x2Mean()*y2Mean();}
  double sigXY3  () { return xy3Mean ()-xMean ()*y3Mean();}
  double sigY2Y2 () { return y4Mean  ()-y2Mean()*y2Mean();}
  double sigYY3  () { return y4Mean  ()-yMean ()*y3Mean();}

  double det     () { return sigXX()*sigYY() -sigXY()*sigXY(); }

                                        // reconstructed parameters of the circle
  double x0    ();
  double y0    ();
  double radius();
  double phi0  ();
  double dfdz  ();                      // straight line slope

  double chi2DofCircle();
  double chi2DofLine  ();

};

#endif
