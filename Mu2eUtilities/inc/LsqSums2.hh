#ifndef murat_LsqSums2
#define murat_LsqSums2
///////////////////////////////////////////////////////////////////////////////
// circle fit
// Author: P. Murat 
// Date:      
///////////////////////////////////////////////////////////////////////////////

class LsqSums2 {

protected:
  double _qn, sw, sx, sy, sx2, sxy, sy2, sx3, sx2y, sxy2, sy3, sx4, sx3y, sx2y2, sxy3, sy4;
  double xmean, ymean, x2mean, xymean, y2mean, x3mean, x2ymean, xy2mean, y3mean;
  double x4mean, x3ymean, x2y2mean,  xy3mean, y4mean;

  double fX0, fY0; // offsets, need to be defined in the very beginning, by default - 0
  

public:
  LsqSums2();
  LsqSums2(const LsqSums2& Sums);

  ~LsqSums2();

  void   clear();
  void   init(const LsqSums2& S);

  void   addPoint   (double X, double Y, double W = 1.);
  void   removePoint(double X, double Y, double W = 1.);

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

  double chi2Dof();
  //  ClassDef(LsqSums2,0)

};

#endif
