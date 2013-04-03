Double_t lCrystalball (Double_t *x, Double_t *par) {
	  // par[0] : norm
	  // par[1] : x0
	  // par[2] : sigma
	  // par[3] : n
	  // par[4] : alpha

	  const double invSqrt2pi=0.398942280401432702863;
	  double DeltaX = (x[0]- par[1]);

	  double tailFval=par[0]*invSqrt2pi/par[2];
	  double tailAbsSigma = fabs(par[2]);
	  DeltaX*=-1.0;
	  if ( DeltaX/tailAbsSigma > -1.*par[4]) {
	    tailFval *= TMath::Gaus(x[0], par[1], par[2]);
	  }
	  else {
	    double tailAbsAlpha = fabs(par[4]);
	    double tailA = pow(par[3]/tailAbsAlpha, par[3])*exp(-0.5*par[4]*par[4]);
	    double tailB = par[3]/tailAbsAlpha - tailAbsAlpha;
	    tailFval *= tailA*pow(tailB-DeltaX/tailAbsSigma, -1.*par[3]);
	  }

	  return tailFval;
}
