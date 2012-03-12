// PoissonCDFsolver.cc

#include "FigureOfMerit/inc/PoissonCDFsolver.h"

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <cassert>

namespace mu2e {

// determine, for a given mu, the first n such that CDF > a given p-value
// This is optimized for pValues close to 1.
int poissonInverseCDF (double mu, double pValue) {
  assert (mu > 0 && pValue > 0 && pValue < 1);
  // guess an m such that the m-th term is about 1-p
  if (pValue < std::exp(-mu)) return 0;
  double lneps = std::log(1-pValue);
  double lnmu  = std::log(mu); 
  //std::cout << "lneps = " << lneps << "  lnmu = " << lnmu << "\n";
  double m0 = lneps/(lnmu-1);
  double m1 = lneps/(1-std::log(m0)+lnmu);
  double m2 = lneps/(1-std::log(m1)+lnmu);
  int m = static_cast<int>(m2);
  //std::cout << "m0 = " << m0 
  //          << "  m1 = " << m1
  //          << "  m2 = " << m2
  //          << "  m = " << m << "\n";
  double mthCDFterm = 1.0;
  for (int k=1 ; k < m+1; ++k) {
    mthCDFterm *= mu/k;
  }
  //std::cout << m << "-th CDF term = " << mthCDFterm << "\n";
  // to get good accuracy, compute the term we will need to go to,
  // then start there and work backwards in accumulating.  
  double tol = 1.0e-15;
  double firstSeriesTerm = mthCDFterm * mu/(m+1);
  double lastSeriesTerm = firstSeriesTerm;
  double endTermMax = tol*firstSeriesTerm; 
  int endk;
  for (endk = m+1; lastSeriesTerm > endTermMax; ++endk) {
    lastSeriesTerm *= mu/(endk+1);
  }
  // Working backwards from that last series term, see when we exceed 1-p
  double target = (1-pValue)*std::exp(mu);
  double seriesTerm = lastSeriesTerm;
  double series = seriesTerm;
  //std::cout << "endk = " << endk << "  lastSeriesTerm = " << lastSeriesTerm
  //          << "  target = " << target << "\n";
  int r;
  for (r = endk; r>0; --r) {
    seriesTerm *= r/mu;
    series += seriesTerm;
  //  std::cout << "r = " << r << "  seriesTerm = " << seriesTerm 
  //            << " series = " << series << "\n"; 
    if (series > target) break;
  }
  return r-1;
}
  
// determine, for a given p-value and n, the mu such that CDF at that n
// is equal to that p-value
double poissonMeanSolver (int n, double pValue) {
#ifdef TRACE_PMS
  std::cout << "poissonMeanSolver (" << n << ",  " << pValue << ")\n"; 
#endif
  assert( n>=0 && pValue > 0 && pValue < 1 );
  if (n==0) {
#ifdef DISPLAY_PMS_RESULTS
    std::cout << "The mu that has a pValue of " << pValue 
              << " at n = " << n << " is " << -std::log(pValue) << "\n";
#endif
    return -std::log(pValue);
  }
  // Use Newton's method, with initial guess mu=1
  double t = 1-pValue;
#ifdef TRACE_PMS
  std::cout <<  " t = " << t << "\n";
#endif
  double oldMuGuess = 0.0;
  double muTolAbs = 1.0e-6;
  double newMuGuess = 1.0;
  if (n > 2) newMuGuess = n;
  for (int i = 0; i < 100; ++i)  {
    if (std::fabs(newMuGuess-oldMuGuess) < muTolAbs) break;
#ifdef TRACE_MUGESSES
    std::cout << "  newMuGuess = " << newMuGuess 
              << "  delta = " << newMuGuess - oldMuGuess;
#endif
    oldMuGuess = newMuGuess;
    double mu = oldMuGuess;
    double ccdf = poissonComplementaryCDF(mu,n);
    double slope = poissonComplementaryCDF(mu,n-1)- ccdf;
#ifdef TRACE_CCDF
    std::cout << "  ccdf = " << ccdf << "  slope = " << slope << "\n";
#endif
    newMuGuess =  std::min(4.0*oldMuGuess, oldMuGuess + (t-ccdf)/slope);
    if (newMuGuess <= 0) newMuGuess = oldMuGuess / 4;
  }
#ifdef DISPLAY_PMS_RESULTS
  std::cout << "The mu that has a pValue of " << pValue 
            << " at n = " << n << " is " << newMuGuess << "\n";
#endif
  return newMuGuess;
}
  

// determine pValue corresponding to given mu and n (for completeness)
double poissonCDF (double mu, int n) 
{
  assert (mu > 0 && n >= 0);
  return 1.0 - poissonComplementaryCDF (mu, n) ;
}
  
// determine 1 - Poisson CDF
double poissonComplementaryCDF (double mu, int n) 
{
  assert (mu > 0 && n >= 0);
  double tol = 1.0e-15;
  double lastCDFterm = 1.0;
  for (int k=1 ; k < n+1; ++k) {
    lastCDFterm *= mu/k;
  }
  if (lastCDFterm == 0) return 0.0;
  // to get good accuracy, compute the term we will need to go to,
  // then start there and work backwards in accumulating.  
  double firstSeriesTerm = lastCDFterm * mu/(n+1);
  double lastSeriesTerm = firstSeriesTerm;
  double endTermMax = tol*firstSeriesTerm; 
  int endk;
  for (endk = n+1; lastSeriesTerm > endTermMax; ++endk) {
    lastSeriesTerm *= mu/(endk+1);
  }
  double seriesTerm = lastSeriesTerm;
  double series = seriesTerm;
  for (int k = endk; k>n+1; --k) {
    seriesTerm *= k/mu;
    series += seriesTerm;
  }
  return series*std::exp(-mu);
}

// determine Feldman/Cousins 90% CL sensitivity given a background mean b
double feldmanCousins90pctSensitivity ( double b ) 
{
  static  FeldmanCousins90pctSensitivity s;
  return s(b);
} 

std::vector<double> FeldmanCousins90pctSensitivity::sensValues_4() {
  // From table XII in HUTP-97/A096 
  // arXiv.physics/971021v2
  // column 3, lines 1-9
  std::vector<double> v;
  v.push_back(2.44);
  v.push_back(2.86);
  v.push_back(3.28);
  v.push_back(3.62);
  v.push_back(3.94);
  v.push_back(4.20);
  v.push_back(4.42);
  v.push_back(4.63);
  v.push_back(4.83);
  return v;
}

std::vector<double> FeldmanCousins90pctSensitivity::sensValues_15() {
  // From table XII in HUTP-97/A096 
  // arXiv.physics/971021v2
  // column 3 (integer values in column 1)
  std::vector<double> v;
  v.push_back(2.44);
  v.push_back(3.28);
  v.push_back(3.94);
  v.push_back(4.42);
  v.push_back(4.83);
  v.push_back(5.18);
  v.push_back(5.53);
  v.push_back(5.90);
  v.push_back(6.18);
  v.push_back(6.49);
  v.push_back(6.76);
  v.push_back(7.02);
  v.push_back(7.28);
  v.push_back(7.51);
  v.push_back(7.75);
  v.push_back(7.99);
  return v;
}

FeldmanCousins90pctSensitivity::FeldmanCousins90pctSensitivity() 
  : grid_4 (0.0,  4.5,  9)
  , grid_15(0.0, 16.0, 16)
  , sensitivity_4  ( sensValues_4(),  grid_4  )
  , sensitivity_15 ( sensValues_15(), grid_15 )
{ }
  
double FeldmanCousins90pctSensitivity::operator() (double b) const 
{ 
  if (b < 0)    return  sensitivity_4(0);
  if (b < 4.0)  return  sensitivity_4(b);
  if (b < 15.0) return  sensitivity_15(b);
  // formula for b > 15 may or may not be particularly accurate
  return 2.132 + 1.404 * std::sqrt(b) + 0.046*b;
} 
  


} // end namespace mu2e

