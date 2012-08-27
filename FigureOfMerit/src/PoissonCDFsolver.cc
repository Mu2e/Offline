// PoissonCDFsolver.cc

#include "FigureOfMerit/inc/PoissonCDFsolver.h"

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <cassert>

//#define DISPLAY_PMS_RESULTS
//#define TRACE_MUGESSES
//#define TRACE_CCDF
//#define TRACE_PMS
//#define TRACE_INV_CDF

namespace mu2e {

// determine, for a given mu, the first n such that CDF > a given p-value
// This is optimized for pValues close to 1.
int poissonInverseCDF (double mu, double pValue) {
  assert (mu > 0 && pValue > 0 && pValue < 1);

  // Math comment - it would be more accurate to determine sums by starting
  // from the smallest term -- which will be high k -- and working downward.
  // However, double precision is adequate for p5sigma (with about 9 digits to 
  // spare) so to avoid algoritm instability for non-small mu, we will keep
  // our calculations simple instead.

  double term = std::exp(-mu);
  if (term > pValue) return 0;
  double sum  = term;
  int i = 1;
  for (; i < 200; ++i) {
    term *= mu/i;
    sum += term;
    if ( sum > pValue ) break;
  }
  return i;

#ifdef REMOVED_OLD_INCORRECT_ALGORITHM
#ifdef TRACE_INV_CDF
  std::cout << "poissonInverseCDF: mu = " << mu << " p = " << pValue 
            << " 1-p = " << 1-pValue << "\n";
#endif
  // guess an m such that the m-th term is about 1-p
  if (pValue < std::exp(-mu)) return 0;
  double lneps = std::log(1-pValue);
  double lnmu  = std::log(mu); 
#ifdef TRACE_INV_CDF
  std::cout << "lneps = " << lneps << "  lnmu = " << lnmu << "\n";
#endif
  double m0 = lneps/(lnmu-1);
  double m1 = lneps/(1-std::log(m0)+lnmu);
  double m2 = lneps/(1-std::log(m1)+lnmu);
  int m = static_cast<int>(m2);
#ifdef TRACE_INV_CDF
   std::cout << "m0 = " << m0 
            << "  m1 = " << m1
            << "  m2 = " << m2
            << "  m = " << m << "\n";
#endif
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
#ifdef TRACE_INV_CDF
     std::cout << "r = " << r << "  seriesTerm = " << seriesTerm 
               << " series = " << series << "\n"; 
#endif
    if (series > target) break;
  }
#ifdef TRACE_INV_CDF
     std::cout << "result = " << r-1 << "\n"; 
#endif  
  return r-1;
#endif // REMOVED_OLD_INCORRECT_ALGORITHM
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
    double slope;
    if (n > 0) {
      slope = poissonComplementaryCDF(mu,n-1)- ccdf;
    } else {
      slope = 1 - ccdf;
    }
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
} // poissonMeanSolver
  
// determine, for a given p-value and (n1,n2), the mu such that the prob
// slice from n1 to n2 is is equal to that p-value.  This takes the hint that
// the solution is in the interval [mu_1,mu_2]
double poissonMeanSliceSolver (int n1, int n2, double pValue, 
                               double mu_1, double mu_2)
{
#ifdef TRACE_PMS
  std::cout << "poissonMeanSliceSolver (" << n1 << ", " << n2 
            << ", " << pValue << ", " << mu_1 << ", " << mu_2 << ")\n"; 
#endif
  assert( n1>=0 && n2>=n1 && pValue > 0 && pValue < 1 );

 // up to here


  if ( (n1==0) && (n2==0) ) {
#ifdef DISPLAY_PMS_RESULTS
    std::cout << "The mu that has a slice probability of " << pValue 
              << " at (" << n1 << ", " << n2  
              << ") is " << -std::log(pValue) << "\n";
#endif
    return -std::log(pValue);
  }
  // Use Newton's method, with initial guess (mu_1+mu_2)/2
  double t = 1-pValue;
#ifdef TRACE_PMS
  std::cout <<  " t = " << t << "\n";
#endif
  double oldMuGuess = mu_1;
  double muTolAbs = 1.0e-6;
  double newMuGuess = (mu_1+mu_2)/2;
  for (int i = 0; i < 100; ++i)  {
    if (std::fabs(newMuGuess-oldMuGuess) < muTolAbs) break;
#ifdef TRACE_MUGESSES
    std::cout << "  newMuGuess = " << newMuGuess 
              << "  delta = " << newMuGuess - oldMuGuess;
#endif
    oldMuGuess = newMuGuess;
    double mu = oldMuGuess;
    double ccdf  = 1 - poissonProbabilitySlice(mu,n1, n2);
    double slope;
    if (n1 > 0) {
      slope  = 1 - ccdf - poissonProbabilitySlice(mu,n1-1,n2-1);
    } else if ( n2 > 0) {
      slope = 1 - ccdf - poissonProbabilitySlice(mu,0,n2-1);
    } else {
      slope = 1 - ccdf;
    }
#ifdef TRACE_CCDF
    std::cout << "  ccdf = " << ccdf << "  slope = " << slope << "\n";
#endif
    newMuGuess =  std::min(.75*mu_2+.25*oldMuGuess, 
                                oldMuGuess + (t-ccdf)/slope);
    if (newMuGuess <= mu_1) newMuGuess = .75*mu_1 + .25*oldMuGuess;
    if ( i == 99 ) {
      std::cout << "???? poissonMeanSliceSolver (" << n1 << ", " << n2 
            << ", " << pValue << ", " << mu_1 << ", " << mu_2 
            << ") did not converge\n"; 
    }
  }
#ifdef DISPLAY_PMS_RESULTS
  std::cout << "The mu that has a slice probability of " << pValue 
              << " at (" << n1 << ", " << n2  
              << ") is " << newMuGuess << "\n";
#endif
  return newMuGuess;
} // poissonMeanSliceSolver

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

// Determine the probability sum from n1 to n2
double poissonProbabilitySlice (double mu, int n1, int n2) {
  if (n1 == 0) {
    return poissonCDF (mu,n2);
  } else {
    return poissonCDF (mu,n2) - poissonCDF (mu,n1-1);
    // Yes, this can be done almost twice as efficiently...
  }
}


} // end namespace mu2e

