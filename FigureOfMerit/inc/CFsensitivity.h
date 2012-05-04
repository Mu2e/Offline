#ifndef CFSENSITIVITY_H
#define CFSENSITIVITY_H

// ----------------------------------------------------------------------
//
// CFsensitivity.h
//
// Functions to calculate the 90% CL (or other p-value) sensitivities in the
// sense of the Cousins/Feldman paper 
//
// ----------------------------------------------------------------------

#include <vector>

namespace mu2e {

// CLrangeBlock represents a horizontal block in (x=n, y=mu) spaces such that:
// The likelihood ratio for all n within the block exceeds that for 
//   any n outside the block.
// For mu in the vertical range, the sum, for n within the block, 
//   of Poisson P(n|mu+b) is at least p.
// For mu in the range, the sum, excluding either end value of n, is less 
//  than p.
struct CLrangeBlock 
{
  int n1;
  int n2;
  double mu_low;
  double mu_high;
  CLrangeBlock (int i, int j, double a, double b) 
    : n1(i), n2(j), mu_low(a), mu_high(b) {}
};

class CFsensitivity
{
public:
  CFsensitivity(double pValue);
  double operator() (double bValue);

public:
  double p; 
  double b;
  int maximum_n;
  std::vector<CLrangeBlock> region;  // Region of Confidence Level p 
private:
  CLrangeBlock zeroBlock () const;
  double likelihoodRatio(int n, double mu) const;
  double criticalMu ( int na, int nb, bool & crossoverHappens ) const;
  CLrangeBlock expandBlockHeight ( CLrangeBlock & block ) const;
  double meanMuMax () const;
  double findUpperMuLimit ( int n, bool & valid ) const;
  enum VerbosityLevels {
      CFV_NONE
    , CFV_ULS
    , CFV_BLOCKS
    , CFV_MU
    , CFV_ALL_BUT_R_CHECKS
    , CFV_ALL
  };
  VerbosityLevels verbosity;
}; // CFsensitivity

int main () {
  CFsensitivity c(.90);
  return 0;
}

} // end namespace mu2e

#endif // CLSENSITIVITY_H
