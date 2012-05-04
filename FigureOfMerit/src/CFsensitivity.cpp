// CFsensitivity.cpp

#include "FigureOfMerit/inc/CFsensitivity.h"
#include "FigureOfMerit/inc/PoissonCDFsolver.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

namespace mu2e {

CFsensitivity::CFsensitivity(double pValue) 
  : p(pValue)
  , verbosity(CFV_ULS)
{
  if (verbosity > CFV_NONE) std::cout << "\nHello from CFsensitivity \n\n";
}

double CFsensitivity::operator() (double bValue)
{
  region.clear();
  b = bValue;
  if (verbosity == CFV_ALL) {
    std::cout << "operator() with b = " << b << "\n";
  }
  maximum_n = b + 20*std::sqrt(b) + 15;
  
  // Start with a valid block bottomed at mu=0
  // (without establishing the top of the block)
  CLrangeBlock block = zeroBlock (); 
  if (verbosity == CFV_ALL) {
    std::cout << "after zeroBlock(), b = " << b << "\n";
  }
  if (verbosity >= CFV_BLOCKS) {
    std::cout << "Block including mu = 0 goes from " << block.n1
              << " to " << block.n2 << "\n";
  }  

  // Determine the top of the block 
  CLrangeBlock nextBlock = expandBlockHeight (block);
  if (verbosity == CFV_ALL) {
    std::cout << "after expandBlockHeight (block), b = " << b << "\n";
  }
  if (verbosity  >= CFV_BLOCKS) {
    std::cout << "Block including mu = 0 expands from mu = 0 to " 
              << block.mu_high << "\n"
              << "nextBlock is [" << nextBlock.n1 << ", " 
                << nextBlock.n2 << "] " << nextBlock.mu_low << "\n";  
  }  
 
  // Keep stacking the blocks until finished
//  std::string sin("...");
  while  ( nextBlock.n2 <= maximum_n )
  {
    if (verbosity  >= CFV_BLOCKS) {
    }
    region.push_back (block);
    if (verbosity  >= CFV_BLOCKS) {
      std::cout << "Block added to CL region is [" << block.n1 
                << ", " << block.n2 << "] with mu range (" 
                << block.mu_low << ", " << block.mu_high << ") \n";
      std::cout << "Probability slice in this block is " 
                << poissonProbabilitySlice(block.mu_low+b, block.n1, block.n2)
                << " --> " 
                << poissonProbabilitySlice(block.mu_high+b, block.n1, block.n2)
                << "\n"; 
    }  
    block = nextBlock;
    nextBlock = expandBlockHeight (block);
    if (verbosity  >= CFV_BLOCKS) {
      std::cout << "nextBlock returned is [" << nextBlock.n1 << ", " 
                << nextBlock.n2 << "] " << nextBlock.mu_low << "  "
                << "\n";
    }  

//    std::cout << sin;
//    std::cin  >> sin;
  
  }
  region.push_back (block);
  
  // Perform probability sum based on the CL region
  if (verbosity  >= CFV_BLOCKS) {
    std::cout << "meanMuMax based on these blocks is " << meanMuMax() << "\n";
  }  
  return meanMuMax(); 
}

CLrangeBlock CFsensitivity::zeroBlock() const 
{
  if (verbosity == CFV_ALL) {
    std::cout << "zeroBlock() \n";
  }
  int n1;
  int n2;
  if ( b <= -std::log(p) ) {
    n1 = 0;
    n2 = 0;
    if (verbosity >= CFV_ALL_BUT_R_CHECKS) {
      std::cout << "zeroBlock() returns " << n1 << "  " << n2 << "\n";
    }
    return CLrangeBlock(n1, n2, 0.0, 0.0);
  }
  int nb = static_cast<int>(std::floor(b));
  if ( nb == 0 ) {
    n1 = 0;
    n2 = 1;
    // TODO - check zero and near-zero cases.  
    // TODO - likelihoodRatio probably should involve mu and maybe b
  } else {
    double rminus = likelihoodRatio(nb-1,0);
    double rplus  = likelihoodRatio(nb+1,0);
    if (rplus <= rminus) {
      n1 = nb-1;
      n2 = nb;
    } else {
      n1 = nb;
      n2 = nb+1;
    }
  }
  // compute term1 = poisson probability at n1;
  double term1 = std::exp(-b);
  for (int j = 1; j<= n1; ++j) {
    term1 *= b/j;
  }
  assert (term1 < p);
  assert ( n2 == n1+1 );
  double term2 = term1 * b/n2;
  double sum =  term1 + term2;
    // Note -- the method below is less numerically precise than alternatives 
    // that calculate 1-p by summing contributions *outside* the range.  But
    // even for p corresponding to 5 sigma, double precision is adequate 
    // for this to method yield precise results. 
  while ( sum < p ) {
    if (n1 == 0) {
      // no choice but to add on the next higher n value
      ++n2;
      term2 *= b/n2;  // Note that we divide AFTER incrementing n2
      sum += term2;
    } else {
      // might need to add the next higher or the next lower n value    
      double rlow  = likelihoodRatio(n1-1,0);
      double rhigh = likelihoodRatio(n2+1,0);
      if ( rlow > rhigh ) {
        term1 *= n1/b;
        --n1;        // Note that we multiply BEFORE decrementing n1
        sum += term1; 
      } else {
        ++n2;
        term2 *= b/n2;  // ... but we divide AFTER incrementing n2
        sum += term2;
      }
    } // n1 != 0 branch
    assert (sum <= 1);
  } // while (sum < p)
  if (verbosity >= CFV_ALL_BUT_R_CHECKS) {
    std::cout << "zeroBlock() returns " << n1 << "  " << n2 << "\n";
  }
  return CLrangeBlock(n1,n2,0.0,0.0);
}

double CFsensitivity::likelihoodRatio(int n, double mu) const
{
  if ( n >= b ) {
    return std::exp(n-mu-b) * std::pow((mu+b)/n, n);
  } else if (n > 0) {
    return std::exp(-mu) * std::pow((mu+b)/b, n);;
  } else {
    return std::exp(-mu);
  }
}

double CFsensitivity::criticalMu 
        ( int na, int nb, bool & crossoverHappens ) const {
  if (verbosity >= CFV_MU) {
    std::cout << "Critical mu for (" << na << ", " << nb << ") ";
  }
  assert (nb > na);
  double r;
  double mu;
  crossoverHappens = true;
  // Cover the "easy" cases first:
  if ( (na < 0) || (nb < b) ) {
    crossoverHappens = false;
    if (verbosity >= CFV_MU) {
      std::cout << "No crossover! \n";
    }
    return 0;
  }
  if ( (b==0) && (na==0) ) {
    if (verbosity >= CFV_MU) {
      std::cout << " (b=0) => " << nb/std::exp(1.0) << "\n";
    }
    return nb/std::exp(1.0);
  }
  
  // Now the meaty cases:
  if (na == 0) {
    r = std::pow(nb,nb)*std::exp(b-nb); 
    if (verbosity >= CFV_MU) {
      std::cout << "(na = 0) r = " << r;
    }
  } else if (na < b) {
    r = std::pow(nb,nb)/std::pow(b,na)*std::exp(b-nb);
    if (verbosity >= CFV_MU) {
      std::cout << "(na < b) r = " << r;
    }
  } else {
    r = std::pow(nb,nb)/std::pow(na,na)*std::exp(na-nb);
    if (verbosity >= CFV_MU) {
      std::cout << "(na >= b) r = " << r;
    }
  }
  mu = std::pow(r,(1.0/(nb-na))) - b;
  if (verbosity >= CFV_MU) {
    std::cout << " => " << mu << "\n";
  }
  return mu;
} // criticalMu

CLrangeBlock CFsensitivity::expandBlockHeight (CLrangeBlock & block) const {
  // Will find how high mu can grow before this block no longer satisfies
  // the conditions for a block.  This can happen in one of two ways:
  // (a)  total probability in the block becomes less than p
  // (b)  likelihood ratio for n2+1 becomes higher than for n1
  // (b') likelihood ratio for n1-1 becomes higher than for n2
  // (b') probably is impossible but we will check for it.   
  // Will also determine n1, n2, and mu_low for the next block. 
  
  // When does mu begin to violate block condition for reason (b) or (b'):
  
  double n1 = block.n1;
  double n2 = block.n2;
  double mu_low = block.mu_low;
  
  if (verbosity > CFV_ALL_BUT_R_CHECKS) {
    if (likelihoodRatio(n2+1,mu_low) > likelihoodRatio(n1,mu_low)) {
      std::cout << "Unexpected: likelihoodRatio( " <<  n2+1
                << ", " << mu_low << ") = " << likelihoodRatio(n2+1,mu_low)
                << " > likelihoodRatio( " << n1 << ", " << mu_low 
                << ") = " << likelihoodRatio(n1,mu_low) << "\n";
    }
  }
  bool crossover_1;
  double totalProb_1 = 0;
  double totalProb_1x = 0;
  double totalProb_1a = 0;
  double totalProb_1s = 0;
  double mu_1 = criticalMu (n1,n2+1, crossover_1);
  if (crossover_1) {
    if (verbosity == CFV_ALL) {
      std::cout << "crossover_1 happens at mu = " << mu_1 
                << " with mu_low = " << block.mu_low << "\n";
    }
    if (mu_1 < block.mu_low + 0.0000001) {
      if (verbosity == CFV_ALL) {
        std::cout << "crossover_1 declared redundant \n";
      }
      crossover_1 = false;
    }
    if (crossover_1) {
      totalProb_1  = poissonProbabilitySlice (mu_1+b, n1+1, n2+1);
      totalProb_1x = poissonProbabilitySlice (mu_1+b, n1, n2+1);
      totalProb_1a = poissonProbabilitySlice (mu_1+b, n1, n2);
      totalProb_1s = poissonProbabilitySlice (mu_1+b, n1+1, n2);
      if (verbosity == CFV_ALL) {
        std::cout << "totalProb_1, based on " << mu_1 << ", [" << n1+1 
                  << ", " << n2+1 << "] is " << totalProb_1 << "\n";
        std::cout << "totalProb_1x, based on " << mu_1 << ", [" << n1 
                  << ", " << n2+1 << "] is " << totalProb_1x << "\n";
        std::cout << "totalProb_1a, based on " << mu_1 << ", [" << n1 
                  << ", " << n2 << "] is " << totalProb_1a << "\n";
        std::cout << "totalProb_1s, based on " << mu_1 << ", [" << n1+1 
                  << ", " << n2 << "] is " << totalProb_1s << "\n";
      }     
    }
  }
  bool crossover_2;
  double totalProb_2 = 0;
  double totalProb_2x = 0;
  double totalProb_2a = 0;
  double totalProb_2s = 0;
  double mu_2 = criticalMu (n1-1,n2, crossover_2);
  if (crossover_2) {
    if (verbosity == CFV_ALL) {
      std::cout << "crossover_2 happens at mu = " << mu_2 
                << " with mu_low = " << block.mu_low << "\n";
    }
    if (mu_2 < block.mu_low + 0.0000001) {
      if (verbosity == CFV_ALL) {
        std::cout << "crossover_2 declared redundant \n";
      }
      crossover_2 = false;
    }
    if (crossover_2) {
      totalProb_2  = poissonProbabilitySlice (mu_2+b, n1-1, n2-1);
      totalProb_2x = poissonProbabilitySlice (mu_2+b, n1-1, n2);
      totalProb_2a = poissonProbabilitySlice (mu_2+b, n1, n2);
      totalProb_2s = poissonProbabilitySlice (mu_2+b, n1, n2-1);
      if (verbosity == CFV_ALL) {
        std::cout << "totalProb_2, based on " << mu_2 << ", [" << n1-1 
                  << ", " << n2-1 << "] is " << totalProb_2 << "\n";
        std::cout << "totalProb_2x, based on " << mu_2 << ", [" << n1-1 
                  << ", " << n2 << "] is " << totalProb_2x << "\n";
        std::cout << "totalProb_2a, based on " << mu_2 << ", [" << n1 
                  << ", " << n2 << "] is " << totalProb_2a << "\n";
        std::cout << "totalProb_2s, based on " << mu_2 << ", [" << n1 
                  << ", " << n2-1 << "] is " << totalProb_2s << "\n";
      }     
    }
  }
  assert ( crossover_1 || crossover_2 );
  bool caseb      = false;
  bool casebprime = false;
  bool casea      = false;
  double mu_top_a = 0;
  if ( crossover_1 && (!crossover_2) ) {
    if ( totalProb_1a > p) {
      assert (totalProb_1x > p);
      caseb = true;
    } else {
      casea = true;
      mu_top_a = mu_1;
    }
  }
  if ( (!crossover_1) && crossover_2 ) {
    if ( totalProb_2a > p ) {
      assert (totalProb_2x > p);
      casebprime = true;
    } else {
      casea = true;
      mu_top_a = mu_2;
    }
  }
  if ( crossover_1 && crossover_2 ) {
    if ( (totalProb_1a >  p) && (mu_1 <= mu_2) ) {
      caseb = true;
    }
    if ( (totalProb_2a >  p) && (mu_1 >  mu_2) ) {
      casebprime = true;
    }
    if ( (totalProb_1a <= p) && (mu_1 <= mu_2) ) {
      casea = true;
      mu_top_a = mu_1;
    }
    if ( (totalProb_2a <= p) && (mu_1 >  mu_2) ) {
      casea = true;
      mu_top_a = mu_2;
    }
  }
  assert (casea || caseb || casebprime);
  double mu = -1; 
  int nn1 = -1;
  int nn2 = -1;
  if (caseb) {
    assert (!(casea || casebprime));
    if (verbosity == CFV_ALL) {
      std::cout << "case b \n";
    }     
    if (totalProb_1 > p) {
      // mu might be less than mu_1, if prob 1s is high enough -- in that case,
      // we want to find where the probability **became** high enough:
      mu = mu_1;
      nn1 = n1+1;
      nn2 = n2+1;
      if (totalProb_1s > p) {
        mu = poissonMeanSliceSolver 
                        (n1+1, n2, p, block.mu_low+b, mu_1+b) - b;
        nn1 = n1+1;
        nn2 = n2;
      }
    } else {
      if (verbosity == CFV_ALL) {
        std::cout << "need to expand rather than shift: totalProb_1 = "
                  << totalProb_1 << "\n";
      }     
      mu  = mu_1;
      nn1 = n1;
      nn2 = n2+1;
    }
  }
  if (casebprime) {
    assert (!(casea || caseb));
    if (verbosity == CFV_ALL) {
      std::cout << "case bprime \n";
    }     
    if (verbosity >= CFV_BLOCKS) {
      std::cout << "Mildly unexpected expansion to left: n1 = " << n1
                << " n2 = " << n2 << " mu_1 = " << mu_1 
                << " mu_2 = " << mu_2 << "\n"; 
    }
    if (totalProb_2 > p) {
      // mu might be less than mu_2, if prob 2s is high enough -- in that case,
      // we want to find where the probability **became** high enough:
      mu = mu_2;
      nn1 = n1-1;
      nn2 = n2-1;
      if (totalProb_2s > p) {
        mu = poissonMeanSliceSolver 
                        (n1, n2-1, p, block.mu_low+b, mu_2+b) - b;
        nn1 = n1;
        nn2 = n2-1;
      }
    } else {
      if (verbosity == CFV_ALL) {
        std::cout << "need to expand rather than shift: totalProb_2 = "
                  << totalProb_2 << "\n";
      }     
      mu  = mu_2;
      nn1 = n1-1;
      nn2 = n2;
    }
  }
  if (caseb || casebprime) {
    assert (!casea);
    if (verbosity == CFV_ALL) {
      std::cout << "crossover happens at mu = " << mu 
                << " at [" << nn1 << ", " << nn2 << "]\n";
    }
    assert (nn1 >= 0 && nn2 >= 0); 
    assert (mu >  block.mu_low);
    block.mu_high = mu;
    if (verbosity == CFV_ALL) {
      std::cout << "totalProb is high enough, so block becomes ["
                << n1 << ", " << n2  << "] (" << block.mu_low
                << ", " << block.mu_high << ")\n";
      std::cout << "... and nextBlock will be [" 
                << nn1 << ", " << nn2  << "] (" << block.mu_high
                << ", " << block.mu_high << ")\n";      
    }
    return CLrangeBlock(nn1, nn2, block.mu_high, block.mu_high); 
  }
  assert (casea);
  assert (!(caseb || casebprime));

  // case (a):
  // Sometimes, the mu limitation comes because probability drops down to p.
  // In that case, the probability at the computed mu for any crossover in
  // likelihood ratios will be below p.  Still, we will have to add one or the 
  // other of the neighboring integers, to boost the slice probability back
  // above p.
  if (casea) { 
    if (verbosity == CFV_ALL) {
      std::cout << "case a \n";
    }     
    if (verbosity >= CFV_BLOCKS) {
      std::cout << "(a) Total probability for (" << n1 << ", " << n2 
                << ") at mu = "
                << mu_top_a << " is " 
                << poissonProbabilitySlice (mu_top_a+b, n1, n2) 
                << "\n";
    }
    double mu_p = poissonMeanSliceSolver 
                        (n1, n2, p, block.mu_low+b, mu_top_a+b) - b;

    if (mu_p <= block.mu_low) {
      std::cout << "Impossible situation:  mu_p = " << mu_p
                << " mu_low = " << block.mu_low << " ???? \n";
      // Here we have a situation where you can't increase mu to get to the
      // boundary of the block.  This looks impossible. 
      assert (mu_p > block.mu_low);
    }

    if (verbosity >= CFV_BLOCKS) {
      std::cout << "Probability slice is " << p << " at mu = " << mu_p << "\n";
  #ifdef NOTDEF
      if (n1 > 0) {
        std::cout << "Check: " <<   poissonCDF (mu_p + b, n2) << " - " 
                  <<  poissonCDF (mu_p + b, n1-1) << " = " 
                  <<  poissonProbabilitySlice(mu_p + b, n1, n2) << "\n";
      } else {
        std::cout << "Check: " <<   poissonCDF (mu_p + b, n2) << " = " 
                  <<  poissonProbabilitySlice(mu_p + b, n1, n2) << "\n";
      }        
  #endif
    }
    block.mu_high = mu_p;

    if (block.n1 == 0) {
      // In order to increase probability we need to add some n; and the 
      // only one availble is n2+1
      if (verbosity >= CFV_BLOCKS) {
        std::cout << "since n1 is 0, nextBlock increases n2 to " 
                  << n2+1 << "\n";
      }
      return CLrangeBlock(n1, n2+1, mu_p, mu_p); 
    }
    double r1 = likelihoodRatio( n1-1, mu_p );
    double r2 = likelihoodRatio( n2+1, mu_p );
    if (r1 > r2) {
      return CLrangeBlock(n1-1, n2, mu_p, mu_p); 
    } else {
      return CLrangeBlock(n1, n2+1, mu_p, mu_p); 
    }
  } // case (a) 
  // Actually, control should never fall through to here.
  std::cout << 
    "expandBlockHeight fell through to returnthat should not be reached!!\n";
  return CLrangeBlock(n1, n2+1, block.mu_high, block.mu_high); 
} // expandBlockHeight 


double CFsensitivity::meanMuMax() const 
{
  if (verbosity >= CFV_ULS) {
    std::cout << "\nBlocks are: \n\n";
    for (unsigned int i=0; i<region.size(); ++i) {
      std:: cout << region[i].n1 << "  " 
                 << region[i].n2 << "  " 
                 << region[i].mu_low << "  " 
                 << region[i].mu_high << "\n";
    }
  }
  if (b == 0) {
    bool zeroCountsIsInRegionForBzero;
    double muMax = findUpperMuLimit ( 0, zeroCountsIsInRegionForBzero );
    assert (zeroCountsIsInRegionForBzero);
    return muMax;
  }
  double tolerance = 1.0e-6;
  double psum = 0;
  double pmusum = 0;
  double poissonProbability_n_b = 1.0;
  for (int n = 0; n <= maximum_n; ++n)  
  {
    if (n!=0) poissonProbability_n_b *= b/n;
        // A factor of exp(-b) is in both numerator and denomintor; we ignore it
    bool valid;
    double muMax = findUpperMuLimit ( n, valid );
    if (verbosity >= CFV_BLOCKS) {
      std::cout << "upper mu limit for " << n << " = " << muMax << "\n";
      if (!valid) std::cout << "--- but valid is false\n";
    }
    if (!valid) continue;
    psum += poissonProbability_n_b;
    double newTerm = poissonProbability_n_b * muMax;
    pmusum += newTerm;
    if (verbosity >= CFV_BLOCKS) {
      std::cout << "P(" << n << ", " << b << ") = " 
                << poissonProbability_n_b << "\n"
                << "newTerm = " << newTerm << " --> " << pmusum << "\n";      
    }
    if (newTerm < tolerance) break;
  }
  assert (psum > 0);
  return pmusum/psum;
}

double CFsensitivity::findUpperMuLimit ( int n, bool & valid ) const 
{
  double muMax = 0;
  valid = false;
  for (unsigned int i = 0; i < region.size(); ++i) {
    if ( (region[i].n1 <= n) && (region[i].n2 >= n) ) {
      valid = true;
      muMax = std::max(muMax, region[i].mu_high);
    }
  }
  return muMax;
}


} // end namespace mu2e

void createTable ( std::ostream & os, double p,
                   double b0, int nb, double bInterval )
{
  mu2e::CFsensitivity s(p);
  os.precision(7);
  os << "  double pValue    = " << p         << "; \n"; 
  os << "  double b0        = " << b0        << "; \n";
  os << "  int    nb        = " << nb        << "; \n";
  os << "  double bInterval = " << bInterval << "; \n";
  os << "  double bInterval = " << bInterval << "; \n";
  os << "  double sensitivityValues[" << nb << "] = \n";
  os << "  { \n";
  os << "   " << s(b0);   
  for (int i = 1; i < nb; ++i) {
    os << ",  " << s(b0+i*bInterval);
    if (i%5 == 4) {
      os << "  // " << b0 + (i-4)*bInterval << "\n";
    }
  }
  os << "\n  }; \n";
}

int main ()
{
  int kind;
  std::cout << "Enter kind of app (0 or 1): ";
  std::cin >> kind;
  if (kind==0) {
    double p;
    std::cout << "Enter p: ";
    std::cin >> p; 
    double b;
    std::cout << "Enter b: ";
    std::cin >> b; 
    mu2e::CFsensitivity s(p);
    double x = s(b);
    std::cout << x << "\n";
    std::cout << "Enter b: ";
    std::cin >> b; 
    x = s(b);
    std::cout << x << "\n";
  } else {
    double p;
    std::cout << "Enter p: ";
    std::cin >> p; 
    double b0;
    double nb;
    double bInterval;
    std:: cout << "Enter b0: ";
    std::cin >> b0; 
    std:: cout << "Enter nb: ";
    std::cin >> nb; 
    std:: cout << "Enter bInterval: ";
    std::cin >> bInterval; 
    std::string fileName;
    std::cout << "Enter file name: ";
    std::cin >> fileName;
    if (fileName == "cout") {
      std::ofstream os(fileName.c_str());
      createTable ( std::cout, p, b0, nb, bInterval );  
    } else {
      std::ofstream os(fileName.c_str());
      createTable ( os, p, b0, nb, bInterval );
    }
  }
  return 0;
} 
