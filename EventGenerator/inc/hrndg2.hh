#ifndef EventGenerator_hrndg2_hh
#define EventGenerator_hrndg2_hh
//
// Interface to the fortran subroutine hrndg2.
//
// $Id: hrndg2.hh,v 1.3 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:35 $
//
// Original author Rob Kutschke
// 
// Notes:
// 1) Inside hrndg2, the working space is:
//     double work[nth][ne];
//    In the default configuration this working space is 3.2 MB
//    so we want it on the heap, not on the stack.  The use of
//    the vector is just an convenient way to do this.
//    The vector must have memory allocated for nth*ne elements.
//

#include <vector>

// The interface to hrndg2
extern "C" {
  void hrndg2_( double*, 
                const long*   const,
                const double* const,
                const double* const,
                const long*   const,
                const double* const,
                const double* const,
                double*,
                double*,
                double*,
                const float* const);

}

namespace mu2e {

  // A convenience function so that user code reads a little cleaner.
  inline void hrndg2( std::vector<double>& work,
                      const long&          ne,
                      const double&        eLow,
                      const double&        eHigh,
                      const long&          nth,
                      const double&        thLow,
                      const double&        thHigh,
                      double&              dimSum,
                      double&              e,
                      double&              th,
                      const float&         pro){
    ::hrndg2_( &work[0], &ne, &eLow, &eHigh, &nth, &thLow, &thHigh, &dimSum, &e, &th, &pro);
  }

}

#endif /* EventGenerator_hrndg2_hh */
