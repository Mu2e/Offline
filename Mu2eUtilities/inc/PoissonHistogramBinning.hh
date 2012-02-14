#ifndef Mu2eUtilities_PoissonHistogramBinning_hh
#define Mu2eUtilities_PoissonHistogramBinning_hh
//
// Choose a suitable range for histogramming a random variate distributed
// according to a mu2e::PoissonOrFixed distribution.
//
// $Id: PoissonHistogramBinning.hh,v 1.1 2012/02/14 17:09:54 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/02/14 17:09:54 $
//
// Contact author Rob Kutschke
//
// The mu2e::PoissonOrFixed distribution is constructed:
//
//   PoissonOrFixed::PoissonOrFixed( double mean );
//
// If (mean >0 ) then this distribution behaves as a Poisson distribution
// with the given mean.  If (mean <= 0 ) then the distribution returns
// a fixed value for all calls:
//
//   int n0 = static_cast<int>(floor(std::abs(mean)));
//
// Notes:
//
// 1) The width argument is only used if the mean is more than 70.
//    For this case, the initial guess at the histogram range is taken
//    to be an interval of +/- width*sqrt(mean) around the mean.  This
//    guess is tweaked give, where possible and a nice choice for bin
//    boundaries.
//

#include <ostream>

namespace mu2e {

  class PoissonHistogramBinning {

  public:

    PoissonHistogramBinning( double mean, double width=4. );

    // Accept compiler written d'tor, copy c'tor and copy assignment.

    int   nbins()  const { return nbins_; }
    double xlow()  const { return xlow_;  }
    double xhigh() const { return xhigh_; }
    double mean()  const { return mean_;  }
    double width() const { return width_; }

  private:

    // Copies of input arguments.
    double mean_;
    double width_;

    // Computed quantities.
    int    nbins_;
    double xlow_;
    double xhigh_;

  }; // end PoissonHistogramBinning

  inline std::ostream& operator<<(std::ostream& ost,
                                  PoissonHistogramBinning const& binning ){
    ost << "Binning for mean: " << binning.mean()
        << " is: "
        << binning.nbins() << " "
        << binning.xlow() << " "
        << binning.xhigh() << " ";
    return ost;
  }

} // end namespace mu2e

#endif
