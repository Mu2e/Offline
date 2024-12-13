// Ed Callaghan
// Interface for art tool to calculate a single-variable functional prescale
// September 2024

#ifndef Blinding_UnivariateLookupKalSeedPrescaleTool_hh
#define Blinding_UnivariateLookupKalSeedPrescaleTool_hh

// stl
#include <vector>

// boost
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

// mu2e
#include "Offline/Blinding/inc/KalSeedPrescaleTool.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"

using boost::math::interpolators::cardinal_cubic_b_spline;

namespace mu2e{
  class UnivariateLookupKalSeedPrescaleTool: public KalSeedPrescaleTool{
    public:
      UnivariateLookupKalSeedPrescaleTool(double,
                                          double,
                                          std::vector<double>);

    protected:
      cardinal_cubic_b_spline<double> _spline;
      double _xmin;
      double _xmax;
      double calculate_acceptance_rate(const KalSeed&) override;
      virtual double calculate_observable(const KalSeed&) = 0;

    private:
      /**/
  };
} // namespace mu2e

#endif
