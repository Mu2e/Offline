#ifndef CosmicDYB_HH
#define CosmicDYB_HH
//
// Muon generator, uses Daya Bay libraries
//
// $Id: CosmicDYB.hh,v 1.3 2010/08/17 15:18:39 wb Exp $
// $Author: wb $
// $Date: 2010/08/17 15:18:39 $
//
// Original author Yury Kolomensky
//

#include <vector>

#include "EventGenerator/inc/GeneratorBase.hh"

// Forward declarations.
class TH1D;
class TH2D;

namespace edm{
  class Run;
}

namespace mu2e {

  // Forward declaratations.
  class SimpleConfig;

  class CosmicDYB: public GeneratorBase {

  public:
    CosmicDYB( edm::Run& run, const SimpleConfig& config );
    virtual ~CosmicDYB();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Histograms
    TH1D* _cosmicMultiplicityH;
    TH1D* _cosmicMomH;
    TH1D* _cosmicChargeH;
    TH1D* _cosmicCosThetaH;
    TH2D* _cosmicCosThetaVsEH;

    // configurable parameters

    // Mean multiplicity. If negative, use -_mean as a fixed number
    double _mean;
    double _muEMin;   // min and max values of muon energy (MeV)
    double _muEMax;
    double _muCosThMin; // min and max zenith angles
    double _muCosThMax;

    // half area to generate events (cm)
    double _dx;
    double _dz;
    double _y0;  // reference depth

    // Dimensions of the 2d working space for hrndg2.
    long _ne;
    long _nth;

    // end of configurable parameters

    // Time range (in ns) over which to generate events.
    double _tmin;
    double _tmax;
    double _dt;

    // Default dimensions
    static const long _default_ne = 4000;
    static const long _default_nth = 100;

    // Working space for hrndg2 ( working space will be on the heap).
    std::vector<double> _workingSpace;

  };  // CosmicDYB

}  // namespace mu2e

#endif  // CosmicDYB_HH
