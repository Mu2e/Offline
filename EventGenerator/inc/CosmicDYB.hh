#ifndef CosmicDYB_HH
#define CosmicDYB_HH
//
// Muon generator, uses Daya Bay libraries
//
// $Id: CosmicDYB.hh,v 1.4 2010/08/18 22:40:15 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/08/18 22:40:15 $
//
// Original author Yury Kolomensky
//

#include <vector>

#include "EventGenerator/inc/GeneratorBase.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

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

    // Do we make the histograms.
    bool _doHistograms;

    // end of configurable parameters

    // Time range (in ns) over which to generate events.
    double _tmin;
    double _tmax;
    double _dt;

    // Random number distributions.
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;

    // Default dimensions
    static const long _default_ne = 4000;
    static const long _default_nth = 100;

    // Working space for hrndg2 ( working space will be on the heap).
    std::vector<double> _workingSpace;


  };  // CosmicDYB

}  // namespace mu2e

#endif  // CosmicDYB_HH
