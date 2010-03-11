#ifndef CosmicDYB_HH
#define CosmicDYB_HH
//
// Muon generator, uses Daya Bay libraries
//
// $Id: CosmicDYB.hh,v 1.1 2010/03/11 02:07:42 yury Exp $
// $Author: yury $ 
// $Date: 2010/03/11 02:07:42 $
//
// Original author Yury Kolomensky
//

#include "EventGenerator/inc/GeneratorBase.hh"

class TH1D;

namespace edm{
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class CosmicDYB: public GeneratorBase{

  public:
    CosmicDYB( edm::Run& run, const SimpleConfig& config );
    virtual ~CosmicDYB();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Histograms
    TH1D* _cosmicMultiplicityH;
    TH1D* _cosmicMomH;
    TH1D* _cosmicChargeH;

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

    // Time range ( in ns) over which to generate events.
    double _tmin;
    double _tmax;
    double _dt;

    // buffer for DYB routines
    double _buffer[100][4000];

  };

} // end namespace mu2e,

#endif


