#ifndef CONVERSION_HH
#define CONVERSION_HH

//
// Generate an electron with the conversion energy
// Uses FoilParticleGenerator to extract a random spot 
// within the target system at
// a random time during the accelerator cycle.
//
// $Id: ConversionGun.hh,v 1.6 2010/10/27 16:42:56 onoratog Exp $
// $Author: onoratog $ 
// $Date: 2010/10/27 16:42:56 $
//

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// Forward references in other namespaces.
namespace edm {
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class ConversionGun: public GeneratorBase{

  public:
    ConversionGun( edm::Run& run, const SimpleConfig& config );
    virtual ~ConversionGun();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // simulation conversions?
    bool _doConvs;

    // Conversion momentum.
    double _p;

    //Particle mass
    double _mass;

    // Limits on the generated direction.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    //Utility to generate direction of the momentum
    RandomUnitSphere    _randomUnitSphere;

  };

} // end namespace mu2e,

#endif


