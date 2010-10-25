#ifndef CONVERSION_HH
#define CONVERSION_HH

//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: ConversionGun.hh,v 1.5 2010/10/25 19:50:21 onoratog Exp $
// $Author: onoratog $ 
// $Date: 2010/10/25 19:50:21 $
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Uniform in time during the requested interval.
//  - Limits on cos(theta) and phi but uniform within the range.
//

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "EventGenerator/inc/FoilParticleGenerator.hh"

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

    // Generator of particle from target
    FoilParticleGenerator fGenerator;


    // simulation conversions?
    bool _doConvs;

    // Conversion momentum.
    double _p;

    // Limits on the generated direction.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    // Limits on the generated time.
    double _tmin;
    double _tmax;

  };

} // end namespace mu2e,

#endif


