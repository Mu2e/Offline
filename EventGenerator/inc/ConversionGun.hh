#ifndef CONVERSION_HH
#define CONVERSION_HH

//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: ConversionGun.hh,v 1.2 2009/10/16 04:20:52 shanahan Exp $
// $Author: shanahan $ 
// $Date: 2009/10/16 04:20:52 $
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Uniform in time during the requested interval.
//  - Limits on cos(theta) and phi but uniform within the range.
//

#include "EventGenerator/inc/GeneratorBase.hh"

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

    // Limits on the generated direction.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    // Range for the above.
    double _dcz;
    double _dphi;
    double _dt;


  };

} // end namespace mu2e,

#endif


