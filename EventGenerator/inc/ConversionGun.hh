#ifndef CONVERSION_HH
#define CONVERSION_HH

//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: ConversionGun.hh,v 1.4 2010/08/18 06:33:00 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/18 06:33:00 $
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Uniform in time during the requested interval.
//  - Limits on cos(theta) and phi but uniform within the range.
//

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"

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

    // Electron mass
    double _mass;

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

    // A distribution used in many places in this class.
    CLHEP::RandFlat _randFlat;

  };

} // end namespace mu2e,

#endif


