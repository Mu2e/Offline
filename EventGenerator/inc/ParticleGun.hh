#ifndef PARTICLEGUN_HH
#define PARTICLEGUN_HH
//
// Shoots a single particle gun and puts its output
// into a generated event.
//
// $Id: ParticleGun.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
// 

#include "EventGenerator/inc/GeneratorBase.hh"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class ParticleGun: public GeneratorBase{

  public:
    ParticleGun( const SimpleConfig& config );
    virtual ~ParticleGun();

    virtual void generate( ToyGenParticleCollection&  );

  private:


  };

} // end namespace mu2e,

#endif


