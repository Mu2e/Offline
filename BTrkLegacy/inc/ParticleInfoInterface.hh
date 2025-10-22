#ifndef ParticleInfoInterface_hh
#define ParticleInfoInterface_hh

//
// Definition of an interface for information that must be
// supplied by external code.
//

// FIMXE:
// This generates a compile time circular dependence - it needs to know the type
// of TrkParticle::type.  It does NOT generate a link time or run time circular depennndence.
#include "BTrkLegacy/inc/TrkParticle.hh"

#include <string>

class ParticleInfoInterface {

public:

  virtual double      mass  ( TrkParticle::type ) const = 0;
  virtual double      charge( TrkParticle::type ) const = 0;
  virtual std::string name  ( TrkParticle::type ) const = 0;

protected:
  ParticleInfoInterface(){}
  virtual ~ParticleInfoInterface(){}

};

#endif // end ParticleInfoInterface_hh
