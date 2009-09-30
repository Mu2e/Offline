#ifndef GENERATORBASE_HH
#define GENERATORBASE_HH
//
//
// Base class to allow generic access to all of the 
// event generator classes.
//
// $Id: GeneratorBase.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
// 

#include "ToyDP/inc/ToyGenParticleCollection.hh"

namespace mu2e {

  class GeneratorBase{

  public:
    GeneratorBase( ){
    }
    virtual ~GeneratorBase(){
    }

    virtual void generate( ToyGenParticleCollection&  ) = 0;

};

} // end namespace mu2e,

#endif

