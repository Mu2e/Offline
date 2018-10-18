// Wrapper around art random number engine that satisfies C++ URBG standard.
// Extracted from Dave Brown's ProtonBunchIntensitySimulator_module.cc

#ifndef Mu2eUtilities_inc_artURBG_hh
#define Mu2eUtilities_inc_artURBG_hh

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// The framework code above depends on CLHEP but does not include the header
// so that we have to.
#include "CLHEP/Random/RandomEngine.h"

namespace mu2e {

  class artURBG {
  public:
    artURBG(art::RandomNumberGenerator::base_engine_t& engine) : _engine(engine) {}
    typedef unsigned int result_type;
    result_type min() { return 0; }
    result_type max() { return UINT_MAX; }
    result_type operator() () { return _engine.operator unsigned int(); }
  private:
    art::RandomNumberGenerator::base_engine_t& _engine; // CLHEP ENGINE
  };

}

#endif/*Mu2eUtilities_inc_artURBG_hh*/
