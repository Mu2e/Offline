#ifndef EventGenerator_CosmicFromTH2_hh
#define EventGenerator_CosmicFromTH2_hh
//
// Cosmic ray muon generator using a TH2 as probability distribution
//
//
// Original author Ralf Ehrlich
//

#include <vector>

#include "EventGenerator/inc/GeneratorBase.hh"

#include "CLHEP/Random/RandFlat.h"

// Forward declarations.
class TH2;
class TFile;
namespace art{
  class Run;
}

namespace mu2e {

  // Forward declaratations.
  class SimpleConfig;

  class CosmicFromTH2: public GeneratorBase {

  public:
    CosmicFromTH2(CLHEP::HepRandomEngine& engine,
                  art::Run& run,
                  const SimpleConfig& config);
    virtual ~CosmicFromTH2();

    virtual void generate( GenParticleCollection&  );

  private:

    TFile *_file;
    TH2   *_histogram;
    double _energy;
    double _p;
    double _time;

    CLHEP::Hep3Vector _cosmicReferencePointInMu2e;
    double _dx;
    double _dy;
    double _dz;

    CLHEP::RandFlat _randFlat;

    bool _createdProductionPlane;
  };  // CosmicDYB

}  // namespace mu2e

#endif /* EventGenerator_CosmicFromTH2_hh */
