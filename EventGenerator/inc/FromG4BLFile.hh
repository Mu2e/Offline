#ifndef FROMG4BLFILE_HH
#define FROMG4BLFILE_HH
//
// Read particles from a file in G4beamline input format.
//
// $Id: FromG4BLFile.hh,v 1.1 2010/08/30 22:50:00 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/30 22:50:00 $
//
// Original author Rob Kutschke
//
// The position is given in the Mu2e coordinate system.
// 

#include <fstream>

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

// External includes
#include "CLHEP/Random/RandPoissonQ.h"

// Forward references.
namespace edm{
  class Run;
}
class TH1F;

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class FromG4BLFile: public GeneratorBase{

  public:
    FromG4BLFile( edm::Run const& run, const SimpleConfig& config );
    virtual ~FromG4BLFile();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Start: Information from the run time configuration.

    // Number of particles per event.
    // If positive, mean of a Poisson distribution.
    // If negative, then exactly that number of particles per event.
    double _mean;

    // The name of the input file.
    std::string _inputFileName;

    // Enable histograms
    bool _doHistograms;

    // End: Information from the run time configuration.

    // Random number distributions.
    CLHEP::RandPoissonQ _randPoissonQ;

    // The input file.
    std::ifstream _inputFile;

    // Histogram information.
    TH1F* _hMultiplicity;
    TH1F* _hMomentum;
    TH1F* _hCz;
    TH1F* _hX0;
    TH1F* _hY0;
    TH1F* _hZ0;
    TH1F* _hT0;

  };

} // end namespace mu2e,

#endif
