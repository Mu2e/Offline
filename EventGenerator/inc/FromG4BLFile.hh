#ifndef EventGenerator_FromG4BLFile_hh
#define EventGenerator_FromG4BLFile_hh
//
// Read particles from a file in G4beamline input format.
//
// $Id: FromG4BLFile.hh,v 1.8 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:35 $
//
// Original author Rob Kutschke
//
// The position is given in the Mu2e coordinate system.
// 

#include <fstream>

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "ToyDP/inc/G4BeamlineInfoCollection.hh"

// External includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandPoissonQ.h"

// Forward references.
namespace art{
  class Run;
}
class TH1F;
class TNtuple;

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class FromG4BLFile: public GeneratorBase{

  public:
    FromG4BLFile( art::Run const& run, const SimpleConfig& config );
    virtual ~FromG4BLFile();

    virtual void generate( ToyGenParticleCollection&  );
    void generate( ToyGenParticleCollection& , G4BeamlineInfoCollection*  );

  private:

    // Start: Information from the run time configuration.

    // Number of particles per event.
    // If positive, mean of a Poisson distribution.
    // If negative, then exactly that number of particles per event.
    double _mean;

    // The midpoint of the target in the coordinates used in the input file.
    CLHEP::Hep3Vector _prodTargetOffset;

    // The center of the production target, in the Mu2e coordinate system.
    CLHEP::Hep3Vector _prodTargetCenter;

    // The origin of the G4beamline coordinate system, specified in the Mu2e system.
    CLHEP::Hep3Vector _g4beamlineOrigin;

    // An offset relative to the g4beamlineOrigin - used to fix minor discrepancies
    // in G4beamline vs Offline geometries.
    CLHEP::Hep3Vector _g4beamlineExtraOffset;

    // The name of the input file.
    std::string _inputFileName;

    // List of pdg Id's in which we are interested.  If empty keep all pdgIds.
    std::vector<int> _pdgIdToKeep;

    // Enable histograms
    bool _doHistograms;

    // Is target frame used in the file?
    bool _targetFrame;

    //Number of particles to skip form the input file.
    //Useful to run grid-jobs reading different segments of the same txt file
    int _nPartToSkip;

    // End: Information from the run time configuration.

    // Random number distributions.
    CLHEP::RandPoissonQ _randPoissonQ;

    // The input file.
    std::ifstream _inputFile;

    // Histogram and ntuple information.
    TH1F* _hMultiplicity;
    TH1F* _hMomentum;
    TH1F* _hCz;
    TH1F* _hX0;
    TH1F* _hY0;
    TH1F* _hZ0;
    TH1F* _hT0;
    TNtuple* _ntup;

  };

} // end namespace mu2e,

#endif /* EventGenerator_FromG4BLFile_hh */
