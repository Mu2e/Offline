#ifndef PiCapture_HH
#define PiCapture_HH
//
//
// Generate photons from pi- capture on Al nuclei.
// Based on Ivano Sarra's work described in Mu2e doc 665-v2
// 
// $Id: PiCapture.hh,v 1.8 2010/09/29 22:59:59 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/09/29 22:59:59 $
//
// Original author Rob Kutschke, P. Shanahan
// 

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandExponential.h"

// Forward declarations outside of namespace mu2e.
class TH1D;
namespace edm{
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class PiCapture: public GeneratorBase{

  public:
    PiCapture( edm::Run& run, const SimpleConfig& config );
    virtual ~PiCapture();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Start: parameters from run time configuration/

    double _mean;         //< mean per event
    double _elow;         //< lower photon energy 
    double _ehi;          //< upper photon energy 
    int    _nbins;        //< number of bins in photon energy pdf
    bool   _doHistograms; // Enable/disable histograms.

    // End: parameters from run time configuration/

    // Random number distributions
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;
    CLHEP::RandExponential _randExponential;
    RandomUnitSphere    _randomUnitSphere;
    CLHEP::RandGeneral  _spectrum;

    //
    // for exponential foil dropoff
    int nFoils;
    double foilMean;

    // Histograms.
    TH1D* _hMultiplicity;
    TH1D* _hEPhot;
    TH1D* _hEPhotZ;
    TH1D* _hzPos;
    TH1D* _hcz;
    TH1D* _hphi;
    TH1D* _ht;
    TH1D* _hFoilNumber;

    // Photon energy spectrum as a continuous function.
    const double energySpectrum(const double e);

    // Compute a binned representation of the photon energy spectrum.
    std::vector<double> PiCapture::binnedEnergySpectrum();

  };

} // end namespace mu2e,

#endif
