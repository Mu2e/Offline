#ifndef EJECTEDPROTONGUN_HH
#define EJECTEDPROTONGUN_HH

//
// Simulate the protons that come from the stopping target when muons capture
// on an Al nucleus.  Use the MECO distribution for the kinetic energy of the
// protons.  Production is uniform across the targets and uniform in time;
// this model needs to be improved.
//
// $Id: EjectedProtonGun.hh,v 1.4 2010/08/19 22:03:11 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/19 22:03:11 $
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Uniform in time during the requested interval.
//  - Limits on cos(theta) and phi but uniform within the range.
//

// Framework includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

// Forward declarations outside of namespace mu2e
class TH1D;
namespace edm {
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class EjectedProtonGun: public GeneratorBase{

  public:
    EjectedProtonGun( edm::Run& run, const SimpleConfig& config );
    virtual ~EjectedProtonGun();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // Start: parameters that can be configured from the config file.

    double _mean;    // Mean number of protons per event
    double _elow;    // Range of proton energy.
    double _ehi;     // 
    double _czmin;   // Range of cos(polar angle)
    double _czmax;  
    double _phimin;  // Range of azimuth
    double _phimax;
    int    _nbins;   // number of bins in proton energy pdf

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    // Histogram control.
    bool _doHistograms;

    // end: parameters that can be configured from the config file.

    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;
    RandomUnitSphere    _randomUnitSphere;
    CLHEP::RandGeneral  _shape;

    // Mass of the proton (taken from particle data table).
    double _mProton;

    TH1D* _hMultiplicity;
    TH1D* _hKE;
    TH1D* _hKEZoom;
    TH1D* _hMomentumMeV;
    TH1D* _hzPosition;
    TH1D* _hcz;
    TH1D* _hphi;
    TH1D* _htime;

    // Compute the value of the energy spectrum at given energy.
    double energySpectrum(double e);

    // Build a binned representation of the energy spectrum.
    std::vector<double> binnedEnergySpectrum();

  };

} // end namespace mu2e,

#endif
