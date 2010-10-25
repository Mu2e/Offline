#ifndef EJECTEDPROTONGUN_HH
#define EJECTEDPROTONGUN_HH

//
// Simulate the protons that come from the stopping target when muons capture
// on an Al nucleus.  Use the MECO distribution for the kinetic energy of the
// protons.  Production is uniform across the targets and uniform in time;
// this model needs to be improved.
//
// $Id: EjectedProtonGun.hh,v 1.5 2010/10/25 21:12:44 onoratog Exp $
// $Author: onoratog $ 
// $Date: 2010/10/25 21:12:44 $
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Uniform in time during the requested interval.
//  - Limits on cos(theta) and phi but uniform within the range.
//

// Framework includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "EventGenerator/inc/FoilParticleGenerator.hh"

// CLHEP includes
#include "CLHEP/Random/RandPoissonQ.h"

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


    // Class for generate particles from target
    FoilParticleGenerator fGenerator;


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

    CLHEP::RandPoissonQ _randPoissonQ;


    TH1D* _hMultiplicity;
    TH1D* _hKE;
    TH1D* _hKEZoom;
    TH1D* _hMomentumMeV;
    TH1D* _hzPosition;
    TH1D* _hcz;
    TH1D* _hphi;
    TH1D* _htime;

  };

} // end namespace mu2e,

#endif
