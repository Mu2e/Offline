#ifndef DECAYINORBIT_HH
#define DECAYINORBIT_HH

//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: DecayInOrbitGun.hh,v 1.2 2010/05/17 21:47:33 genser Exp $
// $Author: genser $ 
// $Date: 2010/05/17 21:47:33 $
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Uniform in time during the requested interval.
//  - Limits on cos(theta) and phi but uniform within the range.
//

#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "CLHEP/Random/RandGeneral.h"

// Framework Includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

//ROOT Includes
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"



class TH1D;

namespace edm {
  class Run;
}

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class DecayInOrbitGun: public GeneratorBase{

  public:
    DecayInOrbitGun( edm::Run& run, const SimpleConfig& config );
    virtual ~DecayInOrbitGun();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    const double EnergyDIOFunc(double x);


    TH1D* _decayInOrbitMultiplicity;
    TH1D* _decayInOrbitEElec;
    TH1D* _decayInOrbitEElecZ;

    RandomUnitSphere _randomUnitSphere;
    std::auto_ptr<CLHEP::RandGeneral> _funcGen;

    double _mean; //< mean per event
    double _elow; //< lower photon energy 
    double _ehi; //< upper photon energy 
    double _bindE; //< energy bin width for generator
    int _nbins; //< number of bins in photon energy pdf


    // simulation conversions?
    bool _doConvs;

    // Conversion momentum.
    double _p;

    // Limits on the generated direction.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    // Range for the above.
    double _dcz;
    double _dphi;
    double _dt;


  };

} // end namespace mu2e,

#endif


