#ifndef EJECTEDPROTONGUN_HH
#define EJECTEDPROTONGUN_HH

//
// Generate an electron with the conversion energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: EjectedProtonGun.hh,v 1.2 2009/12/30 19:14:19 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2009/12/30 19:14:19 $
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


using CLHEP::RandGeneral;

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

    double EnergyEjectedProtonFunc(const double protonKineticEnergyinMeV); 

    TH1D* _ejectedProtonMultiplicity;
    TH1D* _ejectedProtonKE;
    TH1D* _ejectedProtonKEZoom;
    TH1D* _ejectedProtonMomentumMeV;

    RandomUnitSphere _randomUnitSphere;
    std::auto_ptr<RandGeneral> _funcGen;


    double _mean; //< mean per event
    double _elow; //< lower photon energy 
    double _ehi; //< upper photon energy 
    double _bindE; //< energy bin width for generator
    int _nbins; //< number of bins in photon energy pdf


    // simulation of ejected protons?
    bool _doEjectedProton;

    // ejected proton default total energy.
    double _ejectedProtonEnergy;
    double _ejectedProtonMomentum;
    double _ejectedProtonKineticEnergy;

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


