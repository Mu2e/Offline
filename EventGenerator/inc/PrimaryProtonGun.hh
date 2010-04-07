#ifndef PRIMARYPROTONGUN_HH
#define PRIMARYPROTONGUN_HH

//
// Generate a proton with the primary proton energy
// from a random spot within the target system at
// a random time during the accelerator cycle.
//
// $Id: PrimaryProtonGun.hh,v 1.3 2010/04/07 17:43:56 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2010/04/07 17:43:56 $
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

  class PrimaryProtonGun: public GeneratorBase{

  public:
    PrimaryProtonGun( edm::Run& run, const SimpleConfig& config );
    virtual ~PrimaryProtonGun();

    virtual void generate( ToyGenParticleCollection&  );

  private:

    // simulation primary proton gun?
    bool _doPrimaryProt;

    // primary proton momentum.
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


    // primary proton default total energy.
    double _primaryProtonEnergy;
    double _primaryProtonMomentum;
    double _primaryProtonKineticEnergy;


    // primary proton gun coordinates in coordinate system where the center is in the middle of TS
    
    Hep3Vector _beamDisplacementOnTarget; 
    
    double _stdDev;
    double _zOffset;

    
    //histos
    TH1D* _primaryProtonKE;
    TH1D* _primaryProtonKEZoom;
    TH1D* _primaryProtonMomentumMeV;


  };

} // end namespace mu2e,

#endif


