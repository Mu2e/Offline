//
// Class Description:
//
// G4 class that takes produces a conversion electron from an at rest
// muon. Configurable parameters include endpoint of momentum, limits of polar
// and azimuthal angle.
//
// $Id: MuonMinusConversionAtRest.cc,v 1.6 2012/07/15 22:06:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:17 $
//
//-----------------------------------------------------------------------------
// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/MuonMinusConversionAtRest.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// Geant includes
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4GHEKinematicsVector.hh"
#include "G4HadronicProcessStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


namespace mu2e {

static const double pEndPoint = 104.96;

  muMinusConversionAtRest::muMinusConversionAtRest( const fhicl::ParameterSet& config, const G4String& processName, G4ProcessType aType )
    : G4VRestProcess (processName, aType)
  {
    throw cet::exception("CONFIG")<<"muMinusConversionAtRest(ParameterSet,...) is not implemented\n";
  }

muMinusConversionAtRest::muMinusConversionAtRest( const SimpleConfig& config, const G4String& processName, G4ProcessType aType ) :
    G4VRestProcess (processName, aType),
  _p(config.getDouble("conversionGun.p", pEndPoint )),
  _czmin(config.getDouble("conversionGun.czmin",  0.3)),
  _czmax(config.getDouble("conversionGun.czmax",  0.6)),
  _phimin(config.getDouble("conversionGun.phimin", 0. )),
  _phimax(config.getDouble("conversionGun.phimax", CLHEP::twopi )),
  _randomUnitSphere ( *CLHEP::HepRandom::getTheEngine(), _czmin, _czmax, _phimin, _phimax ),
  isInitialised(false)
{
  SetProcessSubType(fHadronAtRest);
  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

muMinusConversionAtRest::~muMinusConversionAtRest()
{
  G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool muMinusConversionAtRest::IsApplicable(const G4ParticleDefinition& p)
{
  return ( &p == G4MuonMinus::MuonMinus() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void muMinusConversionAtRest::PreparePhysicsTable(const G4ParticleDefinition& p)
{
  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this, &p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void muMinusConversionAtRest::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4VParticleChange* muMinusConversionAtRest::AtRestDoIt(const G4Track& track,
                                                        const G4Step&)
{
  //
  // Handles MuonMinus at rest;
  //
  aParticleChange.Initialize(track);
  G4ThreeVector position = track.GetPosition();
  G4double globaltime = track.GetGlobalTime();

  //Pick up momentum vector
  G4ThreeVector p3 = G4ThreeVector(_randomUnitSphere.fire(_p));

  aParticleChange.SetNumberOfSecondaries( 1 );

  G4ParticleDefinition *theElectron = G4Electron::Electron();
  G4double Emass = theElectron->GetPDGMass();

  G4GHEKinematicsVector ConvElectron;

  ConvElectron.SetZero();
  ConvElectron.SetMass( Emass );
  ConvElectron.SetMomentumAndUpdate(p3.x(), p3.y(), p3.z());
  ConvElectron.SetParticleDef( theElectron );

  G4ParticleDefinition* pd = ConvElectron.GetParticleDef();
  if(pd) {

    // G4 takes ownership of this object
    G4DynamicParticle* aNewParticle = new G4DynamicParticle;
    aNewParticle->SetDefinition( pd );
    aNewParticle->SetMomentum( ConvElectron.GetMomentum() );

    // G4 takes ownership of this object
    G4Track* aNewTrack = new G4Track( aNewParticle, globaltime, position );
            aNewTrack->SetTouchableHandle(track.GetTouchableHandle());
        aParticleChange.AddSecondary( aNewTrack );
    }

  aParticleChange.ProposeLocalEnergyDeposit(0.0);
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  return &aParticleChange;
}

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
