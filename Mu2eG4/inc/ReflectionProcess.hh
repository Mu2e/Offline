#ifndef Mu2eG4_ReflectionProcess_hh
#define Mu2eG4_ReflectionProcess_hh
//
// Still in development.
//
// A discrete process that reverses particle momentum and charge.
// It is used for debugging/tuning code for propagation in magnetic fields.
//
//   $Id: ReflectionProcess.hh,v 1.2 2013/10/25 21:46:26 genser Exp $
//   $Author: genser $
//   $Date: 2013/10/25 21:46:26 $
//
// Original author R. Bernstein

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Geant4 includes
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"
#include "cetlib/pow.h"

namespace mu2e{

  class ReflectionProcess : public G4VDiscreteProcess{
  public:
    ReflectionProcess(
                      const G4String& startingVolume, //need to do the string; G4 may not have
                      const G4String& endingVolume,   //made the volume yet
                      const double toleranceForQuitting = 0.001*CLHEP::meter):
      G4VDiscreteProcess("ReflectionProcess",fUserDefined),
      _startingVolume(startingVolume),
      _endingVolume(endingVolume),
      _toleranceForQuitting(toleranceForQuitting),
      alreadyReflected(false)
    {
      //tell the user what I'm doing; could iterate over all volumes and check names exist
      G4cout << "Reflecting through following volumes"
             << G4endl
             << "starting Volume was: " <<startingVolume
             << "\n"
             << "ending Volume was:   " <<endingVolume
             << G4endl;
      pParticleChange = &fMu2eParticleChangeForReflection;

    }
    ~ReflectionProcess(){}

  public:
    const G4String _startingVolume;
    const G4String _endingVolume;
    double _toleranceForQuitting;
    G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
    G4ParticleTable* particleTable;


  protected:
    // GetMeanFreePath returns ctau*beta*gamma for decay in flight
    virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                     G4double   previousStepSize,
                                     G4ForceCondition* condition
                                     );

  private:

    G4ParticleChange fMu2eParticleChangeForReflection;

    G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

    G4VParticleChange* ReflectIt (const G4Track& aTrack, const G4Step& aStep);
    G4VParticleChange* KillIt    (const G4Track& aTrack, const G4Step& aStep);
    G4VParticleChange* DoNothing (const G4Track& aTrack, const G4Step& aStep);

    ReflectionProcess( const ReflectionProcess& right);
    ReflectionProcess & operator=(const ReflectionProcess &right);

    bool alreadyReflected;

    G4ThreeVector startingVertex;

    double Distance(const G4ThreeVector a, const G4ThreeVector b)
    {
      return std::sqrt( cet::sum_of_squares( a.x() - b.x(),
                                             a.y() - b.y(),
                                             a.z() - b.z() ) );
    }

  };


} //end namespace mu2e

#endif /* Mu2eG4_ReflectionProcess_hh */
