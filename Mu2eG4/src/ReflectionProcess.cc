//
// Still in development.
//
// A discrete process that reverses particle momentum and charge.
// It is used for debugging/tuning code for propagation in magnetic fields.
//
//   $Id: ReflectionProcess.cc,v 1.2 2012/07/26 19:01:01 kutschke Exp $
//   $Author: kutschke $
//   $Date: 2012/07/26 19:01:01 $
//
// Original author R. Bernstein
//

// C++ includes
#include <cstdio>
#include <cmath>

// Framework includes
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/ReflectionProcess.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"

using namespace std;

namespace mu2e {
  G4VParticleChange* ReflectionProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
  {
    // we want to do something in two cases:
    //   a) we're at the last step of the ending PhysicalVolume, in which case we should
    //      reflect
    //   b) we are back where we started, in which case we want to make a bunch of histograms
    //      and report on our success or lack thereof.  We also want to stop the track since it's
    //      pointless to continue

    //
    // this is clumsy but IsFirstStepInVolume doesn't seem to work but this does.  If
    // I don't do this, after the track is reflected, the vertex position becomes where the track was reflected from...


    if (aTrack.GetVolume()->GetName() == _startingVolume && aTrack.GetCurrentStepNumber() == 1)
      {
        //
        // need to save where track starts for future reference
        startingVertex = aTrack.GetVertexPosition();
        //
        // and set reflection flag.
        alreadyReflected = false;
      }

    if (aTrack.GetVolume()->GetName() == _endingVolume && !alreadyReflected )
      {
        alreadyReflected = true;
        return ReflectIt(aTrack,aStep);
      }
    else if (aTrack.GetVolume()->GetName() == _startingVolume && alreadyReflected
             && (abs (startingVertex.z() - aTrack.GetPosition().z()) < _toleranceForQuitting ) )
      {
        G4cout << "z difference at KillIt: " << startingVertex.z() << " " << aTrack.GetPosition().z() << " " <<
          abs(startingVertex.z() - aTrack.GetPosition().z())<< G4endl;
        double zDistance = abs(startingVertex.z() - aTrack.GetPosition().z());
        double totalDistanceAtKill = Distance(startingVertex,aTrack.GetPosition());
        //  - square(zDistance) );
        //          double totalDistanceAtKill = safeSqrt(square(Distance(startingVertex,aTrack.GetPosition()))
        //  - square(zDistance) );
        double transverseDistanceAtKill = safeSqrt(totalDistanceAtKill*totalDistanceAtKill - zDistance*zDistance);
        G4cout << "and transverse distance is = " << transverseDistanceAtKill << G4endl;

        //
        // and look at direction cosines

        return KillIt(aTrack,aStep);
      }
    else
      {
        return DoNothing(aTrack,aStep);
      }
  }

  G4VParticleChange* ReflectionProcess::ReflectIt( const G4Track& aTrack,
                                                const G4Step& aStep){
    //
    // reflect momentum and reverse charge since you hit the end of the ending volume

    //
    // first, initialize particle change; all members of G4VParticleChange are set
    // equal to corresponding member in G4Track

    fMu2eParticleChangeForReflection.Initialize(aTrack);

    // set the number of secondaries to 1, i.e. just the reflected particle
    fMu2eParticleChangeForReflection.SetNumberOfSecondaries(1);

    //
    // below is based on Geant4/processes/decay/src/G4UnknownDecay.cc
    G4ThreeVector currentPosition = aTrack.GetPosition();
    G4double currentKineticEnergy = aTrack.GetKineticEnergy();
    G4ThreeVector currentMomentumDirection = aTrack.GetMomentumDirection();
    G4double finalGlobalTime = aTrack.GetGlobalTime();

    const G4TouchableHandle thand = aTrack.GetTouchableHandle();
    //
    //find antiparticle
    G4ParticleDefinition* particleType = aTrack.GetDefinition();
    G4int antiParticlePDGEncoding = - ( particleType->GetPDGEncoding() );
    //
    //FindAntiParticle is not coded, I think; see
    //  http://www-geant4.kek.jp/lxr/source/particles/management/include/G4ParticleTable.hh
    //and that method is not hyperlinked, hence it probably doesn't exist, so do it this way

    G4ParticleDefinition* antiParticle = G4ParticleTable::GetParticleTable()->FindParticle(antiParticlePDGEncoding);
    const G4ThreeVector& startingMomentum = G4ThreeVector();
    G4DynamicParticle* reflectedParticle = new G4DynamicParticle(antiParticle,0.,startingMomentum);
    //
    //and make our reflected particle into a track
    reflectedParticle->SetMomentumDirection(-currentMomentumDirection);
    reflectedParticle->SetKineticEnergy(currentKineticEnergy);
    //    reflectedParticle->DumpInfo();
    G4Track* secondary = new G4Track ( reflectedParticle,finalGlobalTime,currentPosition);
    secondary->SetGoodForTrackingFlag();
    secondary->SetTouchableHandle(thand);
    fMu2eParticleChangeForReflection.AddSecondary(secondary);

    //
    //kill parent
    fMu2eParticleChangeForReflection.ProposeTrackStatus( fStopAndKill );
    fMu2eParticleChangeForReflection.ProposeLocalEnergyDeposit( 0. );
    fMu2eParticleChangeForReflection.ProposeGlobalTime( finalGlobalTime );

    //
    // don't need to clear number of interaction lengths, since this class is only used with physics processes off....
    //delete reflectedParticle;
    return &fMu2eParticleChangeForReflection;
  }


  G4VParticleChange* ReflectionProcess::KillIt(const G4Track& aTrack,const G4Step& aStep){
    //
    //you're back where you started, so stop

    //do I need to do this initialization??
    G4cout << "inside KillIt" << G4endl;
    //     fMu2eParticleChangeForReflection.Initialize(aTrack);
    fMu2eParticleChangeForReflection.ProposeTrackStatus( fStopAndKill );
    return &fMu2eParticleChangeForReflection;
  }


  G4VParticleChange* ReflectionProcess::DoNothing( const G4Track& aTrack,
                                                const G4Step& aStep){
    fMu2eParticleChangeForReflection.Initialize(aTrack); // need to return a proper object
    return &fMu2eParticleChangeForReflection;
  }

  G4double ReflectionProcess::GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
                                           G4ForceCondition* condition)
  {

    *condition = Forced;
    return DBL_MAX;
  }

  G4bool ReflectionProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
  {
    return (aParticleType.GetPDGCharge() != 0);
  }
} //end namespace mu2e
