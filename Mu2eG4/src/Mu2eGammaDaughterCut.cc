// Special process to kill events with low energy photon daughters
//
// Original author M. MacKenzie

// Mu2e includes
#include "Mu2eG4/inc/Mu2eGammaDaughterCut.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"

// G4 includes
#include "G4ios.hh"
#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"

namespace mu2e{

  Mu2eGammaDaughterCut::Mu2eGammaDaughterCut(const G4double minDaughterEnergy, const G4bool killAfterConvert,
                                             G4int verbose,
                                             const G4String& aName)
    : G4VDiscreteProcess(aName,fUserDefined), minDaughterEnergy_(minDaughterEnergy),
      killAfterConvert_(killAfterConvert),
      verbose_(verbose), accepted_(0), photonEnergy_(-1.)
  {
    if(verbose_ > 0) {
      G4cout << GetProcessName() << " is created " << G4endl;
    }
  }

  Mu2eGammaDaughterCut::~Mu2eGammaDaughterCut()
  {}

  G4bool Mu2eGammaDaughterCut::IsApplicable(const G4ParticleDefinition& particle) {
    bool retval = (particle.GetParticleName() == "gamma" ||
                   particle.GetParticleName() == "e+" ||
                   particle.GetParticleName() == "e-");
    if(verbose_ > 1 && retval) G4cout << "Mu2eGammaDaughterCut::" << __func__ << ": Adding particle: "
                                      << particle.GetParticleName()
                                      << G4endl;
    return retval;
  }

  G4double Mu2eGammaDaughterCut::PostStepGetPhysicalInteractionLength( const G4Track& track,
                                                                           G4double,
                                                                           G4ForceCondition* condition) {

    *condition = NotForced;

    if(track.GetTrackID() == 1 && track.GetCurrentStepNumber() == 1) { // new event
      accepted_ = 0; //reset accepted flag
      photonEnergy_ = track.GetTotalEnergy();
    }


    //check if the event has been marked to accept or kill
    if(accepted_ < 0) { //kill event
      return 0.0;
    } else if(accepted_ > 0) { //accepted event
      return DBL_MAX;
    }

    if(verbose_ > 9) G4cout << "Mu2eGammaDaughterCut::" << __func__ << ": Track seen: ID = "
                            << track.GetTrackID() << " Parent ID = " << track.GetParentID()
                            << " step number = " << track.GetCurrentStepNumber()
                            << " CreationCode = " << Mu2eG4UserHelpers::findCreationCode(&track)
                            << " E = " << track.GetTotalEnergy()
                            << " E_gamma = " << photonEnergy_ << " accepted = " << accepted_
                            << G4endl;

    if(track.GetTrackID() == 1) { //generated photon, should always be the first track in the event
      double energy = track.GetTotalEnergy();
      //update the photon energy if no energy is stored or not too small of an energy
      if(energy > 1.) {
        photonEnergy_ = energy;
      }
      if(verbose_ > 1) G4cout << "Mu2eGammaDaughterCut::" << __func__ << ": Updating photon energy" << G4endl;
      return DBL_MAX;
    }


    if(track.GetParentID() == 1) { //daughter of generated photon
      double energy = track.GetTotalEnergy();
      bool isCompt = Mu2eG4UserHelpers::findCreationCode(&track) == ProcessCode(ProcessCode::compt);
      bool isConvn = Mu2eG4UserHelpers::findCreationCode(&track) == ProcessCode(ProcessCode::conv);
      if(!isCompt && !isConvn) {
        if(verbose_ > 1) G4cout << "Mu2eGammaDaughterCut::" << __func__ << ": Not a relevant process" << G4endl;
        return DBL_MAX;
      }
      if(energy > minDaughterEnergy_) {
        if(verbose_ > 0) G4cout << "Mu2eGammaDaughterCut::" << __func__ << ": Event passed the test" << G4endl;
        accepted_ = (isConvn&&killAfterConvert_) ? -1 : 1; //kill the event if a conversion and asked to
        return DBL_MAX;
      } else if(photonEnergy_ > 0. && photonEnergy_ - energy < minDaughterEnergy_) {
        accepted_ = -1; //fails
        if(verbose_ > 0) G4cout << "Mu2eGammaDaughterCut::" << __func__ << ": Event failed the test" << G4endl;
        return 0.0;
      }
      if(verbose_ > 1) G4cout << "Mu2eGammaDaughterCut::" << __func__ << ": Event continues through test" << G4endl;
      return DBL_MAX; //still could have produced a passing daughter --> continue
    }
    if(verbose_ > 9) G4cout << "Mu2eGammaDaughterCut::" << __func__ << ": Event not tested" << G4endl;
    return DBL_MAX; //not the primary photon or a daughter of it
  }

  G4double Mu2eGammaDaughterCut::GetMeanFreePath(const G4Track&,G4double,
                                                     G4ForceCondition*){
    return DBL_MAX;
  }

  G4VParticleChange* Mu2eGammaDaughterCut::PostStepDoIt(const G4Track& trk,
                                                            const G4Step& step) {
    pParticleChange->Initialize(trk);
    pParticleChange->ProposeTrackStatus(fStopAndKill);
    return pParticleChange;
  }

}
