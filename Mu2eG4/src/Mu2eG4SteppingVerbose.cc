//
// Verbose version of the stepping action.
//
//
// Original author Rob Kutschke
//
// The intial release is just a copy of the G4 Novice N02 example.
//

#include "Offline/Mu2eG4/inc/Mu2eG4SteppingVerbose.hh"

#include "Geant4/G4SteppingManager.hh"
#include "Geant4/G4UnitsTable.hh"

namespace mu2e {

  Mu2eG4SteppingVerbose::Mu2eG4SteppingVerbose(){}

  Mu2eG4SteppingVerbose::~Mu2eG4SteppingVerbose(){}

  void Mu2eG4SteppingVerbose::StepInfo(){
    CopyState();

    G4int prec = G4cout.precision(3);

    if( verboseLevel >= 1 ){
      if( verboseLevel >= 4 ) VerboseTrack();
      if( verboseLevel >= 3 ){
        G4cout << G4endl;
        G4cout << std::setw( 5) << "#Step#"     << " "
               << std::setw( 6) << "X"          << "    "
               << std::setw( 6) << "Y"          << "    "
               << std::setw( 6) << "Z"          << "    "
               << std::setw( 9) << "KineE"      << " "
               << std::setw( 9) << "dEStep"     << " "
               << std::setw(10) << "StepLeng"
               << std::setw(10) << "TrakLeng"
               << std::setw(10) << "Volume"    << "  "
               << std::setw(10) << "Process"   << G4endl;
      }

      G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
             << std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
             << std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
             << std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
             << std::setw(6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
             << std::setw(6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
             << std::setw(6) << G4BestUnit(fStep->GetStepLength(),"Length")
             << std::setw(6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
             << "  ";

      // if( fStepStatus != fWorldBoundary){
      if( fTrack->GetNextVolume() != 0 ) {
        G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
      } else {
        G4cout << std::setw(10) << "OutOfWorld";
      }

      if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != nullptr){
        G4cout << "  "
               << std::setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
          ->GetProcessName();
      } else {
        G4cout << "   UserLimit";
      }

      G4cout << G4endl;

      if( verboseLevel == 2 ){
        G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
          fN2ndariesAlongStepDoIt +
          fN2ndariesPostStepDoIt;
        if(tN2ndariesTot>0){
          G4cout << "    :----- List of 2ndaries - "
                 << "#SpawnInStep=" << std::setw(3) << tN2ndariesTot
                 << "(Rest="  << std::setw(2) << fN2ndariesAtRestDoIt
                 << ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
                 << ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
                 << "), "
                 << "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
                 << " ---------------"
                 << G4endl;

          for(size_t lp1=(*fSecondary).size()-tN2ndariesTot;
              lp1<(*fSecondary).size(); lp1++){
            G4cout << "    : "
                   << std::setw(6)
                   << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
                   << std::setw(6)
                   << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
                   << std::setw(6)
                   << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
                   << std::setw(6)
                   << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
                   << std::setw(10)
                   << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
            G4cout << G4endl;
          }

          G4cout << "    :-----------------------------"
                 << "----------------------------------"
                 << "-- EndOf2ndaries Info ---------------"
                 << G4endl;
        }
      }

    }
    G4cout.precision(prec);
  }


  void Mu2eG4SteppingVerbose::TrackingStarted()
  {

    CopyState();
    G4int prec = G4cout.precision(3);
    if( verboseLevel > 0 ){

      G4cout << std::setw( 5) << "Step#"      << " "
             << std::setw( 6) << "X"          << "    "
             << std::setw( 6) << "Y"          << "    "
             << std::setw( 6) << "Z"          << "    "
             << std::setw( 9) << "KineE"      << " "
             << std::setw( 9) << "dEStep"     << " "
             << std::setw(10) << "StepLeng"
             << std::setw(10) << "TrakLeng"
             << std::setw(10) << "Volume"     << "  "
             << std::setw(10) << "Process"    << G4endl;

      G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
             << std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
             << std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
             << std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
             << std::setw(6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
             << std::setw(6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
             << std::setw(6) << G4BestUnit(fStep->GetStepLength(),"Length")
             << std::setw(6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
             << "  ";

      if(fTrack->GetNextVolume()){
        G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
      } else {
        G4cout << std::setw(10) << "OutOfWorld";
      }
      G4cout  << "    initStep" << G4endl;
    }
    G4cout.precision(prec);
  }

} // end namespace mu2e
