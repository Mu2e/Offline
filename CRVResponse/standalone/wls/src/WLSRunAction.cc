#include "WLSRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"

#include "WLSDetectorConstruction.hh"
#include "WLSSteppingAction.hh"
#include "WLSEventAction.hh"

#include <algorithm>
#include <ctime>

WLSRunAction::WLSRunAction()
{
}

WLSRunAction::~WLSRunAction()
{
}

void WLSRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");

  CLHEP::HepRandom::showEngineStatus();
}

void WLSRunAction::EndOfRunAction(const G4Run* )
{
}
