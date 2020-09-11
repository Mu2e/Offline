//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// Original author KLG
//
// Modeled after G4VisCommandSceneHandler...
//
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined  G4VIS_USE_OPENGLQT ) 
#include "Mu2eG4/inc/Mu2eVisCommands.hh"

#include "G4VVisCommand.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

Mu2eVisCommandSceneHandlerDrawEvent::Mu2eVisCommandSceneHandlerDrawEvent() {
  
  fpCommand = new G4UIcommand("/vis/sceneHandler/drawEvent", this);
  fpCommand -> SetGuidance("Draws current event at the current scene/viewer");
  fpCommand -> SetGuidance("If a viewer does not exists, will try to create\n"
			    "one for the current sceneHandler");

  G4UIparameter* parameter;
  G4bool omitable;
  parameter = new G4UIparameter("viewer-name", 's', omitable = true);
  parameter -> SetCurrentAsDefault(true);
  fpCommand -> SetParameter(parameter);

}

Mu2eVisCommandSceneHandlerDrawEvent::~Mu2eVisCommandSceneHandlerDrawEvent() {
  delete fpCommand;
}

G4String Mu2eVisCommandSceneHandlerDrawEvent::GetCurrentValue(G4UIcommand*) {
  G4VSceneHandler* pSceneHandler = fpVisManager->GetCurrentSceneHandler();
  const G4ViewerList& viewerList = pSceneHandler->GetViewerList();
  G4int nViewers = viewerList.size();
  if (nViewers > 0) {
    G4VViewer* pViewer = pSceneHandler->GetCurrentViewer();
    return pViewer ? pViewer->GetName() : G4String("");
  } else {
    return G4String("");
  }
}

void Mu2eVisCommandSceneHandlerDrawEvent::SetNewValue(G4UIcommand*,
						  G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<
      "ERROR: No valid scenes available.  Please create one."
	     << G4endl;
    }
    return;
  }

  G4VSceneHandler* pSceneHandler = fpVisManager->GetCurrentSceneHandler();
  if (!pSceneHandler) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<
      "ERROR: Current scene handler not defined.  Please select or create one."
	     << G4endl;
    }
    return;
  }

  G4VViewer* pViewer;

  // pSceneHandler->GetViewCount() and pSceneHandler->GetCurrentViewer()
  // do not return expected results, pSceneHandler->GetViewerList() does
  const G4ViewerList& viewerList = pSceneHandler->GetViewerList();
  G4int nViewers = viewerList.size();
  if (nViewers == 0 || !(pSceneHandler->GetCurrentViewer()) ) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
      "WARNING: Current viewer not defined or no viewers.  Will try to create one"
	     << G4endl;
    }
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/create ! "+newValue);
    pViewer = pSceneHandler->GetCurrentViewer();
    if (!pViewer) {    
      if (verbosity >= G4VisManager::errors) {
	G4cout <<
	  "ERROR: Could not create a viewer, returning"
	       << G4endl;
      }
      return;
    }
  } 

  const G4Event* currentEvent = G4RunManager::GetRunManager()->GetCurrentEvent();

  if (!currentEvent) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<
	"ERROR: No current event, returning"
	     << G4endl;
    }
    return;
  }

  // fpVisManager->ClearTransientStoreIfMarked(); // this "clears the old tracks"
  // it is private and it looks like we may not need it

  pSceneHandler->DrawEvent(currentEvent);

  // const G4ViewParameters& viewParams = pViewer->GetViewParameters();
  // if (viewParams.IsAutoRefresh()) {
  //   G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
  // }
  // else {
  //   if (verbosity >= G4VisManager::warnings) {
  //     G4cout << "Issue /vis/viewer/refresh or flush to see the effect."
  // 	     << G4endl;
  //   }
  // }

  return;

}
#endif
