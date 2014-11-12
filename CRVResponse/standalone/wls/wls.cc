#ifndef WIN32
#include <unistd.h>
#endif

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "WLSPhysicsList.hh"
#include "WLSDetectorConstruction.hh"
#include "WLSPrimaryGeneratorAction.hh"

#include "WLSRunAction.hh"
#include "WLSEventAction.hh"
#include "WLSSteppingAction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

// argc holds the number of arguments (including the name) on the command line
// -> it is ONE when only the name is  given !!!
// argv[0] is always the name of the program
// argv[1] points to the first argument, and so on

int main(int argc,char** argv) 
{
  if(argc<2) return(-1);
  int mode=atoi(argv[1]);
  if(mode==-1 && argc!=6) return(-1);
  if(mode==0 && argc!=3) return(-1);

  G4String physName = "QGSP_BERT_EMV";
//  G4String physName = "QGSP_BERT_HP";  //for neutrons
  G4int seed = 123;

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(seed);

  G4RunManager *runManager = new G4RunManager;
  WLSDetectorConstruction* detector = new WLSDetectorConstruction();

  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new WLSPhysicsList(physName));

  WLSPrimaryGeneratorAction *generator = new WLSPrimaryGeneratorAction(mode);
  WLSRunAction* runAction = new WLSRunAction(mode);
  WLSEventAction* eventAction = new WLSEventAction(mode, atoi(argv[5]));

  if(mode==-1)
  {
    int binx=atoi(argv[2]);
    int biny=atoi(argv[3]);
    int binz=atoi(argv[4]);
    generator->SetBins(binx,biny,binz);
  }

  runManager->SetUserAction(generator);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);

  if(mode==-1) runManager->SetUserAction( new WLSSteppingAction(mode) );
  else runManager->SetUserAction( new WLSSteppingAction(mode, argv[2]) );

#if 0 
  G4VisManager *visManager = new G4VisExecutive();
  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  G4UIExecutive *ui = new G4UIExecutive(1,argv);

  visManager->Initialize();
  UImanager->ApplyCommand("/run/initialize");
  UImanager->ApplyCommand("/control/execute vis.mac");
  UImanager->ApplyCommand("/run/beamOn");

  ui->SessionStart();

  delete ui;
  delete visManager;

#else

  runManager->Initialize();
  if(mode==-1)
  {
    int binx=atoi(argv[2]);
    int biny=atoi(argv[3]);
    int binz=atoi(argv[4]);

    int numberOfEvents=2; //one event for scintillation, and one for cerenkov
    if(binx==1 && biny==1 && binz==1) WLSEventAction::Instance()->doStoreConstants(true);
    if(binx==1 && biny==1) numberOfEvents=4; //two additional event for Cerenkov in the fibers

    runManager->BeamOn(numberOfEvents);
  }
  else if(mode==0)
  {
    runManager->BeamOn(1000);
  }
#endif

  delete runManager;

  return 0;
}
