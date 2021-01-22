#ifndef WIN32
#include <unistd.h>
#endif

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "WLSMaterials.hh"
#include "WLSPhysicsList.hh"
#include "WLSDetectorConstruction.hh"
#include "WLSPrimaryGeneratorAction.hh"

#include "WLSRunAction.hh"
#include "WLSEventAction.hh"
#include "WLSSteppingAction.hh"
#include "WLSStackingAction.hh"

#include "G4StepLimiterPhysics.hh"
#include "G4TransportationManager.hh"
#include "G4GDMLParser.hh"

#include "CLHEP/Random/Randomize.h"
#include "MakeCrvPhotons.hh"

#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>

bool findArgs(int argc, char** argv, const char* c)
{
  for(int i=1; i<argc; i++)
  {
    if(strcmp(argv[i],c)==0) return true;
  }
  return false; 
}

bool findArgs(int argc, char** argv, const char* c, std::string &value)
{
  for(int i=1; i<argc; i++)
  {
    if(strcmp(argv[i],c)==0)
    {
      if(i+1<argc)
      {
        value=argv[i+1];
        return true;
      }
      return false; 
    }
  }
  return false; 
}
bool findArgs(int argc, char** argv, const char* c, int &value)
{
  for(int i=1; i<argc; i++)
  {
    if(strcmp(argv[i],c)==0)
    {
      if(i+1<argc)
      {
        value=atoi(argv[i+1]);
        return true;
      }
      return false; 
    }
  }
  return false; 
}
bool findArgs(int argc, char** argv, const char* c, double &value)
{
  for(int i=1; i<argc; i++)
  {
    if(strcmp(argv[i],c)==0)
    {
      if(i+1<argc)
      {
        value=atof(argv[i+1]);
        return true;
      }
      return false; 
    }
  }
  return false; 
}

void DrawHistograms(const std::string &lookupFilename)
{
  std::unique_ptr<mu2eCrv::MakeCrvPhotons> crvPhotons;
  CLHEP::HepJamesRandom  engine(0);
  CLHEP::RandFlat randFlat(engine);
  CLHEP::RandGaussQ randGausQ(engine);
  CLHEP::RandPoissonQ randPoissonQ(engine);
  crvPhotons = std::unique_ptr<mu2eCrv::MakeCrvPhotons>(new mu2eCrv::MakeCrvPhotons(randFlat, randGausQ, randPoissonQ));
  crvPhotons->LoadLookupTable(lookupFilename);
  crvPhotons->DrawHistograms();
}

int main(int argc, char** argv) 
{
  WLSSteppingAction::simulationMode mode=WLSSteppingAction::Undefined;
  int simType=-1;
  int minBin=0;
  int maxBin=-1;
  int n=1000;
  double lengthOption=0;
  int    reflectorOption=0;
  std::string lookupFilename="";

  bool verbose = findArgs(argc, argv, "-v");

  if(findArgs(argc, argv, "-h"))
  {
    std::cout<<"Usage ./wls [OPTION]"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Available options:"<<std::endl;
    std::cout<<"-v           Verbose"<<std::endl;
    std::cout<<"-h           Help"<<std::endl;
    std::cout<<"-c           Create lookup table"<<std::endl;
    std::cout<<"-s           Run a simulation with full GEANT"<<std::endl;
    std::cout<<"-S           Run a simulation with lookup table"<<std::endl;
    std::cout<<"-d           Draw histograms of lookup table"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Options for creating the lookup table:"<<std::endl;
    std::cout<<"-t simtype        Simulation type:"<<std::endl;
    std::cout<<"                  0  scintillation in scintillator"<<std::endl;
    std::cout<<"                  1  cerenkov in scintillator"<<std::endl;
    std::cout<<"                  2  cerenkov in fiber"<<std::endl;
    std::cout<<"-l length option  Length of the scintillator counter in mm"<<std::endl;
    std::cout<<"-R reflector option"<<std::endl;
    std::cout<<"                  0  no reflector (default)"<<std::endl;
    std::cout<<"                 -1  reflector at negative side"<<std::endl;
    std::cout<<"                  1  reflector at postiive side"<<std::endl;
    std::cout<<"                 -2  black tape at negative side"<<std::endl;
    std::cout<<"                  2  black tape at positive side"<<std::endl;
    std::cout<<"-m minbin         Minimum bin in lookup table (default is 0)."<<std::endl;
    std::cout<<"-M maxbin         Maximum bin in lookup table (default is"<<std::endl;
    std::cout<<"                  the maximum number of bins for this simulation type)."<<std::endl;
    std::cout<<"-n photons        Number of photons to simulate for each bin (default 1000)."<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Options for running the simulation:"<<std::endl;
    std::cout<<"-l length option  Length of the scintillator counter in mm"<<std::endl;
    std::cout<<"-R reflector option"<<std::endl;
    std::cout<<"                  0  no reflector (default)"<<std::endl;
    std::cout<<"                 -1  reflector at negative side"<<std::endl;
    std::cout<<"                  1  reflector at positive side"<<std::endl;
    std::cout<<"                 -2  black tape at negative side"<<std::endl;
    std::cout<<"                  2  black tape at positive side"<<std::endl;
    std::cout<<"-n events    Number of events to simulate (default 1000)."<<std::endl;
    std::cout<<"-r seed      seed for random number generator (default: 0)."<<std::endl;
    std::cout<<"-y pos       y coordinate of starting point in mm (default: 0 = center between fibers)."<<std::endl;
    std::cout<<"-z pos       z coordinate of starting point in mm (default: 1000 = 1m away from left side of counter)."<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Options for running the simulation with lookup table:"<<std::endl;
    std::cout<<"-f filename  File with lookup table used for running a simulation"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Options for drawing histograms of the lookup table:"<<std::endl;
    std::cout<<"-f filename  File with lookup table used for running a simulation"<<std::endl;
    std::cout<<std::endl;
    return 0;
  }

  if(findArgs(argc, argv, "-c")+findArgs(argc, argv, "-s")+findArgs(argc, argv, "-S")+findArgs(argc, argv, "-d")!=1)
  {
    std::cout<<"Need to specify exactly one of the following options: -c, -s, -S, or -d."<<std::endl;
    std::cout<<"Use -h for help."<<std::endl;
    return -1;
  }

  if(findArgs(argc, argv, "-S") || findArgs(argc, argv, "-d"))
  {
    if(!findArgs(argc, argv, "-f", lookupFilename))
    {
      std::cout<<"Filename for lookup table needs to be specified"<<std::endl;
      std::cout<<"Use -h for help"<<std::endl;
      return -1;
    }
  }

  if(findArgs(argc, argv, "-d"))
  {
    DrawHistograms(lookupFilename);
    return 0;
  }

  if(findArgs(argc, argv, "-c")) mode=WLSSteppingAction::CreateLookupTables;
  if(findArgs(argc, argv, "-s")) mode=WLSSteppingAction::UseGeantOnly;
  if(findArgs(argc, argv, "-S")) mode=WLSSteppingAction::UseGeantAndLookupTables;

  if(mode==WLSSteppingAction::CreateLookupTables)
  {
    if(!findArgs(argc, argv, "-t", simType))
    {
      std::cout<<"Simulation type needs to be specified"<<std::endl;
      std::cout<<"Use -h for help"<<std::endl;
      return -1;
    }
    findArgs(argc, argv, "-m", minBin);
    findArgs(argc, argv, "-M", maxBin);
    findArgs(argc, argv, "-n", n);
  }

  if(mode==WLSSteppingAction::UseGeantOnly || mode==WLSSteppingAction::UseGeantAndLookupTables)
  {
    findArgs(argc, argv, "-n", n);
  }

  if(!findArgs(argc, argv, "-l", lengthOption))
  {
    std::cout<<"Option for scintillator counter length needs to be specified"<<std::endl;
    std::cout<<"Use -h for help"<<std::endl;
    return -1;
  }
  findArgs(argc, argv, "-R", reflectorOption);

  double posY=0;
  double posZ=1000;
  int r=0;
  findArgs(argc, argv, "-y", posY);
  findArgs(argc, argv, "-z", posZ);
  findArgs(argc, argv, "-r", r);

  G4int seed = r;

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(seed);

  G4RunManager *runManager = new G4RunManager;

  WLSMaterials::GetInstance();
  runManager->SetUserInitialization(new WLSDetectorConstruction(lengthOption, reflectorOption));
  runManager->SetUserInitialization(new WLSPhysicsList());

  WLSPrimaryGeneratorAction *generator = new WLSPrimaryGeneratorAction(mode, n, simType, minBin, verbose, posY, posZ);   
                                                                       //n,simType,minBin not needed in modes UseGeantOnly, and UseGeantAndLookupTables
                                                                       //posY, posZ has no effect in mode CreateLookupTables
  WLSRunAction* runAction = new WLSRunAction();
  std::string singlePEWaveformFilename="singlePEWaveform_v3.txt";
  std::string photonMapFilename="photonMap.root";
  WLSEventAction* eventAction = new WLSEventAction(mode, singlePEWaveformFilename, photonMapFilename, n, simType, minBin, verbose); 
  WLSSteppingAction* steppingAction = new WLSSteppingAction(mode, lookupFilename); //lookupFilename not needed in modes CreateLookupTables, and UseGeantOnly
  WLSStackingAction* stackingAction = new WLSStackingAction();

  runManager->SetUserAction(generator);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(steppingAction);
  runManager->SetUserAction(stackingAction);
  runManager->Initialize();

/*
  G4GDMLParser parser;
  parser.SetRegionExport(true);
  parser.Write("out.gdml", G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume()->GetLogicalVolume());
*/

  if(mode==WLSSteppingAction::CreateLookupTables)
  {
    const std::vector<double> &xBins     = WLSDetectorConstruction::Instance()->GetXBins();
    const std::vector<double> &yBins     = WLSDetectorConstruction::Instance()->GetYBins();
    const std::vector<double> &zBins     = WLSDetectorConstruction::Instance()->GetZBins();
    const std::vector<double> &betaBins  = WLSDetectorConstruction::Instance()->GetBetaBins();
    const std::vector<double> &thetaBins = WLSDetectorConstruction::Instance()->GetThetaBins();
    const std::vector<double> &phiBins   = WLSDetectorConstruction::Instance()->GetPhiBins();
    const std::vector<double> &rBins     = WLSDetectorConstruction::Instance()->GetRBins();

    int nXBins=xBins.size()-1;   //e.g. 3 bins need 4 entries in the vector (for 3 bin boundaries)
    int nYBins=yBins.size()-1;
    int nZBins=zBins.size()-1;
    int nBetaBins=betaBins.size()-1;
    int nThetaBins=thetaBins.size()-1;
    int nPhiBins=phiBins.size()-1;
    int nRBins=rBins.size()-1;

    int nBins=0;
    switch(simType)
    {
      case 0: //scintillation in scintillator
              nBins=nZBins*nYBins*nXBins;
              break;
      case 1: //Cerenkov in scintillator
              nBins=nBetaBins*nZBins*nYBins*nXBins;
              break;
      case 2: //Cerenkov in fiber
              nBins=nZBins*nRBins*nPhiBins*nThetaBins*nBetaBins;
              break;
    }

    if(maxBin==-1 || maxBin>=nBins) maxBin=nBins-1;

    int numberOfEvents=maxBin-minBin+1;
    if(numberOfEvents>0)
    {
      std::cout<<std::endl<<std::endl;
      std::cout<<"About to simulate "<<numberOfEvents<<" bins from "<<minBin<<" to "<<maxBin<<"."<<std::endl;
      runManager->BeamOn(numberOfEvents);
      std::cout<<std::endl<<std::endl;
    }
    else
    {

      std::cout<<std::endl<<std::endl;
      std::cout<<"Exceeded the maximum number of bins. no simulation done."<<std::endl;
      std::cout<<std::endl<<std::endl;

      std::stringstream filename;
      filename<<"LookupTable_"<<simType<<"_";
      filename.fill('0');
      filename.width(6);
      filename<<minBin;
      std::remove(filename.str().c_str());
      std::ofstream lookupfile(filename.str(),std::ios::binary|std::ios::app);  //create an empty file to avoid anomalies in the script
      lookupfile.close();
    }
  }
  else if(mode==WLSSteppingAction::UseGeantOnly || mode==WLSSteppingAction::UseGeantAndLookupTables)
  {
    runManager->BeamOn(n);
  }

  delete runManager;

  return 0;
}
