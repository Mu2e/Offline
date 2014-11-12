#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"

#include "WLSEventAction.hh"
#include "WLSSteppingAction.hh"

#include "CrvPEresponse.hh"

#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

#include <TH3D.h>
#include <TNtuple.h>

WLSSteppingAction* WLSSteppingAction::fgInstance = NULL;

WLSSteppingAction::WLSSteppingAction(int mode, const std::string &lookupFileName) : _crvPEresponse(NULL), _mode(mode)
{
  fgInstance = this;

  if(_mode==0)
  {
    _crvPEresponse = new CrvPEresponse();
    _crvPEresponse->LoadLookupTable(lookupFileName);
  }
}

WLSSteppingAction::~WLSSteppingAction()
{
  if(_crvPEresponse) delete _crvPEresponse;
  _crvPEresponse=NULL;
}

void WLSSteppingAction::UserSteppingAction(const G4Step* theStep)
{
  G4OpBoundaryProcessStatus theStatus = Undefined;

  G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  if (OpManager) 
  {
     G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
     G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

     for ( G4int i=0; i<MAXofPostStepLoops; i++) 
     {
         G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
         G4OpBoundaryProcess *opProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
         if (opProcess) { theStatus = opProcess->GetStatus(); break;}
     }
  }

  if(theStatus==Detection)
  {
     G4VPhysicalVolume* thePostPV = theStep->GetPostStepPoint()->GetPhysicalVolume();
     G4String thePostPVname = " ";
     if(thePostPV)
     {
         if(thePostPV->GetName()=="PhotonDet")
         {
           PEs[0][thePostPV->GetCopyNo()]++;
           ArrivalTimes[0][thePostPV->GetCopyNo()].push_back(theStep->GetPostStepPoint()->GetGlobalTime());
           if(_mode==-1)
           {
             int numberOfFiberEmissions=0;
             int trackID = theStep->GetTrack()->GetTrackID();
             while(1)
             {
                std::map<int,int>::const_iterator wlsIter=_wlsTracks.find(trackID);
                if(wlsIter!=_wlsTracks.end())
                {
                  numberOfFiberEmissions++;
                  trackID=wlsIter->second;  //parent ID
                }
                else break;
             }
             FiberEmissions[thePostPV->GetCopyNo()].push_back(numberOfFiberEmissions);
           }
         }
     }
  }

  if(_mode==-1)
  {
     if(theStep->GetTrack()->GetCreatorProcess()!=NULL)
     {
       if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()=="OpWLS")
       {
         int trackID=theStep->GetTrack()->GetTrackID();
         int parentID=theStep->GetTrack()->GetParentID();
         _wlsTracks[trackID]=parentID;
       }
     }
  }

  if(_mode==0)
  {
    const G4ThreeVector &p1 = theStep->GetPreStepPoint()->GetPosition();
    const G4ThreeVector &p2 = theStep->GetPostStepPoint()->GetPosition();
    const double t1 = theStep->GetPreStepPoint()->GetGlobalTime();
    const double t2 = theStep->GetPostStepPoint()->GetGlobalTime();
    int PDGcode = theStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
    double beta = (theStep->GetPreStepPoint()->GetBeta() + theStep->GetPostStepPoint()->GetBeta())/2.0;
    double charge = theStep->GetTrack()->GetParticleDefinition()->GetPDGCharge();
    double energyDepositedTotal= theStep->GetTotalEnergyDeposit();
    double energyDepositedNonIonizing = theStep->GetNonIonizingEnergyDeposit();

    static bool first=true;
    if(first)
    {
      first=false;
      _crvPEresponse->SetScintillationYield(_scintillationYield);
      _crvPEresponse->SetScintillatorDecayTimeFast(_scintillatorDecayTimeFast);
      _crvPEresponse->SetScintillatorDecayTimeSlow(_scintillatorDecayTimeSlow);
      _crvPEresponse->SetFiberDecayTime(_fiberDecayTime);
    }

    if(PDGcode!=0)  //ignore optical photons
    {
      _crvPEresponse->MakePEs(p1, p2, t1, t2,  
                            PDGcode, beta, charge,
                            energyDepositedTotal,
                            energyDepositedNonIonizing);
 
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        PEs[1][SiPM]=_crvPEresponse->GetPEs(SiPM);
        std::vector<double> times=_crvPEresponse->GetArrivalTimes(SiPM);
        ArrivalTimes[1][SiPM].assign(times.begin(),times.end());
      }
    }
  }
}

int WLSSteppingAction::GetPEs(int i, int SiPM)
{
  return PEs[i][SiPM];
}

const std::vector<double> &WLSSteppingAction::GetArrivalTimes(int i, int SiPM)
{
  return ArrivalTimes[i][SiPM];
}

const std::vector<int> &WLSSteppingAction::GetFiberEmissions(int SiPM)
{
  return FiberEmissions[SiPM];
}

void WLSSteppingAction::Reset()
{
  for(int i=0; i<2; i++)
  for(int SiPM=0; SiPM<4; SiPM++)
  {
    PEs[i][SiPM]=0;
    ArrivalTimes[i][SiPM].clear();
    if(i==0) FiberEmissions[SiPM].clear();
  }

  _wlsTracks.clear();

  if(_mode==0)
  {
    _crvPEresponse->Reset();
  }
}

